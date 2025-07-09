pub(crate) use crate::structures::{Region, Strand, Transcript, TranscriptRegion};
use anyhow::Result;
use bio::alphabets::dna;
use bio::data_structures::interval_tree::IntervalTree;
use bio::io::fasta;
use std::collections::HashMap;
use log::error;

impl Transcript {
    pub fn new(id: String, chromosome: String, mut regions: Vec<Region>) -> Result<Self> {
        if regions.is_empty() {
            anyhow::bail!("Transcript {} has no regions.", id);
        }

        // Verify all regions have the same strand
        let first_strand = regions[0].strand;
        if regions.iter().any(|r| r.strand != first_strand) {
            anyhow::bail!("Transcript {} has mixed strands.", id);
        }

        // Sort regions depending on the strand
        match first_strand {
            Strand::Plus => regions.sort_by_key(|r| r.start),
            Strand::Minus => regions.sort_by(|a, b| b.start.cmp(&a.start)),
        }

        // Check for overlapping regions
        let mut interval_tree: IntervalTree<usize, &Region> = IntervalTree::new();

        for region in &regions {
            if region.start > region.end {
                panic!("Negative width region {}..{}, region {} strand {}", region.start, region.end, region.id, region.strand);
            } else if region.end - region.start + 1 <= 3 {
                println!("  Suspicious: {} is only {} nucleoptide length: {} .. {}"
                    , region.id, region.end - region.start + 1, region.start, region.end);
            }
            let interval = region.start..region.end + 1; // bio uses half-open intervals
            if let Some(overlap) = interval_tree.find(interval.clone()).next() {
                panic!(
                    "Transcript {} in chromosome {} has overlapping regions: {} and {} overlap with interval {:?}.",
                    id,
                    chromosome,
                    region.id,
                    overlap.data().id,
                    overlap.interval()
                );
            }
            interval_tree.insert(interval, region);
        }

        Ok(Self {
            id,
            chromosome,
            regions,
        })
    }
}

/// Build transcripts from a vector of TranscriptRegion structs.
pub fn build_transcripts_from_regions(
    transcript_regions: Vec<TranscriptRegion>,
) -> Result<Vec<Transcript>> {
    // Collect regions grouped by transcript ID
    let mut transcript_map: HashMap<String, (String, Vec<Region>)> = HashMap::new();

    for tr in transcript_regions {
        let entry = transcript_map
            .entry(tr.transcript_id.clone())
            .or_insert_with(|| (tr.chromosome.clone(), Vec::new()));

        // Sanity-check chromosome consistency
        if entry.0 != tr.chromosome {
            panic!(
                "Transcript {} has regions from multiple chromosomes: {} vs {}",
                tr.transcript_id,
                entry.0,
                tr.chromosome
            );
        }

        entry.1.push(Region {
            id: tr.region_id.clone(),
            start: tr.start,
            end: tr.end,
            strand: tr.strand,
        });
    }

    // Now build validated transcripts
    let mut transcripts = Vec::new();

    for (id, (chromosome, regions)) in transcript_map {
        let transcript = Transcript::new(id, chromosome, regions)?;
        transcripts.push(transcript);
    }

    Ok(transcripts)
}

/// Load genome sequences into memory from FASTA file.
fn load_genome_to_memory(fasta_path: &str) -> Result<HashMap<String, Vec<u8>>> {
    let reader = fasta::Reader::from_file(fasta_path)?;
    let mut genome = HashMap::new();

    for record in reader.records() {
        let record = record?;
        genome.insert(record.id().to_owned(), record.seq().to_owned());
    }
    Ok(genome)
}

/// Extract sequence for a single transcript.
fn extract_transcript_sequence(
    genome: &HashMap<String, Vec<u8>>,
    transcript: &Transcript,
) -> Result<Vec<u8>> {
    let chromosome_seq = genome.get(&transcript.chromosome).ok_or_else(|| {
        anyhow::anyhow!(
            "Chromosome '{}' not found in genome.",
            transcript.chromosome
        )
    })?;

    let mut sequence = Vec::with_capacity(transcript.size());

    // Sort regions according to strand orientation
    let mut sorted_regions = transcript.regions.clone();
    sorted_regions.sort_by_key(|r| r.start);

    // Extract sequences:
    for region in &sorted_regions {
        let start = region.start - 1;
        let end = region.end;

        sequence.extend_from_slice(&chromosome_seq[start..end]);
    }

    // Reverse complement entire sequence for minus strand:
    if transcript.regions[0].strand == Strand::Minus {
        sequence = dna::revcomp(sequence);
    }

    Ok(sequence)
}

fn _extract_transcript_sequence(
    genome: &HashMap<String, Vec<u8>>,
    transcript: &Transcript,
) -> Result<Vec<u8>> {
    let chromosome_seq = genome.get(&transcript.chromosome).ok_or_else(|| {
        anyhow::anyhow!(
            "Chromosome '{}' not found in genome.",
            transcript.chromosome
        )
    })?;

    let mut sorted_regions = transcript.regions.clone();

    match transcript.regions[0].strand {
        Strand::Plus => sorted_regions.sort_by_key(|r| r.start),
        Strand::Minus => sorted_regions.sort_by(|a, b| b.start.cmp(&a.start)),
    }

    let mut sequence = Vec::new();

    // Extract sequences with debug annotations
    for (index, region) in sorted_regions.iter().enumerate() {
        let start = region.start - 1; // 0-based indexing
        let end = region.end;         // exclusive at end

        if end > chromosome_seq.len() {
            anyhow::bail!(
                "Region ({}-{}) in transcript {} exceeds chromosome length.",
                region.start,
                region.end,
                transcript.id
            );
        }

        // Insert debug start marker
        let marker = match region.strand {
            Strand::Plus => format!("[{}+:", index + 1),
            Strand::Minus => format!("[{}-:", index + 1),
        };
        sequence.extend_from_slice(marker.as_bytes());

        // Actual exon sequence
        sequence.extend_from_slice(&chromosome_seq[start..end]);

        // Insert debug end marker
        let marker_end = match region.strand {
            Strand::Plus => format!("{}+]", index + 1),
            Strand::Minus => format!("{}-]", index + 1),
        };
        sequence.extend_from_slice(marker_end.as_bytes());
    }

    // Reverse complement if needed
    if transcript.regions[0].strand == Strand::Minus {
        sequence = dna::revcomp(sequence);
    }

    Ok(sequence)
}

/// Build transcriptome sequences and write to FASTA file.
pub fn build_transcriptome_sequences(
    transcripts: &[Transcript],
    genome_fasta_path: &str,
    output_fasta_path: &str,
) -> Result<()> {
    // Load genome into memory
    let genome = load_genome_to_memory(genome_fasta_path)?;

    // Open FASTA writer for output
    let mut writer = fasta::Writer::to_file(output_fasta_path)?;

    // Extract and write each transcript
    for transcript in transcripts {
        let seq = extract_transcript_sequence(&genome, transcript)?;
        writer.write(&transcript.id, None, &seq)?;
    }

    Ok(())
}
