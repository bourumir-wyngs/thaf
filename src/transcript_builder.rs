pub(crate) use crate::structures::{Region, Strand, Transcript, TranscriptRegion};
use anyhow::Result;
use bio::alphabets::dna;
use bio::data_structures::interval_tree::IntervalTree;
use bio::io::fasta;
use std::collections::HashMap;
use crate::error::Error;

impl Transcript {
    pub fn new(id: String, chromosome: String, mut regions: Vec<Region>, errors: &mut Vec<Error>) -> Option<Self> {
        if regions.is_empty() {
            errors.push(Error::fatal(format!("Transcript {} has no regions.", id)));
            return None;
        }

        // Verify all regions have the same strand
        let first_strand = regions[0].strand;
        if regions.iter().any(|r| r.strand != first_strand) {
            errors.push(Error::fatal(format!("Transcript {} has mixed strands.", id)));
            return None;
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
                errors.push(Error::fatal(format!("Negative width region {}..{}, region {} strand {}", region.start, region.end, region.id, region.strand)));
                return None;
            } else if region.end - region.start + 1 <= 3 {
                errors.push(Error::warning(format!("Suspicious: {} is only {} nucleoptide length: {} .. {}", region.id, region.end - region.start + 1, region.start, region.end)));
            }
            let interval = region.start..region.end + 1; // bio uses half-open intervals
            if let Some(overlap) = interval_tree.find(interval.clone()).next() {
                errors.push(Error::fatal(format!(
                    "Transcript {} in chromosome {} has overlapping regions: {} and {} overlap with interval {:?}.",
                    id,
                    chromosome,
                    region.id,
                    overlap.data().id,
                    overlap.interval()
                )));
                return None;
            }
            interval_tree.insert(interval, region);
        }

        Some(Self {
            id,
            chromosome,
            regions,
        })
    }
}

/// Build transcripts from a vector of TranscriptRegion structs.
pub fn build_transcripts_from_regions(
    transcript_regions: Vec<TranscriptRegion>,
    errors: &mut Vec<Error>,
) -> Vec<Transcript> {
    // Collect regions grouped by transcript ID
    let mut transcript_map: HashMap<String, (String, Vec<Region>)> = HashMap::new();

    for tr in transcript_regions {
        let entry = transcript_map
            .entry(tr.transcript_id.clone())
            .or_insert_with(|| (tr.chromosome.clone(), Vec::new()));

        // Sanity-check chromosome consistency
        if entry.0 != tr.chromosome {
            errors.push(Error::fatal(format!(
                "Transcript {} has regions from multiple chromosomes: {} vs {}",
                tr.transcript_id,
                entry.0,
                tr.chromosome
            )));
            continue;
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
        if let Some(transcript) = Transcript::new(id, chromosome, regions, errors) {
            transcripts.push(transcript);
        }
    }

    transcripts
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

#[cfg(test)]
mod tests {
    use super::*;
    use crate::error::Severity;

    fn build_region(id: &str, start: usize, end: usize, strand: Strand) -> Region {
        Region { id: id.to_string(), start, end, strand }
    }

    #[test]
    fn test_transcript_sort_plus() {
        let r1 = build_region("r1", 20, 25, Strand::Plus);
        let r2 = build_region("r2", 10, 15, Strand::Plus);
        let mut errors = Vec::new();
        let t = Transcript::new("tx1".into(), "chr1".into(), vec![r1, r2], &mut errors).unwrap();
        assert!(errors.is_empty());
        assert!(t.regions[0].start < t.regions[1].start);
    }

    #[test]
    fn test_transcript_sort_minus() {
        let r1 = build_region("r1", 10, 15, Strand::Minus);
        let r2 = build_region("r2", 20, 25, Strand::Minus);
        let mut errors = Vec::new();
        let t = Transcript::new("tx1".into(), "chr1".into(), vec![r1, r2], &mut errors).unwrap();
        assert!(errors.is_empty());
        assert!(t.regions[0].start > t.regions[1].start);
    }

    #[test]
    fn test_transcript_overlap_fails() {
        let r1 = build_region("r1", 10, 20, Strand::Plus);
        let r2 = build_region("r2", 15, 25, Strand::Plus);
        let mut errors = Vec::new();
        let res = Transcript::new("tx1".into(), "chr1".into(), vec![r1, r2], &mut errors);
        assert!(res.is_none());
        assert!(!errors.is_empty());
    }

    #[test]
    fn test_extract_transcript_sequence_plus() {
        let genome = HashMap::from([("chr1".to_string(), b"ACGTAACCGGTT".to_vec())]);
        let regions = vec![build_region("r1", 1, 4, Strand::Plus), build_region("r2", 5, 8, Strand::Plus)];
        let mut errors = Vec::new();
        let t = Transcript::new("tx1".into(), "chr1".into(), regions, &mut errors).unwrap();
        assert!(errors.is_empty());
        let seq = extract_transcript_sequence(&genome, &t).unwrap();
        assert_eq!(seq, b"ACGTAACC");
    }

    #[test]
    fn test_extract_transcript_sequence_minus() {
        let genome = HashMap::from([("chr1".to_string(), b"ACGTAACCGGTT".to_vec())]);
        let regions = vec![build_region("r1", 1, 4, Strand::Minus), build_region("r2", 5, 8, Strand::Minus)];
        let mut errors = Vec::new();
        let t = Transcript::new("tx1".into(), "chr1".into(), regions, &mut errors).unwrap();
        assert!(errors.is_empty());
        let seq = extract_transcript_sequence(&genome, &t).unwrap();
        assert_eq!(seq, b"GGTTACGT");
    }

    #[test]
    fn test_build_transcripts_from_regions() {
        let trs = vec![TranscriptRegion { chromosome: "chr1".into(), start: 1, end: 3, strand: Strand::Plus, transcript_id: "tx1".into(), region_id: "r1".into(), gene_id: None },
                        TranscriptRegion { chromosome: "chr1".into(), start: 5, end: 6, strand: Strand::Plus, transcript_id: "tx1".into(), region_id: "r2".into(), gene_id: None }];
        let mut errors = Vec::new();
        let ts = build_transcripts_from_regions(trs, &mut errors);
        assert_eq!(errors.len(), 2);
        assert!(errors.iter().all(|e| matches!(e.severity, Severity::Warning)));
        assert_eq!(ts.len(), 1);
        assert_eq!(ts[0].regions.len(), 2);
    }
}
