use std::collections::{HashMap, HashSet};
use std::fs::File;
use std::io::{BufRead, BufReader, BufWriter};
use crate::structures::{Strand, TranscriptRegion};
use std::io::Write;

pub fn parse_gff3_to_regions(gff3_path: &str, feature_types: &Vec<String>) -> anyhow::Result<Vec<TranscriptRegion>> {
    let feature_set: HashSet<&str> = feature_types.iter().map(|s| s.as_str()).collect();
    let reader = BufReader::new(File::open(gff3_path)?);
    let mut regions = Vec::new();

    let mut transcript_to_gene: HashMap<String, String> = HashMap::new();

    for line in reader.lines() {
        let line = line?;
        if line.starts_with('#') || line.trim().is_empty() {
            continue;
        }

        let cols: Vec<&str> = line.split('\t').collect();
        if cols.len() != 9 {
            continue;
        }

        let chromosome = cols[0].to_owned();
        let feature_type = cols[2];
        let start = cols[3].parse::<usize>()?;
        let end = cols[4].parse::<usize>()?;
        let strand = cols[6].chars().next().unwrap_or('.');

        let attributes = parse_attributes(cols[8]);

        match feature_type {
            "gene" => {
                if let Some(gene_id) = attributes.get("ID") {
                    transcript_to_gene.insert(gene_id.clone(), gene_id.clone());
                }
            }
            "mRNA" | "transcript" => {
                if let (Some(transcript_id), Some(gene_id)) =
                    (attributes.get("ID"), attributes.get("Parent"))
                {
                    transcript_to_gene.insert(transcript_id.clone(), gene_id.clone());
                }
            }
            feat if feature_set.contains(feat) => {
                if let Some(transcript_id) = attributes.get("Parent") {
                    let gene_id = transcript_to_gene.get(transcript_id).cloned();
                    let region_id = if let Some(id) = attributes.get("ID") {
                        id.clone()
                    } else {
                        panic!("Missing transcript id");
                    };

                    regions.push(TranscriptRegion {
                        chromosome,
                        start,
                        end,
                        region_id,
                        strand: Strand::from_char(strand),
                        transcript_id: transcript_id.clone(),
                        gene_id,
                    });
                }
            }
            _ => (),
        }
    }

    Ok(regions)
}

#[allow(dead_code)]
pub fn write_regions_to_tsv(regions: &[TranscriptRegion], out_path: &str) -> anyhow::Result<()> {
    let mut writer = BufWriter::new(File::create(out_path)?);
    writeln!(writer, "chromosome\tstart\tend\tstrand\ttranscript_id")?;

    for region in regions {
        writeln!(
            writer,
            "{}\t{}\t{}\t{}\t{}",
            region.chromosome, region.start, region.end, region.strand, region.transcript_id
        )?;
    }

    Ok(())
}

pub fn write_genemap(regions: &[TranscriptRegion], out_path: &str) -> anyhow::Result<()> {
    let mut writer = BufWriter::new(File::create(out_path)?);
    let mut seen = HashMap::new();

    writeln!(writer, "transcript_id\tgene_id")?;

    for region in regions {
        if let Some(gene_id) = &region.gene_id {
            let transcript_id = &region.transcript_id;
            if !seen.contains_key(transcript_id) {
                writeln!(writer, "{}\t{}", transcript_id, gene_id)?;
                seen.insert(transcript_id, ());
            }
        }
    }

    Ok(())
}

fn parse_attributes(attr_str: &str) -> HashMap<String, String> {
    attr_str
        .trim_end_matches(';')
        .split(';')
        .filter_map(|attr| attr.split_once('='))
        .map(|(k, v)| (k.trim().to_owned(), v.trim().to_owned()))
        .collect()
}

