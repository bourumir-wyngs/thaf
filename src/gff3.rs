use crate::error::Error;
use crate::structures::{Strand, TranscriptRegion};
use std::collections::{HashMap, HashSet};
use std::fs::File;
use std::io::Write;
use std::io::{BufRead, BufReader, BufWriter};

pub fn parse_gff3_to_regions(
    gff3_path: &str,
    feature_types: &Vec<String>,
    errors: &mut Vec<Error>,
) -> anyhow::Result<Vec<TranscriptRegion>> {
    let feature_set: HashSet<&str> = feature_types.iter().map(|s| s.as_str()).collect();
    let reader = BufReader::new(File::open(gff3_path)?);
    let mut regions = Vec::new();

    let mut transcript_to_gene: HashMap<String, String> = HashMap::new();
    let mut warn_missing_tx_parent = false;
    let mut warn_missing_feature_parent = false;

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
        let strand_char = cols[6].chars().next().unwrap_or('.');

        let attributes = parse_attributes(cols[8]);

        match feature_type {
            "gene" => {
                if let Some(gene_id) = attributes.get("ID") {
                    transcript_to_gene.insert(gene_id.clone(), gene_id.clone());
                }
            }
            "mRNA" | "transcript" => {
                if let Some(transcript_id) = attributes.get("ID") {
                    let gene_id = if let Some(parent) = attributes.get("Parent") {
                        parent.clone()
                    } else {
                        if !warn_missing_tx_parent {
                            errors.push(Error::warning(
                                "Transcript entry missing Parent attribute; using transcript ID as gene ID",
                            ));
                            warn_missing_tx_parent = true;
                        }
                        transcript_id.clone()
                    };
                    transcript_to_gene.insert(transcript_id.clone(), gene_id);
                }
            }
            feat if feature_set.contains(feat) => {
                let region_id = if let Some(id) = attributes.get("ID") {
                    id.clone()
                } else {
                    errors.push(Error::fatal("Missing transcript id".to_string()));
                    continue;
                };

                let transcript_id = if let Some(parent) = attributes.get("Parent") {
                    parent.clone()
                } else {
                    if !warn_missing_feature_parent {
                        errors.push(Error::warning(
                            "Feature missing Parent attribute; using feature ID as transcript and gene ID",
                        ));
                        warn_missing_feature_parent = true;
                    }
                    transcript_to_gene.insert(region_id.clone(), region_id.clone());
                    region_id.clone()
                };

                let gene_id = transcript_to_gene.get(&transcript_id).cloned();

                if let Some(strand) = Strand::from_char(strand_char, errors) {
                    regions.push(TranscriptRegion {
                        chromosome,
                        start,
                        end,
                        region_id,
                        strand,
                        transcript_id,
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

#[cfg(test)]
mod tests {
    use super::*;
    use crate::error::Severity;
    #[test]
    fn test_parse_attributes_basic() {
        let attrs = parse_attributes("ID=exon1;Parent=tx1;");
        assert_eq!(attrs.get("ID"), Some(&"exon1".to_string()));
        assert_eq!(attrs.get("Parent"), Some(&"tx1".to_string()));
    }

    #[test]
    fn test_parse_gff3_to_regions_simple() {
        use std::io::Write;
        let mut file = tempfile::NamedTempFile::new().unwrap();
        writeln!(file, "chr1	src	gene	1	10	.	+	.	ID=gene1;").unwrap();
        writeln!(file, "chr1	src	mRNA	1	10	.	+	.	ID=tx1;Parent=gene1;").unwrap();
        writeln!(file, "chr1	src	exon	1	5	.	+	.	ID=ex1;Parent=tx1;").unwrap();
        writeln!(file, "chr1	src	exon	6	10	.	+	.	ID=ex2;Parent=tx1;").unwrap();
        let path = file.path().to_str().unwrap().to_string();
        let mut errors = Vec::new();
        let regions = parse_gff3_to_regions(&path, &vec!["exon".to_string()], &mut errors).unwrap();
        assert!(errors.is_empty());
        assert_eq!(regions.len(), 2);
        assert_eq!(regions[0].transcript_id, "tx1");
    }

    #[test]
    fn test_missing_parents() {
        use std::io::Write;
        let mut file = tempfile::NamedTempFile::new().unwrap();
        writeln!(file, "chr1\tsrc\tmRNA\t1\t10\t.\t+\t.\tID=tx1").unwrap();
        writeln!(file, "chr1\tsrc\texon\t1\t5\t.\t+\t.\tID=ex1;Parent=tx1").unwrap();
        writeln!(file, "chr1\tsrc\texon\t6\t10\t.\t+\t.\tID=ex2").unwrap();
        let path = file.path().to_str().unwrap().to_string();
        let mut errors = Vec::new();
        let regions = parse_gff3_to_regions(&path, &vec!["exon".to_string()], &mut errors).unwrap();
        assert_eq!(regions.len(), 2);
        assert!(
            errors
                .iter()
                .any(|e| matches!(e.severity, Severity::Warning))
        );
        assert_eq!(regions[0].transcript_id, "tx1");
        assert_eq!(regions[1].transcript_id, "ex2");
    }
}
