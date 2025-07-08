use anyhow::Result;
use clap::{Arg, Command};
use server::pipeline::tasks::transcripts::transcript_builder::{build_transcripts_from_regions, Strand, TranscriptRegion};
use std::{
    collections::HashMap,
    fs::File,
    io::{BufRead, BufReader, BufWriter, Write},
};

fn feature_synonyms() -> Vec<&'static str> {
    if false {
        vec![
            "CDS", "five_prime_UTR", "three_prime_UTR",
            "5UTR", "3UTR", "UTR5", "UTR3", "5'UTR", "3'UTR", "UTR5'", "UTR3'",
            "five_prime_untranslated_region", "three_prime_untranslated_region",
        ]
    } else {
        vec!["exon"]
    }
}

fn main() -> Result<()> {
    let matches = Command::new("gff3_parser")
        .version("1.3")
        .about("Extract transcript regions and gene maps from GFF3.")
        .arg(
            Arg::new("input")
                .short('i')
                .long("input")
                .value_name("INPUT_GFF3")
                .help("Input GFF3 annotation file")
                .required(true),
        )
        .arg(
            Arg::new("output")
                .short('o')
                .long("output")
                .value_name("OUTPUT_FILE")
                .help("Output TSV file for transcript regions")
                .required(true),
        )
        .arg(
            Arg::new("genemap")
                .short('g')
                .long("genemap")
                .value_name("GENEMAP_FILE")
                .help("Output TSV file for transcript-to-gene mapping")
                .required(false),
        )
        .arg_required_else_help(true)
        .get_matches();

    let input_file = matches.get_one::<String>("input").unwrap();
    let output_file = matches.get_one::<String>("output").unwrap();
    let genemap_file = matches.get_one::<String>("genemap");

    let regions = parse_gff3_to_regions(input_file)?;
    write_regions_to_tsv(&regions, output_file)?;

    if let Some(genemap_path) = genemap_file {
        write_genemap(&regions, genemap_path)?;
    }

    let transcripts = build_transcripts_from_regions(regions)?;    

    Ok(())
}

pub fn parse_gff3_to_regions(gff3_path: &str) -> Result<Vec<TranscriptRegion>> {
    let reader = BufReader::new(File::open(gff3_path)?);
    let feature_types = feature_synonyms();
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
            feat if feature_types.contains(&feat) => {
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

pub fn write_regions_to_tsv(regions: &[TranscriptRegion], out_path: &str) -> Result<()> {
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

pub fn write_genemap(regions: &[TranscriptRegion], out_path: &str) -> Result<()> {
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

