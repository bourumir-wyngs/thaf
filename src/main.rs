mod transcript_builder;
mod structures;
mod gff3;
mod error;

use crate::gff3::{parse_gff3_to_regions, write_genemap};
use crate::transcript_builder::{build_transcriptome_sequences, build_transcripts_from_regions};
use crate::error::{Error, Severity};
use anyhow::Result;
use clap::{Arg, Command};

fn main() -> Result<()> {
    let matches = Command::new("thaf")
        .version("0.0.2")
        .about("Extract transcript sequences from genome FASTA based on GFF3 annotations.")
        .arg(
            Arg::new("gff3")
                .short('f')
                .long("gff3")
                .value_name("INPUT_GFF3")
                .help("Input GFF3 annotation file")
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
        .arg(
            Arg::new("dna")
                .short('d')
                .long("dna")
                .value_name("DNA_FASTA")
                .help("Genome FASTA file for extracting sequences")
                .required(true),
        )
        .arg(
            Arg::new("transcriptome")
                .long("transcriptome")
                .short('t')
                .value_name("OUTPUT_FASTA")
                .help("Output FASTA file for transcript sequences")
                .required(true),
        )
        .arg(
            Arg::new("features")
                .short('e')
                .long("features")
                .value_name("FEATURES")
                .help("Features to extract (comma-separated, defaults to 'exon')")
                .required(false),
        )
        .arg(
            Arg::new("error")
                .short('r')
                .long("error")
                .value_name("ERROR_LOG")
                .help("Write warnings and errors to this file")
                .required(false),
        )
        .arg_required_else_help(true)
        .get_matches();

    let input_file = matches.get_one::<String>("gff3").unwrap();
    let dna_fasta = matches.get_one::<String>("dna").unwrap();
    let transcriptome_fasta = matches.get_one::<String>("transcriptome").unwrap();
    let genemap_file = matches.get_one::<String>("genemap");
    let error_file = matches.get_one::<String>("error");
    let features: Vec<String> = matches
        .get_one::<String>("features")
        .map(|s| s.split(',').map(|item| item.trim().to_string()).collect())
        .unwrap_or_else(|| vec!["exon".to_string()]);

    let mut errors: Vec<Error> = Vec::new();
    
    println!("  Features: {:?}", features);

    // Parsing regions from GFF3
    let regions = parse_gff3_to_regions(input_file,
                                        &features,
                                        &mut errors)?;

    // Optionally write genemap
    if let Some(genemap_path) = genemap_file {
        write_genemap(&regions, genemap_path)?;
    }

    // Build transcripts from regions
    let transcripts = build_transcripts_from_regions(regions, &mut errors);

    // Extract and write transcript sequences
    build_transcriptome_sequences(&transcripts, dna_fasta, transcriptome_fasta)?;

    if let Some(path) = error_file {
        use std::fs::File;
        use std::io::Write;
        let mut f = File::create(path)?;
        for e in &errors {
            writeln!(f, "[{:?}] {}", e.severity, e.message)?;
        }
    } else {
        for e in &errors {
            println!("[{:?}] {}", e.severity, e.message);
        }
    }

    if errors.iter().any(|e| matches!(e.severity, Severity::Fatal)) {
        std::process::exit(1);
    }

    Ok(())
}


