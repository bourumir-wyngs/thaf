mod transcript_builder;
mod structures;
mod gff3;

use crate::gff3::{parse_gff3_to_regions, write_genemap};
use crate::transcript_builder::{build_transcriptome_sequences, build_transcripts_from_regions};
use anyhow::Result;
use clap::{Arg, Command};

fn feature_synonyms() -> Vec<String> {
    // We could support various other synonyms here as well or introduce a CLI parameter.
    // Will be seen.
    let v = vec!["exon"];
    //let v = vec!["CDS"];
    v.iter().map(|s| s.to_string()).collect()
}

fn main() -> Result<()> {
    let matches = Command::new("gff3_parser")
        .version("1.4")
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
        .arg_required_else_help(true)
        .get_matches();

    let input_file = matches.get_one::<String>("gff3").unwrap();
    let dna_fasta = matches.get_one::<String>("dna").unwrap();
    let transcriptome_fasta = matches.get_one::<String>("transcriptome").unwrap();
    let genemap_file = matches.get_one::<String>("genemap");
    let features: Vec<String> = matches
        .get_one::<String>("features")
        .map(|s| s.split(',').map(|item| item.trim().to_string()).collect())
        .unwrap_or_else(|| vec!["exon".to_string()]);
    
    println!("  Features: {:?}", features);

    // Parsing regions from GFF3
    let regions = parse_gff3_to_regions(input_file, 
                                        &features)?;

    // Optionally write genemap
    if let Some(genemap_path) = genemap_file {
        write_genemap(&regions, genemap_path)?;
    }

    // Build transcripts from regions
    let transcripts = build_transcripts_from_regions(regions)?;

    // Extract and write transcript sequences
    build_transcriptome_sequences(&transcripts, dna_fasta, transcriptome_fasta)?;

    Ok(())
}


