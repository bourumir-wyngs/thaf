use std::fs::File;
use std::io::Write;
use tempfile::tempdir;
use bio::io::fasta;
use thaf::gff3::parse_gff3_to_regions;
use thaf::transcript_builder::{build_transcriptome_sequences, build_transcripts_from_regions};

#[test]
fn minimal_transcript_extraction() -> anyhow::Result<()> {
    // Create temporary directory for test files
    let dir = tempdir()?;
    let gff3_path = dir.path().join("test.gff3");
    let genome_path = dir.path().join("genome.fa");
    let transcriptome_path = dir.path().join("trans.fa");

    // Prepare minimal GFF3 with five exon features on two chromosomes
    {
        let mut gff = File::create(&gff3_path)?;
        writeln!(gff, "##gff-version 3")?;
        writeln!(gff, "chr1\tsrc\tgene\t1\t20\t.\t+\t.\tID=g1")?;
        writeln!(gff, "chr1\tsrc\tmRNA\t1\t20\t.\t+\t.\tID=tx1;Parent=g1")?;
        writeln!(gff, "chr1\tsrc\texon\t1\t3\t.\t+\t.\tID=ex1;Parent=tx1")?;
        writeln!(gff, "chr1\tsrc\texon\t5\t8\t.\t+\t.\tID=ex2;Parent=tx1")?;
        writeln!(gff, "chr1\tsrc\tgene\t9\t20\t.\t-\t.\tID=g2")?;
        writeln!(gff, "chr1\tsrc\tmRNA\t9\t20\t.\t-\t.\tID=tx2;Parent=g2")?;
        writeln!(gff, "chr1\tsrc\texon\t14\t16\t.\t-\t.\tID=ex3;Parent=tx2")?;
        writeln!(gff, "chr1\tsrc\texon\t18\t20\t.\t-\t.\tID=ex4;Parent=tx2")?;
        writeln!(gff, "chr2\tsrc\tgene\t1\t10\t.\t+\t.\tID=g3")?;
        writeln!(gff, "chr2\tsrc\tmRNA\t1\t10\t.\t+\t.\tID=tx3;Parent=g3")?;
        writeln!(gff, "chr2\tsrc\texon\t2\t4\t.\t+\t.\tID=ex5;Parent=tx3")?;
    }

    // Minimal genome with two chromosomes
    {
        let mut genome = File::create(&genome_path)?;
        writeln!(genome, ">chr1")?;
        writeln!(genome, "AAACCCGGGTTTAAACCCGG")?;
        writeln!(genome, ">chr2")?;
        writeln!(genome, "GGTTAACCAA")?;
    }

    // Parse regions and build transcripts
    let regions = parse_gff3_to_regions(gff3_path.to_str().unwrap(), &vec!["exon".into()])?;
    let transcripts = build_transcripts_from_regions(regions)?;

    // Build transcript sequences
    build_transcriptome_sequences(&transcripts, genome_path.to_str().unwrap(), transcriptome_path.to_str().unwrap())?;

    // Read produced FASTA and collect sequences
    let reader = fasta::Reader::from_file(transcriptome_path)?;
    let mut seqs = std::collections::HashMap::new();
    for rec in reader.records() {
        let rec = rec?;
        seqs.insert(rec.id().to_owned(), String::from_utf8(rec.seq().to_vec())?);
    }

    assert_eq!(seqs.get("tx1").unwrap(), "AAACCGG");
    assert_eq!(seqs.get("tx2").unwrap(), "CCGGTT");
    assert_eq!(seqs.get("tx3").unwrap(), "GTT");

    Ok(())
}