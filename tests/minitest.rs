use std::fs::File;
use std::io::Write;
use tempfile::tempdir;
use bio::io::fasta;
use thaf::gff3::parse_gff3_to_regions;
use thaf::transcript_builder::{build_transcriptome_sequences, build_transcripts_from_regions};
use thaf::error::Error;
use thaf::error::Severity;

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
        writeln!(gff, "Super-Scaffold_197 SLAB mRNA 213511889 213514721 . + . ID=STRG.36498.1.p3;\
Parent=STRG.36498;Name=ORF type:complete len:528 (-)%2Cscore%3D112.45%2Ctr|YYY|YYY|98.712|\
0.0%2A|PF0YY93.17|8.5e-138%2CACP_syn_III_C|YYY|1.4e+04%2CACP_syn_III_C|PF08541.14|1.6e+04%2\
CACP_syn_III_C|PF08541.14|4.9e+03%2CACP_syn_III_C|PF08541.14|5.9e-12%2CChal_sti_synt_C\
|PF02797.19|1.7e-10%2CChal_sti_synt_N|PF00195.23|0.00047%2CACP_syn_III|PF08545.14|0.00078")?;
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
    let mut errors = Vec::<Error>::new();
    let regions = parse_gff3_to_regions(gff3_path.to_str().unwrap(), &vec!["exon".into()], &mut errors)?;
    let transcripts = build_transcripts_from_regions(regions, &mut errors);
    assert_eq!(errors.len(), 4);
    assert!(errors.iter().all(|e| matches!(e.severity, Severity::Warning)));

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