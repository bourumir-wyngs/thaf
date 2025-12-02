<p>
    <a href="https://verdanta.info#oss"><img src="https://verdanta.tech/verdanta_logo_small.jpg"
    alt="Relative median time vs baseline"
    width="100px"/>
    </a>
</p>

[![GitHub Workflow Status](https://img.shields.io/github/actions/workflow/status/bourumir-wyngs/thaf/rust.yml)](https://github.com/bourumir-wyngs/thaf/actions)
[![crates.io](https://img.shields.io/crates/v/thaf.svg)](https://crates.io/crates/thaf)
[![crates.io](https://img.shields.io/crates/l/thaf.svg)](https://crates.io/crates/thaf)
[![crates.io](https://img.shields.io/crates/d/thaf.svg)](https://crates.io/crates/thaf)
# GFF3 Parser and Transcriptome Extractor

## Overview

`thaf` is a command-line tool to extract transcript sequences from a genome FASTA file based on GFF3 annotation files. It can also generate transcript-to-gene mapping files compatible with tools such as Salmon.

## Features

* Parses GFF3 annotation files to identify transcript regions.
* The default feature to be extracted is 'exon', but this is easy to change with the -e switch
* Extracts transcript sequences directly from genome FASTA files.
* Handles forward and reverse strands automatically.
* Generates transcript-to-gene mapping files.

## Usage

### Command-Line Arguments

```bash
thaf \
  -f <INPUT_GFF3> \
  -d <DNA_FASTA> \
  -t <OUTPUT_FASTA> \
  [-g <GENEMAP_FILE>]
  [-e <FEATURES>]
```

### Required Arguments

* `-f, --gff3 <INPUT_GFF3>`: Path to the input GFF3 annotation file.
* `-d, --dna <DNA_FASTA>`: Path to the input genome FASTA file.
* `-t, --transcriptome <OUTPUT_FASTA>`: Path to the output transcriptome FASTA file.

### Optional Arguments

* `-g, --genemap <GENEMAP_FILE>`: Path to the output TSV file for transcript-to-gene mapping.
* `-e, --features <FEATURES>`: Comma-separated list of GFF3 features to extract (default: exon).
* `-r, --error <ERROR_LOG>`: Write warnings and errors to this file instead of standard output.

## Example

```bash
thaf \
  -f annotations.gff3 \
  -d genome.fa \
  -t transcriptome.fa \
  -g genemap.tsv
  -r run.log \
  -e CDS
```

This will produce:

* `transcriptome.fa`: FASTA file containing extracted transcript sequences.
* `genemap.tsv`: Tab-separated file mapping transcripts to genes.

This small project was inspired by a segmentation fault encountered while using one of the popular tools, and the lack of any readily available tool capable of producing even a simple `genemap` table.

Sequence boundaries, exon order, and reverse-complementation have been validated against outputs from `gffread`, which unfortunately does not produce a `genemap`. `thaf` checks for obvious inconsistencies, such as overlapping exons or exons belonging to different strands or chromosomes.

Unlike `gffread`, `thaf` loads the entire genome into memory. As a result, it cannot handle extremely large genomes, such as that of the fern *Tmesipteris oblanceolata* (~160 Gb). However, a typical 32 Gb workstation is enough for processing the crop and plant genomes we commonly work with, and the simpler algorithm should make the code easier to maintain.

We are grateful to the [**rust-bio**](https://crates.io/crates/bio) package, which provides exon overlap detection and reverse-complement functionality.

