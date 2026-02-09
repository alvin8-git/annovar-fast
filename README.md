# ANNOVAR-fast

High-performance variant annotation using Python with tabix-indexed databases. Drop-in replacement for the ANNOVAR `convert2annovar.pl` + `table_annovar.pl` pipeline.

## Overview

ANNOVAR-fast replaces the original Perl-based ANNOVAR workflow with a Python implementation that uses `pysam` for O(1) tabix lookups and `cyvcf2` for fast VCF parsing. It produces output identical to the standard ANNOVAR `table_annovar.pl` multianno format (170 columns).

## Requirements

- Python 3.7+
- [cyvcf2](https://github.com/brentp/cyvcf2) - fast VCF parsing
- [pysam](https://github.com/pysam-developers/pysam) - tabix database access
- Tabix-indexed ANNOVAR databases in `humandb-tbi/`
- Gene mRNA FASTA files (`hg38_refGeneMrna.fa`, `hg38_ensGeneMrna.fa`) in `humandb/`

```bash
pip install cyvcf2 pysam
```

## Usage

```bash
python3 annovar-fast.py input.vcf [-o output.txt] [-v]
```

- `input.vcf` - Input VCF file
- `-o` - Output file (default: `<sample>.annovar.txt`)
- `-v` - Verbose/debug logging

If a reference file `<sample>.annovar.txt` exists alongside the input VCF, the tool automatically compares output against it.

## Architecture

```
annovar-fast.py          Main entry point and CLI
config.py                Database paths, column mappings, protocol order
vcf_processor.py         VCF parsing and ANNOVAR format conversion (cyvcf2)
database_manager.py      Tabix queries for filter and region databases (pysam)
gene_annotator.py        refGene/ensGene annotation with AA change computation
annotation_engine.py     Orchestrates all annotators, writes output
```

## Databases

Supports 39 annotation databases matching the `annovar-hg38.awk` protocol:

| Type | Databases |
|------|-----------|
| Gene | refGene, ensGene |
| Filter | esp6500siv2, exac03, 1000g2015aug (6 pops), intervar, kaviar, mcap, revel, gnomad211_exome, gnomad30_genome, snp141, ljb26_all, clinvar, dbscsnv11, hrcr1, gme, icgc28, TCGA, civic, tumorportal, cg69, nci60 |
| Region | cytoBand, genomicSuperDups, cosmic, gwasCatalog, encRegTfbsClustered, wgEncodeRegDnaseClustered, wgRna, dgvMerged, phastConsElements30way, phastConsElements100way, targetScanS |

## Configuration

Database paths are set in `config.py`:

- `HUMANDB_TBI_DIR` - directory containing `*.txt.gz` + `*.txt.gz.tbi` tabix files
- `HUMANDB_DIR` - directory containing gene mRNA FASTA files

Edit these paths to match your installation.

## Output Format

Tab-separated with 170 columns: 5 variant fields (Chr, Start, End, Ref, Alt), annotation columns from all databases in protocol order, and 11 Otherinfo columns (zygosity + original VCF fields).
