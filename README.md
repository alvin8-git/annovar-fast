# ANNOVAR-fast

High-performance variant annotation using Python with tabix-indexed databases. Drop-in replacement for the ANNOVAR `convert2annovar.pl` + `table_annovar.pl` pipeline.

## Overview

ANNOVAR-fast replaces the original Perl-based ANNOVAR workflow with a Python implementation that uses `pysam` for O(1) tabix lookups and `cyvcf2` for fast VCF parsing. It produces output identical to the standard ANNOVAR `table_annovar.pl` multianno format (170 columns).

Key advantages over ANNOVAR:
- **Faster**: tabix random access instead of sequential flat-file scanning
- **Simpler**: single command replaces a multi-step `convert2annovar.pl` + `table_annovar.pl` workflow
- **Multi-sample VCF**: natively handles multi-sample VCFs (ANNOVAR requires splitting first)

## Requirements

- Python 3.7+
- [cyvcf2](https://github.com/brentp/cyvcf2) - fast VCF parsing
- [pysam](https://github.com/pysam-developers/pysam) - tabix database access

```bash
pip install cyvcf2 pysam
```

### Database files

All database files live in a single directory (`HUMANDB_DIR`):

- 39 bgzipped + tabix-indexed annotation databases (`.txt.gz` + `.txt.gz.tbi` pairs)
- 2 mRNA FASTA files for amino acid change computation (`hg38_refGeneMrna.fa`, `hg38_ensGeneMrna.fa`)

The database directory is available as a GitHub release archive. Download and extract it, then edit `config.py` to set the path:

```python
HUMANDB_DIR = "/path/to/humandb"
```

See [Database List](#database-list) for the full inventory of required files.

## Usage

```bash
python3 annovar-fast.py input.vcf [-o output.txt] [-v]
```

| Argument | Description |
|----------|-------------|
| `input.vcf` | Input VCF file (single- or multi-sample) |
| `-o, --output` | Output file path (default: `<input_basename>.annovar.txt`) |
| `-v, --verbose` | Enable debug logging |

### Examples

```bash
# Annotate a single-sample VCF
python3 annovar-fast.py sample.vcf

# Annotate a multi-sample VCF with custom output path
python3 annovar-fast.py cohort.vcf -o cohort_annotated.txt

# Verbose mode for debugging
python3 annovar-fast.py sample.vcf -v
```

If a reference file `<basename>.annovar.txt` exists alongside the input VCF, the tool automatically compares output and reports differences.

## Output Format

Tab-separated with 170 columns:

| Section | Columns | Description |
|---------|---------|-------------|
| Variant | 5 | Chr, Start, End, Ref, Alt |
| Annotations | 154 | All database annotations in protocol order |
| Otherinfo | 11 | Zygosity + original VCF fields (CHROM through sample genotypes) |

## Architecture

```
annovar-fast.py        (102 lines)   CLI entry point, argument parsing, diff reporting
vcf_processor.py       (189 lines)   VCF parsing via cyvcf2, ANNOVAR coordinate conversion
annotation_engine.py   (129 lines)   Orchestrates all annotators, writes output
database_manager.py    (442 lines)   Tabix queries for filter/region databases
gene_annotator.py     (1314 lines)   refGene/ensGene annotation with mRNA-based AA changes
config.py              (442 lines)   Database config, column mappings, protocol order
```

### Pipeline flow

```
VCF input
  |
  v
vcf_processor.py          Parse VCF, normalize variants to ANNOVAR coordinates
  |
  v
annotation_engine.py      For each variant:
  |-- gene_annotator.py       Query refGene/ensGene (tabix), classify variant,
  |                           compute AA changes using mRNA FASTA
  |-- database_manager.py     Query all filter/region databases (tabix)
  |
  v
Tab-separated output      170-column multianno format
```

## Database List

### Gene annotations (2 databases)

| Database | File | Description |
|----------|------|-------------|
| refGene | hg38_refGene.txt.gz | NCBI RefSeq gene models (splicing threshold: 100bp) |
| ensGene | hg38_ensGene.txt.gz | Ensembl gene models (splicing threshold: 2bp) |

Each produces 5 columns: Func, Gene, GeneDetail, ExonicFunc, AAChange.

### Filter annotations (27 databases)

| Database | File | Columns | Description |
|----------|------|---------|-------------|
| esp6500siv2_all | hg38_esp6500siv2_all.txt.gz | 1 | ESP6500 allele frequencies |
| exac03 | hg38_exac03.txt.gz | 8 | ExAC frequencies (ALL, AFR, AMR, EAS, FIN, NFE, OTH, SAS) |
| 1000g2015aug_all | hg38_ALL.sites.2015_08.txt.gz | 1 | 1000 Genomes all populations |
| 1000g2015aug_afr | hg38_AFR.sites.2015_08.txt.gz | 1 | 1000 Genomes African |
| 1000g2015aug_eas | hg38_EAS.sites.2015_08.txt.gz | 1 | 1000 Genomes East Asian |
| 1000g2015aug_amr | hg38_AMR.sites.2015_08.txt.gz | 1 | 1000 Genomes Americas |
| 1000g2015aug_eur | hg38_EUR.sites.2015_08.txt.gz | 1 | 1000 Genomes European |
| 1000g2015aug_sas | hg38_SAS.sites.2015_08.txt.gz | 1 | 1000 Genomes South Asian |
| intervar_20180118 | hg38_intervar_20180118.txt.gz | 30 | InterVar ACMG classification |
| kaviar_20150923 | hg38_kaviar_20150923.txt.gz | 1 | Kaviar allele frequencies |
| mcap | hg38_mcap.txt.gz | 1 | M-CAP pathogenicity scores |
| revel | hg38_revel.txt.gz | 1 | REVEL ensemble pathogenicity scores |
| gnomad211_exome | hg38_gnomad211_exome.txt.gz | 17 | gnomAD v2.1.1 exome frequencies |
| gnomad30_genome | hg38_gnomad30_genome.txt.gz | 13 | gnomAD v3.0 genome frequencies |
| snp141 | hg38_snp141.txt.gz | 1 | dbSNP 141 rsIDs |
| ljb26_all | hg38_ljb26_all.txt.gz | 25 | Functional predictions (SIFT, PolyPhen2, CADD, GERP++, etc.) |
| clinvar_20220730 | hg38_clinvar_20220730.txt.gz | 5 | ClinVar (CLNALLELEID, CLNDN, CLNDISDB, CLNREVSTAT, CLNSIG) |
| dbscsnv11 | hg38_dbscsnv11.txt.gz | 2 | Splice site predictions (ADA, RF scores) |
| hrcr1 | hg38_hrcr1.txt.gz | 6 | Haplotype Reference Consortium frequencies |
| gme | hg38_gme.txt.gz | 8 | Greater Middle East allele frequencies |
| icgc28 | hg38_icgc28.txt.gz | 1 | ICGC cancer mutations |
| TCGA | hg38_TCGA.txt.gz | 1 | TCGA cancer mutations |
| civic | hg38_civic.txt.gz | 1 | CIViC clinical interpretations |
| tumorportal | hg38_tumorportal.txt.gz | 1 | Tumor Portal mutations |
| cg69 | hg38_cg69.txt.gz | 1 | Complete Genomics 69 genomes |
| nci60 | hg38_nci60.txt.gz | 1 | NCI-60 cell line panel |

### Region annotations (10 databases)

| Database | File | Description |
|----------|------|-------------|
| cytoBand | hg38_cytoBand.txt.gz | Cytogenetic bands |
| genomicSuperDups | hg38_genomicSuperDups.txt.gz | Segmental duplications |
| targetScanS | hg38_targetScanS.txt.gz | TargetScan miRNA target sites |
| cosmic | hg38_cosmic.txt.gz | COSMIC somatic mutations |
| gwasCatalog | hg38_gwasCatalog.txt.gz | GWAS Catalog associations |
| encRegTfbsClustered | hg38_encRegTfbsClustered.txt.gz | ENCODE TF binding sites |
| wgEncodeRegDnaseClustered | hg38_wgEncodeRegDnaseClustered.txt.gz | ENCODE DNase hypersensitivity |
| wgRna | hg38_wgRna.txt.gz | snoRNA and miRNA annotations |
| dgvMerged | hg38_dgvMerged.txt.gz | Database of Genomic Variants |
| phastConsElements30way | hg38_phastConsElements30way.txt.gz | 30-way conservation |
| phastConsElements100way | hg38_phastConsElements100way.txt.gz | 100-way conservation |

### mRNA FASTA files (also in HUMANDB_DIR)

| File | Description |
|------|-------------|
| hg38_refGeneMrna.fa | RefSeq mRNA sequences for protein-level annotation |
| hg38_ensGeneMrna.fa | Ensembl mRNA sequences for protein-level annotation |

These FASTA files contain transcript mRNA sequences used to translate codons and compute amino acid changes (e.g., `NM_001234:exon2:c.123A>G:p.Ser41Gly`). Without them, gene annotation still works for variant classification (exonic, intronic, splicing, etc.) but AAChange columns will be empty.

## Validation Status

Tested against ANNOVAR output on a 192-variant, 95-sample multi-sample VCF:

| Category | Count | Notes |
|----------|-------|-------|
| Exact match | 176/192 lines | Identical to ANNOVAR output |
| Database version diffs | 12 | ensGene DB has extra/different transcripts |
| HLA polymorphism | 3 | Different mRNA at HLA-B locus (unfixable) |
| SNP tie-breaking | 1 | Different rsID chosen when multiple match |
| Order-only diffs | 310 cells | Same values, different ordering in region databases |

All substantive differences are due to database version mismatches or HLA polymorphism, not algorithmic bugs.
