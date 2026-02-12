# annovar-fast

High-performance variant annotation for hg38 whole-genome and whole-exome sequencing data. Drop-in replacement for ANNOVAR's `convert2annovar.pl` + `table_annovar.pl` pipeline.

## Features

- **Fast** -- Fork-based parallel annotation across all CPU cores with tabix random access. WGS (4.7M variants) completes in ~14 minutes on 64 cores vs ~160 minutes with ANNOVAR.
- **Simple** -- Single command replaces a multi-step ANNOVAR workflow.
- **Multi-sample VCF** -- Natively handles multi-sample VCFs without splitting.
- **39 databases** -- Gene annotation (refGene, ensGene), population frequencies (gnomAD 4.1, 1000G, ExAC, ESP), pathogenicity scores (CADD, REVEL, SIFT, PolyPhen2), clinical (ClinVar 2025, InterVar 2025), cancer (COSMIC, ICGC, TCGA), conservation (phastCons), and more.

## Quick start

```bash
git clone https://github.com/alvin8-git/annovar-fast.git
cd annovar-fast
pip install cyvcf2 pysam
# Download and configure the database (see Database Setup below)
python3 annovar-fast.py input.vcf.gz
```

## Requirements

- Python 3.7+
- Linux (fork-based parallelism requires POSIX `fork()`)
- [cyvcf2](https://github.com/brentp/cyvcf2) -- Fast VCF parsing via htslib
- [pysam](https://github.com/pysam-developers/pysam) -- Tabix-indexed database access

```bash
pip install cyvcf2 pysam
```

## Installation

```bash
git clone https://github.com/alvin8-git/annovar-fast.git
cd annovar-fast
```

No compilation or `setup.py` needed. The tool runs directly from the cloned directory.

## Database setup

annovar-fast requires a set of pre-indexed annotation databases (29 GB total). All files must be in a single directory.

### 1. Download

Download the database archive:

```bash
# TODO: Replace with actual download URL
wget -O humandb-tbi.tar.gz <DOWNLOAD_URL>
```

### 2. Extract

```bash
tar xzf humandb-tbi.tar.gz -C /path/to/databases
```

This creates the database files (`.txt.gz` + `.txt.gz.tbi` pairs and `.fa` files) in the target directory.

### 3. Configure

Edit `config.py` and set the `HUMANDB_DIR` path to your extracted database directory:

```python
HUMANDB_DIR = "/path/to/databases"
```

### 4. Verify

Run the tool with `--verbose` to confirm all 39 databases load successfully:

```bash
python3 annovar-fast.py --help
python3 annovar-fast.py sample.vcf -v 2>&1 | grep "Loaded .* databases"
# Expected: "Loaded 37 databases" (37 filter/region + 2 gene handled separately)
```

## Usage

```bash
python3 annovar-fast.py <input.vcf> [-o output.txt] [-v]
```

| Argument | Description |
|----------|-------------|
| `input.vcf` | Input VCF file (`.vcf` or `.vcf.gz`, single- or multi-sample) |
| `-o, --output` | Output file path (default: `<input_stem>.annovar.txt`) |
| `-v, --verbose` | Enable debug logging |

### Examples

```bash
# Annotate a single-sample WES VCF
python3 annovar-fast.py sample.vcf.gz

# Custom output path
python3 annovar-fast.py cohort.vcf -o cohort_annotated.txt

# Debug mode
python3 annovar-fast.py sample.vcf -v
```

### Sequential mode

Parallel annotation is enabled by default. To force sequential processing (useful for debugging), call from Python:

```python
from annotation_engine import AnnotationEngine
engine = AnnotationEngine()
engine.process_vcf("input.vcf.gz", "output.txt", parallel=False)
engine.close()
```

## Output format

Tab-separated file with 170+ columns:

| Section | Columns | Description |
|---------|---------|-------------|
| Variant | 5 | Chr, Start, End, Ref, Alt |
| Annotations | 154 | All database annotations in protocol order |
| Otherinfo | 11+ | Zygosity + original VCF fields (CHROM through sample genotypes) |

The Otherinfo column count scales with the number of samples in the input VCF (11 for single-sample, 11 + N-1 for multi-sample).

## Performance

Benchmarks on a 64-core machine:

| Dataset | Variants | ANNOVAR (Perl) | annovar-fast | Speedup |
|---------|----------|---------------|--------------|---------|
| WES (single-sample) | 74,755 | ~10 min | 23 sec | ~26x |
| WGS (single-sample) | 4,699,828 | 160 min | 14 min | ~11x |
| Multi-sample (95 samples) | 192 | ~1 min | 2 sec | ~30x |

Memory usage is approximately 600 MB for the parent process (mRNA sequences) plus ~10 MB per worker (tabix handles only; mRNA data is shared via copy-on-write after fork).

## Architecture

```
annovar-fast.py          CLI entry point
vcf_processor.py         VCF parsing (cyvcf2), ANNOVAR coordinate conversion
annotation_engine.py     Orchestrator, parallel dispatch, output writing
database_manager.py      Tabix queries for filter/region databases
gene_annotator.py        refGene/ensGene annotation, AA change computation
config.py                Database config, column mappings, protocol order
```

### How parallelism works

```
Parent process:
  1. Load all databases + mRNA sequences (~600 MB)
  2. Read all variants from VCF
  3. Split into equal-size chunks (one per CPU core)
  4. Fork worker pool (workers inherit mRNA via copy-on-write)

Each worker:
  1. Re-open tabix file handles (file descriptors can't be shared across fork)
  2. Annotate its chunk of variants independently
  3. Return (index, formatted_line) tuples

Parent process:
  5. Merge results in original variant order
  6. Write output file
```

## Databases

### Gene annotations

| Database | File | Columns | Description |
|----------|------|---------|-------------|
| refGene | hg38_refGene.txt.gz | 5 | NCBI RefSeq gene models (splicing threshold: 100 bp) |
| ensGene | hg38_ensGene_new.txt.gz | 5 | Ensembl gene models (splicing threshold: 2 bp) |

Each produces: Func, Gene, GeneDetail, ExonicFunc, AAChange.

mRNA FASTA files (`hg38_refGeneMrna.fa`, `hg38_ensGeneMrna.fa`) provide transcript sequences for amino acid change computation (e.g., `p.Ser41Gly`).

### Population frequency databases

| Database | File | Columns | Description |
|----------|------|---------|-------------|
| gnomad41_exome | hg38_gnomad41_exome.txt.gz | 17 | gnomAD v4.1 exome frequencies |
| gnomad41_genome | hg38_gnomad41_genome.txt.gz | 13 | gnomAD v4.1 genome frequencies |
| exac03 | hg38_exac03.txt.gz | 8 | ExAC v0.3 frequencies |
| esp6500siv2_all | hg38_esp6500siv2_all.txt.gz | 1 | ESP6500 frequencies |
| 1000g2015aug_* | hg38_*.sites.2015_08.txt.gz | 1 each | 1000 Genomes (ALL, AFR, EAS, AMR, EUR, SAS) |
| kaviar_20150923 | hg38_kaviar_20150923.txt.gz | 3 | Kaviar frequencies |
| hrcr1 | hg38_hrcr1.txt.gz | 6 | Haplotype Reference Consortium |
| gme | hg38_gme.txt.gz | 8 | Greater Middle East frequencies |

### Pathogenicity and functional predictions

| Database | File | Columns | Description |
|----------|------|---------|-------------|
| ljb26_all | hg38_ljb26_all.txt.gz | 25 | SIFT, PolyPhen2, LRT, MutationTaster, FATHMM, MetaSVM, MetaLR, VEST4, CADD, GERP++, phyloP, SiPhy |
| intervar_20250721 | hg38_intervar_20250721.txt.gz | 29 | InterVar ACMG classification (2025) |
| clinvar_20250721 | hg38_clinvar_20250721.txt.gz | 5 | ClinVar clinical significance (2025) |
| mcap | hg38_mcap.txt.gz | 1 | M-CAP pathogenicity scores |
| revel | hg38_revel.txt.gz | 1 | REVEL ensemble scores |
| dbscsnv11 | hg38_dbscsnv11.txt.gz | 2 | Splice site predictions (ADA, RF) |

### Cancer databases

| Database | File | Columns | Description |
|----------|------|---------|-------------|
| cosmic | hg38_cosmic.txt.gz | 1 | COSMIC somatic mutations |
| icgc28 | hg38_icgc28.txt.gz | 2 | ICGC cancer mutations |
| TCGA | hg38_TCGA.txt.gz | 1 | TCGA cancer mutations |
| civic | hg38_civic.txt.gz | 1 | CIViC clinical interpretations |
| tumorportal | hg38_tumorportal.txt.gz | 1 | Tumor Portal mutations |

### Region and variant databases

| Database | File | Description |
|----------|------|-------------|
| snp141 | hg38_snp141.txt.gz | dbSNP 141 rsIDs |
| cytoBand | hg38_cytoBand.txt.gz | Cytogenetic bands |
| genomicSuperDups | hg38_genomicSuperDups.txt.gz | Segmental duplications |
| gwasCatalog | hg38_gwasCatalog.txt.gz | GWAS Catalog associations |
| encRegTfbsClustered | hg38_encRegTfbsClustered.txt.gz | ENCODE TF binding sites |
| wgEncodeRegDnaseClustered | hg38_wgEncodeRegDnaseClustered.txt.gz | ENCODE DNase hypersensitivity |
| wgRna | hg38_wgRna.txt.gz | snoRNA and miRNA annotations |
| dgvMerged | hg38_dgvMerged.txt.gz | Database of Genomic Variants |
| phastConsElements30way | hg38_phastConsElements30way.txt.gz | 30-way conservation |
| phastConsElements100way | hg38_phastConsElements100way.txt.gz | 100-way conservation |
| targetScanS | hg38_targetScanS.txt.gz | TargetScan miRNA target sites |
| cg69 | hg38_cg69.txt.gz | Complete Genomics 69 genomes |
| nci60 | hg38_nci60.txt.gz | NCI-60 cell line panel |

### Updating databases

To update a database (e.g., ClinVar to a newer release):

1. Prepare the new file in ANNOVAR's tabix-compatible format (bgzipped, tabix-indexed):
   ```bash
   bgzip new_database.txt
   tabix -s 1 -b 2 -e 3 new_database.txt.gz
   ```

2. Place the `.txt.gz` and `.txt.gz.tbi` files in your `HUMANDB_DIR`.

3. Update the entry in `config.py`:
   ```python
   'clinvar_20260101': {
       'file': 'hg38_clinvar_20260101.txt.gz',
       'type': 'filter',
       'operation': 'f',
       'format': 'generic',
       'columns': ['CLNALLELEID', 'CLNDN', 'CLNDISDB', 'CLNREVSTAT', 'CLNSIG'],
       'data_start_col': 5,
       'has_chr_prefix': False,
   },
   ```

4. Update the `DATABASE_ORDER` list and `OUTPUT_COLUMNS` list in `config.py` to reference the new key name.

## License

This project is provided as-is for research use. The annotation databases are subject to their respective licenses and terms of use.
