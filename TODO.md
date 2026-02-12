# TODO

## Pending

### Database archive
- [ ] Upload `humandb-tbi-required.tar.gz` (29 GB) to file hosting
  - Archive location: `/data/alvin/annovar/humandb-tbi-required.tar.gz`
  - Contains all 39 databases (`.txt.gz` + `.tbi` pairs) + 2 mRNA FASTA files
  - Once hosted, update `README.md` — replace `<DOWNLOAD_URL>` placeholder in the Database Setup section

### Remaining annotation diffs (16 substantive, 310 order-only)

Tested against ANNOVAR output on Merge.vcf (192 variants, 95 samples):

| Category | Count | Details |
|----------|-------|---------|
| Database version (ensGene) | 12 | Gene.ensGene (11) + GeneDetail.ensGene (1) — different transcript set in hg38_ensGene_new vs original |
| HLA-B mRNA polymorphism | 3 | AAChange.refGene — unfixable, mRNA differs at HLA locus |
| snp141 tie-breaking | 1 | Both rsIDs valid, different strand entry chosen |
| Order-only (cosmetic) | 310 | Region DBs return same values in different order (dgvMerged 131, encRegTfbsClustered 129, cosmic 35, gwasCatalog 13, Gene.refGene 2) |

None of these are algorithmic bugs. The ensGene diffs would disappear with the same DB version; the rest are inherent to different tools querying the same data.

### Potential improvements
- [ ] Add `--sequential` CLI flag to `annovar-fast.py` (currently only available via Python API)
- [ ] Add `--workers N` CLI flag to override `MAX_WORKERS`
- [ ] Add `--db` CLI flag to override `HUMANDB_DIR` without editing `config.py`
- [ ] Stream output during parallel processing instead of buffering all results in memory (reduces peak RAM for WGS)
- [ ] Add progress logging (e.g., variants/sec, ETA) for long-running WGS jobs

## Completed

### Parallel annotation (2025-02-12)
- [x] Fork-based parallelism with equal-size chunking across all CPU cores
- [x] Workers inherit mRNA data via copy-on-write, re-open tabix handles
- [x] WES 75K variants: 23s on 64 cores (was ~10 min in ANNOVAR)
- [x] WGS 4.7M variants: 14 min on 64 cores (was 160 min in ANNOVAR)
- [x] Output identical to sequential mode (verified on WES + multi-sample)

### Database updates (2025-02-12)
- [x] gnomAD v2.1.1/v3.0 -> gnomAD v4.1 (exome + genome) with prefixed column names
- [x] ClinVar 20220730 -> 20250721
- [x] InterVar 20180118 -> 20250721
- [x] ljb26_all column names updated (RadialSVM->MetaSVM, LR->MetaLR, VEST3->VEST4, phyloP46way->phyloP30way)
- [x] ensGene -> ensGene_new

### Core annotation engine
- [x] VCF parsing via cyvcf2 with gzip support
- [x] ANNOVAR coordinate conversion (SNP, insertion, deletion, block substitution)
- [x] Gene annotation: exonic, intronic, splicing, UTR, upstream/downstream, intergenic, ncRNA
- [x] Amino acid change computation using mRNA FASTA
- [x] 37 filter/region database queries via tabix
- [x] Multi-sample VCF support

## Reference files
- ANNOVAR source: `/data/alvin/annovar/Software/annovar/annotate_variation.pl`
- Test VCFs: `/data/alvin/annovar/vcf-hg38/` (iSeq, T7_WES, Merge, NA12878_WGS)
- Database directory: `/data/alvin/annovar/humandb-tbi/`

## Debugging
```bash
# Run annotation and compare against ANNOVAR reference
python3 annovar-fast.py /data/alvin/annovar/vcf-hg38/Merge.vcf -o /tmp/merge_test.txt

# Detailed diff analysis with order detection
python3 << 'PYEOF'
import re
ref_lines = open('/data/alvin/annovar/vcf-hg38/Merge.annovar.txt').readlines()
got_lines = open('/tmp/merge_test.txt').readlines()
header = ref_lines[0].strip().split('\t')
from collections import Counter
diffs = Counter()
for i in range(1, min(len(ref_lines), len(got_lines))):
    rf = ref_lines[i].strip().split('\t')
    gf = got_lines[i].strip().split('\t')
    for j in range(min(len(rf), len(gf))):
        if rf[j] == gf[j]:
            continue
        col = header[j] if j < len(header) else f'col{j}'
        def normalize(s):
            s = re.sub(r'(Name|Score)=', '', s)
            return set(re.split(r'[,;]', s))
        if normalize(rf[j]) == normalize(gf[j]):
            diffs[col + ' (order)'] += 1
        else:
            diffs[col] += 1
for col, count in sorted(diffs.items(), key=lambda x: -x[1]):
    print(f'{col}: {count}')
PYEOF
```
