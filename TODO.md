# ANNOVAR-fast TODO - Remaining Diffs from Merge.vcf Test

## Test Command
```bash
cd /data/alvin/annovar/annovar-fast
python3 annovar-fast.py /data/alvin/annovar/vcf-hg38/Merge.vcf -o /tmp/merge_test.txt
```

## Current Status (16 substantive diffs, 310 order-only)
**Progress: 188 → 34 → 16 substantive diffs**

The iSeq single-sample test passes perfectly.
The Merge.vcf multi-sample test (192 variants, 95 samples) has remaining issues:

### Substantive Diffs
| Column | Diffs | Notes |
|--------|-------|-------|
| Gene.ensGene | 11 | Database version mismatch (extra CTD-2540B15.7/MIR6891 genes) |
| AAChange.refGene | 3 | All HLA-B mRNA polymorphism (unfixable) |
| GeneDetail.ensGene | 1 | Database version (different exon structure for ENST00000610292.4) |
| snp141 | 1 | rs147324178 vs rs41563412 at chr6:31356965 (tie-breaking) |

### Order-only Diffs (not bugs, just different ordering of equivalent values)
| Column | Diffs |
|--------|-------|
| dgvMerged | 131 |
| encRegTfbsClustered | 129 |
| cosmic | 35 |
| gwasCatalog | 13 |
| Gene.refGene | 2 |

### Categorization of remaining 16 substantive diffs
- **Database version diffs (12)**: Gene.ensGene (11) + GeneDetail.ensGene (1) — would be fixed by
  using the same ensGene database version as the expected output
- **HLA-B polymorphism (3)**: AAChange.refGene — unfixable, different mRNA sequence at HLA locus
- **snp141 tie-breaking (1)**: Both rsIDs are valid at the same position, just different strand entries

### Order-only diffs (310)
Region databases return overlapping regions in different order than ANNOVAR's internal
bin-based traversal. ANNOVAR uses UCSC binning which visits smaller/more-specific regions
first, while tabix returns by start position. These are cosmetic — same values, different order.

## What Was Fixed

### Session 1
1. **SNP141 strand handling** (35 diffs fixed): Reverse complement alleles when strand='-'
2. **Exonic function priority** (MYD88 stoploss): Implemented exonic_buckets with priority system
3. **NFS deletion position** (CEBPA p.G139del): Protein comparison approach instead of formula
4. **Downstream/upstream classification** (7 F2 gene diffs): Added neargene=1000 check
5. **UTR cDNA annotations** (UTR3 c.* and UTR5 c.-): Added _make_utr_annotation method
6. **Splicing annotation ordering**: Sort by exon number
7. **Splicing+UTR combined detail**: Append UTR3/UTR5 after cDNA splicing annotations
8. **ncRNA_exonic;splicing classification**: Cross-transcript combined annotation

### Session 2
9. **ensGene database restored**: Switched back to original hg38_ensGene, fixing AAChange.ensGene (5→0)
   and GeneDetail.ensGene (11→8→1)
10. **Upstream/downstream minimum distance**: When multiple transcripts of the same gene are
    near the variant, report the minimum distance (fixed 7 chr11 dist= diffs)
11. **Exonic;splicing same-transcript detection**: Added `_check_exonic_splicing` method that
    checks if an exonic variant is also near a splice site of an adjacent exon, using
    ANNOVAR's strand-aware loop order (+ strand: check lower-k exons; - strand: check higher-k exons).
    Fixed HLA-B chr6:31355219 without false positive on TP53 chr17:7675070.
12. **phastConsElements100way**: Fixed two bugs — tracked best_normscore alongside best_score,
    and report only the best-scoring region instead of all overlapping regions.

## Reference Files
- **ANNOVAR source**: `/data/alvin/annovar/Software/annovar/annotate_variation.pl`
  - Exonicsplicing check (exonicsplicing flag): lines 786-807
  - Exonicsplicing check (no flag, intronic only): lines 819-837
  - Minus strand exon loop: lines 1002-1174
  - Insertion logic: lines 1618-1700
  - Deletion logic: lines 1701-1905
  - Block substitution: lines 1731-1754, 1952-1967
  - SNV logic: lines 1755-1800
  - SNP141 matching: lines 2426-2500
  - Exonic function priority: lines 2018-2116
  - mRNA loading with multimap: lines 3321-3399, key line 3496
- **Test input**: `/data/alvin/annovar/vcf-hg38/Merge.vcf`
- **Expected output**: `/data/alvin/annovar/vcf-hg38/Merge.annovar.txt`
- **mRNA files**: `/data/alvin/Databases/humandb/hg38_refGeneMrna.fa`, `hg38_ensGeneMrna.fa`
  - Both refGene and ensGene mRNA files are used for computing AA changes (protein annotations)
  - Splice site classification uses only exon boundaries from the gene model, not mRNA sequences

## Debugging Tips
```bash
# Full diff count with order detection
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
