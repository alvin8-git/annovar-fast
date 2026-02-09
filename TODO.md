# ANNOVAR-fast TODO - Remaining Diffs from Merge.vcf Test

## Test Command
```bash
cd /data/alvin/annovar/annovar-fast
python3 annovar-fast.py /data/alvin/annovar/vcf-hg38/Merge.vcf -o /tmp/merge_test.txt
```

## Current Status (34 substantive diffs, 310 order-only)
**Progress: 188 → 34 substantive diffs** (reduced by 154)

The iSeq single-sample test (5 variants) passes perfectly (only Name= ordering diffs).
The Merge.vcf multi-sample test (192 variants, 95 samples) has remaining issues:

### Substantive Diffs
| Column | Diffs | Notes |
|--------|-------|-------|
| Gene.ensGene | 11 | Database version mismatch (skip for now) |
| GeneDetail.ensGene | 11 | Database version mismatch (skip for now) |
| AAChange.ensGene | 5 | Database version mismatch (skip for now) |
| AAChange.refGene | 3 | All HLA-B mRNA polymorphism (unfixable) |
| GeneDetail.refGene | 1 | HLA-B exonic;splicing detail missing |
| Func.refGene | 1 | HLA-B same-transcript exonic;splicing |
| snp141 | 1 | rs147324178 vs rs41563412 at chr6:31356965 |
| phastConsElements100way | 1 | Score/Name formatting (Score=345 vs 356) |

### Order-only Diffs (not bugs, just different ordering of equivalent values)
| Column | Diffs |
|--------|-------|
| dgvMerged | 131 |
| encRegTfbsClustered | 129 |
| cosmic | 35 |
| gwasCatalog | 13 |
| Gene.refGene | 2 |

### ensGene Database Version Mismatch (27 diffs - SKIP FOR NOW)
The tabix-compressed hg38_ensGene database has different transcript versions than the
flat file ANNOVAR used to generate expected output. Examples:
- KRAS: tabix has ENST00000556131.2 / ENST00000557334.6, expected uses .1 / .5
- MYD88: ENST00000396334.8 missing from tabix database entirely
- HLA-B: tabix has ENST00000640094.2, expected uses .1
- CEBPA region: extra ENSG00000267727 gene in tabix
**Resolution**: Regenerate tabix from same ensGene source, or accept as database version diffs.

## Remaining Issues to Fix

### 1. HLA-B same-transcript exonic;splicing (1 Func + 1 GeneDetail diff)
**Problem**: chr6:31355219 C>T is in exon4 of NM_005514 AND 97bp from exon5 (within
splicing_threshold=100). ANNOVAR classifies as "exonic;splicing" with GeneDetail showing
the splicing cDNA annotation. Our code classifies as just "exonic" with GeneDetail=".".

**Root cause**: `_classify_variant` returns 'exonic' when a variant is inside an exon,
without checking if it's also near a splice site of a DIFFERENT exon in the same transcript.

**Previous attempt**: Added per-transcript also_splicing check, but it caused a false positive
for TP53 chr17:7675070 (17bp inside exon7, but within 100bp of exon6 end). ANNOVAR doesn't
classify this as exonic;splicing because the variant is too deep in the exon (ANNOVAR uses
the exonicsplicing logic which only triggers when the variant is near the exon boundary it's IN).

**Investigation needed**: Study ANNOVAR lines 786-807 for exact exonicsplicing condition:
```perl
# ANNOVAR: if ($exession++ and ...) where $exession counts exonic matches
# The exonicsplicing check appears to be: if the variant is in an exon AND
# also within splicing_threshold of the SAME exon's boundary
```
This suggests ANNOVAR checks if the variant is within splicing_threshold of the boundary
of the exon it's IN, not a different exon. For chr6:31355219, the variant is 4bp from
exon4 end (31355223), which is within the 2bp canonical splice site. But ANNOVAR may also
check the 97bp distance to exon5.

### 2. AAChange.refGene HLA-B (3 diffs - likely UNFIXABLE)
All 3 are HLA-B block substitution protein annotations that differ because the mRNA
sequence for NM_005514#6#31353874 in our database differs from the alt contig mRNAs
that ANNOVAR may have used. The chr6 main assembly has different sequence content at
these positions due to HLA polymorphism.
- Line 51: chr6:31356181-31356183 TTG>GTC → p.D201_K202delinsET vs p.E201E
- Line 60: chr6:31356423-31356424 GC>AT → p.S121N vs p.R121N
- Line 89: chr6:31356825-31356827 TCT>ATC → p.E69M vs p.T69M

### 3. snp141 (1 diff)
chr6:31356965 T>C: GOT=rs147324178, EXP=rs41563412. May be database ordering - both
rsIDs may match but we return the first one found while ANNOVAR returns a different one.

### 4. phastConsElements100way (1 diff)
chrX:48791208-48791257 (50bp deletion): GOT=Score=345;Name=lod=39,lod=35,
EXP=Score=356;Name=lod=39. We're finding an extra overlapping region (lod=35) and the
score differs. May be a query range issue for large deletions.

### 5. Order-only diffs (310 total)
Region databases return overlapping regions in different order than ANNOVAR. ANNOVAR
processes the raw flat file sequentially, while our tabix queries may return records in
a different order. These are cosmetic differences - the same values are present.

**Possible fix**: Sort region results by genomic position or by the order they appear
in the database file.

## What Was Fixed (this session)
1. **SNP141 strand handling** (35 diffs fixed): Reverse complement alleles when strand='-'
2. **Exonic function priority** (MYD88 stoploss): Implemented exonic_buckets with priority system
3. **NFS deletion position** (CEBPA p.G139del): Protein comparison approach instead of formula
4. **Downstream/upstream classification** (7 F2 gene diffs): Added neargene=1000 check
5. **UTR cDNA annotations** (UTR3 c.* and UTR5 c.-): Added _make_utr_annotation method
6. **Splicing annotation ordering**: Sort by exon number
7. **Splicing+UTR combined detail**: Append UTR3/UTR5 after cDNA splicing annotations
8. **ncRNA_exonic;splicing classification**: Cross-transcript combined annotation

## Reference Files
- **ANNOVAR source**: `/data/alvin/annovar/Software/annovar/annotate_variation.pl`
  - Exonicsplicing check: lines 786-807
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

```bash
# Column-specific diff analysis (change column name as needed)
python3 /tmp/diff_analysis.py func    # Func/Gene/GeneDetail/ExonicFunc
python3 /tmp/diff_analysis.py aachange  # AAChange
python3 /tmp/diff_analysis.py snp141    # snp141
```
