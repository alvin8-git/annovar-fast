# ANNOVAR-fast TODO - Remaining Diffs from Merge.vcf Test

## Test Command
```bash
cd /data/alvin/annovar/annovar-fast
python3 annovar-fast.py /data/alvin/annovar/vcf-hg38/Merge.vcf -o /tmp/merge_test.txt
```

## Current Status (188 substantive diffs, 308 order-only)
The iSeq single-sample test (5 variants) passes perfectly (only Name= ordering diffs).
The Merge.vcf multi-sample test (192 variants, 95 samples) has remaining issues:

| Column | Diffs | Notes |
|--------|-------|-------|
| AAChange.ensGene | 38 | Same issues as refGene but for ensGene transcripts |
| snp141 | 36 | False positive matches still occurring |
| AAChange.refGene | 33 | Insertion cDNA position, frameshift annotation |
| GeneDetail.refGene | 21 | Missing splicing annotations (large splicing_threshold=100) |
| Gene.ensGene | 18 | Likely downstream of Func/GeneDetail issues |
| GeneDetail.ensGene | 12 | Same as refGene issues |
| Func.refGene | 9 | Some variants not classified correctly |
| Gene.refGene | 9 | Same as Func issues |
| Func.ensGene | 7 | Same as refGene issues |
| ExonicFunc.refGene | 2 | Probably stoploss vs nonsynonymous edge cases |
| ExonicFunc.ensGene | 2 | Same |
| phastConsElements100way | 1 | Minor region annotation issue |

## Priority Issues to Fix

### 1. SNP141 false positives (36 diffs)
**Problem**: The new class-based matching is closer but still has false positives.
**Root cause**: Need to investigate specific cases. The `_encode_obs_for_query` function
may not handle all variant types correctly, or there may be strand-related issues in
the snp141 database (some entries have strand='-' requiring revcom of alleles).
**Investigation**:
- Check chr1:169549635 A>- (got rs199919250, expected rs55717622) - deletion matching issue
- Check if strand field (column 6) affects allele encoding for snp141

### 2. AAChange diffs (33 refGene + 38 ensGene)
Several sub-issues:

**a. Insertion cDNA position off by 1** (e.g., c.1493_1494ins vs c.1494_1495ins)
- The new code added strand-specific insertion position logic but needs verification
- Test with chr1:43349288 ->GT on NM_005373 (+ strand)
- Expected: c.1494_1495insGT (cds_start+1 to cds_start+2 for + strand)

**b. Insertion protein annotation** (e.g., p.L498delins vs p.H499Vfs*3)
- Frameshift insertions should show: `p.{next_aa}{varpos+1}fs*{distance_to_stop}`
- ANNOVAR computes the variant codon by inserting obs into wtnt3, translates it,
  then for frameshift: `p.{wtaa}{varpos}fs` but the expected shows next AA
- Need to re-examine the varpos calculation for insertions on + strand

**c. Stoploss detection** (line 29: got nonsynonymous SNV, expected stoploss)
- When wildtype AA is 'X' (stop codon), should return 'stoploss'
- Issue may be due to wrong mRNA lookup (multimap key issue)

### 3. GeneDetail splicing diffs (21 refGene + 12 ensGene)
**Problem**: Missing splicing annotations for variants near splice sites.
**Sub-issues**:

**a. Missing donor/acceptor annotations for large splicing_threshold**
- refGene uses splicing_threshold=100, so variants up to 100bp from exon boundaries
  should be detected as splicing
- Some variants are near two exons simultaneously (in a short intron)
- The fix to collect all splicing annotations should help, but some may still be missed
- Example: chr6:31355288 is both 29bp from exon4 end (donor) and 65bp from exon5 start (acceptor)

**b. Variant overlapping CDS boundary gives UTR annotation instead of splicing**
- The _splicing_plus/minus methods check `start >= tx.cdsstart` (+ strand) or
  `start <= tx.cdsend` (- strand) to decide between cDNA annotation vs UTR5
- This may be incorrectly classifying some CDS-adjacent splicing variants

### 4. Func/Gene diffs (9+9 refGene, 7+18 ensGene)
**Problem**: Some variants classified differently (e.g., exonic vs splicing, intronic vs splicing)
**Root cause**: The `_classify_variant` method processes exon overlap first, then splicing.
ANNOVAR may classify differently for variants that partially overlap both an exon and a
splice region. Need to investigate specific cases.

### 5. phastConsElements100way (1 diff)
Minor issue - probably a formatting or score rounding difference.

## Reference Files
- **ANNOVAR source**: `/data/alvin/annovar/Software/annovar/annotate_variation.pl`
  - Insertion logic: lines 1618-1700
  - Deletion logic: lines 1701-1905
  - Block substitution: lines 1731-1754, 1952-1967
  - SNV logic: lines 1755-1800
  - SNP141 matching: lines 2426-2500
  - mRNA loading with multimap: lines 3321-3399, key line 3496
- **Test input**: `/data/alvin/annovar/vcf-hg38/Merge.vcf`
- **Expected output**: `/data/alvin/annovar/vcf-hg38/Merge.annovar.txt`
- **mRNA files**: `/data/alvin/Databases/humandb/hg38_refGeneMrna.fa`, `hg38_ensGeneMrna.fa`

## Debugging Tips
```python
# Quick diff analysis
python3 -c "
ref_lines = open('/data/alvin/annovar/vcf-hg38/Merge.annovar.txt').readlines()
got_lines = open('/tmp/merge_test.txt').readlines()
header = ref_lines[0].strip().split('\t')
j = header.index('AAChange.refGene')  # Change column name as needed
for i in range(1, min(len(ref_lines), len(got_lines))):
    rf = ref_lines[i].strip().split('\t')
    gf = got_lines[i].strip().split('\t')
    if rf[j] != gf[j]:
        print(f'Line {i+1}: VAR={gf[0]}:{gf[1]}-{gf[2]} {gf[3]}>{gf[4]}')
        print(f'  GOT={gf[j][:120]}')
        print(f'  EXP={rf[j][:120]}')
        print()
"
```
