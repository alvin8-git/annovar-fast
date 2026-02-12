#!/usr/bin/env python3
"""
Gene annotation module for ANNOVAR-fast.
Handles refGene and ensGene annotation including:
- Functional classification (exonic, intronic, splicing, UTR, intergenic)
- Amino acid change computation for exonic variants
- HGVS notation for splicing variants
"""

import os
import re
import logging
from typing import Dict, List, Optional, Tuple

import pysam

from config import HUMANDB_DIR, NA_STRING, DATABASE_CONFIG

logger = logging.getLogger(__name__)

# Standard genetic code
CODON_TABLE = {
    'TTT': 'F', 'TTC': 'F', 'TTA': 'L', 'TTG': 'L',
    'CTT': 'L', 'CTC': 'L', 'CTA': 'L', 'CTG': 'L',
    'ATT': 'I', 'ATC': 'I', 'ATA': 'I', 'ATG': 'M',
    'GTT': 'V', 'GTC': 'V', 'GTA': 'V', 'GTG': 'V',
    'TCT': 'S', 'TCC': 'S', 'TCA': 'S', 'TCG': 'S',
    'CCT': 'P', 'CCC': 'P', 'CCA': 'P', 'CCG': 'P',
    'ACT': 'T', 'ACC': 'T', 'ACA': 'T', 'ACG': 'T',
    'GCT': 'A', 'GCC': 'A', 'GCA': 'A', 'GCG': 'A',
    'TAT': 'Y', 'TAC': 'Y', 'TAA': 'X', 'TAG': 'X',
    'CAT': 'H', 'CAC': 'H', 'CAA': 'Q', 'CAG': 'Q',
    'AAT': 'N', 'AAC': 'N', 'AAA': 'K', 'AAG': 'K',
    'GAT': 'D', 'GAC': 'D', 'GAA': 'E', 'GAG': 'E',
    'TGT': 'C', 'TGC': 'C', 'TGA': 'X', 'TGG': 'W',
    'CGT': 'R', 'CGC': 'R', 'CGA': 'R', 'CGG': 'R',
    'AGT': 'S', 'AGC': 'S', 'AGA': 'R', 'AGG': 'R',
    'GGT': 'G', 'GGC': 'G', 'GGA': 'G', 'GGG': 'G',
}

AA_3TO1 = {
    'A': 'A', 'R': 'R', 'N': 'N', 'D': 'D', 'C': 'C',
    'E': 'E', 'Q': 'Q', 'G': 'G', 'H': 'H', 'I': 'I',
    'L': 'L', 'K': 'K', 'M': 'M', 'F': 'F', 'P': 'P',
    'S': 'S', 'T': 'T', 'W': 'W', 'Y': 'Y', 'V': 'V',
    'X': 'X',
}


def revcom(seq):
    """Reverse complement a DNA sequence."""
    comp = str.maketrans('ACGTacgt', 'TGCAtgca')
    return seq.translate(comp)[::-1]


def translate_dna(seq):
    """Translate DNA to protein."""
    seq = seq.upper()
    protein = []
    for i in range(0, len(seq) - 2, 3):
        codon = seq[i:i+3]
        aa = CODON_TABLE.get(codon, '?')
        protein.append(aa)
    return ''.join(protein)


class Transcript:
    """Represents a gene transcript."""
    __slots__ = ['name', 'name2', 'chrom', 'strand', 'txstart', 'txend',
                 'cdsstart', 'cdsend', 'exonstarts', 'exonends', 'exoncount']

    def __init__(self, fields):
        """Parse a refGene/ensGene record (already split by tab).
        Format: bin, name, chr, strand, txStart, txEnd, cdsStart, cdsEnd, exonCount, exonStarts, exonEnds, score, name2, ...
        """
        self.name = fields[1]
        self.chrom = fields[2]
        self.strand = fields[3]
        # Convert 0-based to 1-based (matching ANNOVAR: $txstart++; $cdsstart++; map {$_++} @exonstart)
        self.txstart = int(fields[4]) + 1
        self.txend = int(fields[5])  # exonEnds are already correct (1-based inclusive)
        self.cdsstart = int(fields[6]) + 1
        self.cdsend = int(fields[7])
        self.exoncount = int(fields[8])
        # Parse and convert exon starts (0-based -> 1-based)
        starts = [int(x) + 1 for x in fields[9].rstrip(',').split(',') if x]
        ends = [int(x) for x in fields[10].rstrip(',').split(',') if x]
        self.exonstarts = starts
        self.exonends = ends
        self.name2 = fields[12] if len(fields) > 12 else self.name


class GeneAnnotator:
    """Annotates variants with gene-level information."""

    def __init__(self, db_name: str):
        """Initialize for refGene or ensGene."""
        self.db_name = db_name
        config = DATABASE_CONFIG[db_name]
        self.columns = config['columns']
        self.splicing_threshold = config.get('splicing_threshold', 2)
        self.mrna_seqs = {}

        # Open tabix file
        gene_file = os.path.join(HUMANDB_DIR, config['file'])
        self.tabix = None
        if os.path.exists(gene_file) and os.path.exists(gene_file + '.tbi'):
            self.tabix = pysam.TabixFile(gene_file)

        # Load mRNA sequences
        mrna_file = os.path.join(HUMANDB_DIR, f"hg38_{db_name}Mrna.fa")
        if os.path.exists(mrna_file):
            self._load_mrna(mrna_file)

    def _load_mrna(self, fasta_path: str):
        """Load mRNA sequences from FASTA file.

        Uses composite keys matching ANNOVAR's multimap convention:
        For entries with 'leftmost exon at chrX:POS', the key becomes 'NAME#CHR#POS'
        (chr prefix stripped). This handles transcripts that map to multiple locations.
        """
        logger.info(f"Loading mRNA sequences from {fasta_path}")
        seqid = None
        curseq = []
        count = 0

        import re
        multimap_re = re.compile(r'leftmost exon at (?:chr)?([\w.]+):(\d+)')

        with open(fasta_path) as f:
            for line in f:
                line = line.rstrip('\n\r')
                if line.startswith('>'):
                    if seqid and curseq:
                        self.mrna_seqs[seqid] = ''.join(curseq)
                        count += 1
                    # Parse sequence ID with multimap support
                    parts = line[1:].split(None, 1)
                    seqid = parts[0]
                    # Check for multimap annotation
                    m = multimap_re.search(line)
                    if m:
                        seqid = f"{seqid}#{m.group(1)}#{m.group(2)}"
                    curseq = []
                else:
                    curseq.append(line.upper())

        if seqid and curseq:
            self.mrna_seqs[seqid] = ''.join(curseq)
            count += 1

        logger.info(f"Loaded {count} mRNA sequences for {self.db_name}")

    def annotate(self, chrom: str, start: int, end: int, ref: str, alt: str) -> Dict[str, str]:
        """Annotate a variant with gene information."""
        na = {col: NA_STRING for col in self.columns}
        if not self.tabix:
            return na

        # Get overlapping transcripts
        transcripts = self._get_transcripts(chrom, start, end)
        if not transcripts:
            return self._format_result('intergenic', '.', '.', '.', '.')

        # Classify variant for each transcript
        # Priority: exonic > splicing > ncRNA_exonic > UTR > intronic > upstream/downstream
        best_func = None
        gene_names = set()
        gene_detail_parts = []
        utr_detail_parts = []
        exonic_found = False
        splicing_found = False
        intronic_found = False
        downstream_genes = {}  # gene_name -> detail
        upstream_genes = {}
        utr_found = False

        # ANNOVAR buckets exonic annotations by function type and picks highest priority
        # Priority: fsins > fsdel > fssub > stopgain > stoploss > nfsins > nfsdel > nfssub > nssnv > ssnv
        _EXONIC_PRIORITY = {
            'frameshift insertion': 0, 'frameshift deletion': 1, 'frameshift substitution': 2,
            'stopgain': 3, 'stoploss': 4,
            'nonframeshift insertion': 5, 'nonframeshift deletion': 6, 'nonframeshift substitution': 7,
            'nonsynonymous SNV': 8, 'synonymous SNV': 9, 'unknown': 10,
        }
        exonic_buckets = {}  # exonic_func -> list of aa_change strings

        for tx in transcripts:
            func, detail = self._classify_variant(tx, start, end, ref, alt)

            if func == 'exonic':
                exonic_found = True
                gene_names.add(tx.name2)
                # Compute AA change
                ef, ac = self._compute_exonic_annotation(tx, start, end, ref, alt)
                if ef != '.' and ac != '.':
                    exonic_buckets.setdefault(ef, []).append(ac)
                # ANNOVAR also checks if an exonic variant is near a splice site
                # of another exon in the same transcript (exonicsplicing)
                splice_detail = self._check_exonic_splicing(tx, start, end, ref, alt)
                if splice_detail:
                    splicing_found = True
                    gene_detail_parts.append(splice_detail)
            elif func == 'splicing':
                splicing_found = True
                gene_names.add(tx.name2)
                if detail:
                    gene_detail_parts.append(detail)
            elif func == 'intronic':
                intronic_found = True
                gene_names.add(tx.name2)
            elif func in ('UTR5', 'UTR3'):
                utr_found = True
                gene_names.add(tx.name2)
                if detail:
                    utr_detail_parts.append(detail)
                if best_func is None or best_func in ('intronic',):
                    best_func = func
            elif func == 'ncRNA_exonic':
                gene_names.add(tx.name2)
                if best_func is None:
                    best_func = func
            elif func == 'downstream':
                # Keep the minimum distance across transcripts of the same gene
                if tx.name2 not in downstream_genes or self._parse_dist(detail) < self._parse_dist(downstream_genes[tx.name2]):
                    downstream_genes[tx.name2] = detail
            elif func == 'upstream':
                if tx.name2 not in upstream_genes or self._parse_dist(detail) < self._parse_dist(upstream_genes[tx.name2]):
                    upstream_genes[tx.name2] = detail

        # Determine final function based on priority
        # ANNOVAR: genic annotations suppress upstream/downstream
        foundgenic = exonic_found or splicing_found or intronic_found or utr_found or (best_func is not None)

        if exonic_found and splicing_found:
            final_func = 'exonic;splicing'
        elif exonic_found:
            final_func = 'exonic'
        elif best_func == 'ncRNA_exonic' and splicing_found:
            final_func = 'ncRNA_exonic;splicing'
        elif splicing_found:
            final_func = 'splicing'
        elif best_func:
            final_func = best_func
        elif intronic_found:
            final_func = 'intronic'
        elif not foundgenic and downstream_genes:
            final_func = 'downstream'
            gene_names = set(downstream_genes.keys())
            gene_detail_parts = list(downstream_genes.values())
        elif not foundgenic and upstream_genes:
            final_func = 'upstream'
            gene_names = set(upstream_genes.keys())
            gene_detail_parts = list(upstream_genes.values())
        else:
            final_func = 'intergenic'

        gene_name = ';'.join(sorted(gene_names)) if gene_names else '.'

        # UTR details only appear in GeneDetail when final func is UTR (not when exonic/splicing)
        if final_func in ('UTR5', 'UTR3'):
            gene_detail_parts.extend(utr_detail_parts)
        gene_detail = ';'.join(gene_detail_parts) if gene_detail_parts else '.'

        # Pick highest-priority exonic function and its AA changes
        exonic_func = '.'
        aa_change = '.'
        if exonic_buckets:
            best_ef = min(exonic_buckets.keys(), key=lambda x: _EXONIC_PRIORITY.get(x, 99))
            exonic_func = best_ef
            aa_change = ','.join(exonic_buckets[best_ef])

        if final_func in ('splicing', 'intronic', 'downstream', 'upstream',
                          'ncRNA_exonic;splicing'):
            exonic_func = '.'
            aa_change = '.'

        # Sort AA changes by exon number then transcript name (matching ANNOVAR's sortExonicAnnotation)
        if aa_change != '.':
            aa_change = self._sort_exonic_annotation(aa_change)

        return self._format_result(final_func, gene_name, gene_detail, exonic_func, aa_change)

    def _get_transcripts(self, chrom: str, start: int, end: int) -> List[Transcript]:
        """Get all transcripts overlapping the variant region."""
        transcripts = []
        try:
            # Query with padding for splicing threshold
            padding = max(self.splicing_threshold, 1000)
            qstart = max(0, start - padding - 1)
            qend = end + padding

            for record in self.tabix.fetch(chrom, qstart, qend):
                fields = record.split('\t')
                if len(fields) < 13:
                    continue
                tx = Transcript(fields)
                # Check if variant is within or near this transcript
                if start <= tx.txend + padding and end >= tx.txstart - padding:
                    transcripts.append(tx)
        except (ValueError, StopIteration):
            pass
        return transcripts

    def _classify_variant(self, tx: Transcript, start: int, end: int, ref: str, alt: str) -> Tuple[str, str]:
        """Classify a variant relative to a transcript.
        Returns (function_type, detail_annotation).
        """
        neargene = 1000  # ANNOVAR default

        # Check if variant is outside transcript boundaries
        if end < tx.txstart:
            dist = tx.txstart - end
            if dist <= neargene:
                if tx.strand == '+':
                    return ('upstream', f"dist={dist}")
                else:
                    return ('downstream', f"dist={dist}")
            return ('intergenic', '')
        if start > tx.txend:
            dist = start - tx.txend
            if dist <= neargene:
                if tx.strand == '+':
                    return ('downstream', f"dist={dist}")
                else:
                    return ('upstream', f"dist={dist}")
            return ('intergenic', '')

        # Check if this is a non-coding transcript (cdsstart == cdsend)
        is_ncrna = (tx.cdsstart >= tx.cdsend)

        # Iterate through exons to classify
        exon_result = None  # Track if variant is in an exon
        exon_overlap_k = -1  # Which exon it overlaps
        for k in range(tx.exoncount):
            exon_start = tx.exonstarts[k]
            exon_end = tx.exonends[k]

            # Check if variant overlaps this exon
            if start >= exon_start and end <= exon_end:
                exon_overlap_k = k
                # Variant is within this exon
                if is_ncrna:
                    exon_result = ('ncRNA_exonic', '')
                    break

                # Check if in CDS or UTR
                if start >= tx.cdsstart and end <= tx.cdsend:
                    exon_result = ('exonic', '')
                    break
                elif end < tx.cdsstart:
                    utr_type = 'UTR5' if tx.strand == '+' else 'UTR3'
                    detail = self._make_utr_annotation(tx, start, ref, alt, utr_type)
                    return (utr_type, detail)
                elif start > tx.cdsend:
                    utr_type = 'UTR3' if tx.strand == '+' else 'UTR5'
                    detail = self._make_utr_annotation(tx, start, ref, alt, utr_type)
                    return (utr_type, detail)
                else:
                    # Spans CDS boundary - treat as exonic
                    exon_result = ('exonic', '')
                    break

            # Check if variant partially overlaps exon (for indels)
            if start <= exon_end and end >= exon_start:
                exon_overlap_k = k
                if is_ncrna:
                    exon_result = ('ncRNA_exonic', '')
                    break
                # Check CDS overlap
                if start <= tx.cdsend and end >= tx.cdsstart:
                    exon_result = ('exonic', '')
                    break
                elif end < tx.cdsstart:
                    utr_type = 'UTR5' if tx.strand == '+' else 'UTR3'
                    return (utr_type, '')
                elif start > tx.cdsend:
                    utr_type = 'UTR3' if tx.strand == '+' else 'UTR5'
                    return (utr_type, '')
                else:
                    exon_result = ('exonic', '')
                    break

        if exon_result is not None:
            return exon_result

        # Variant is in intron - check for splicing
        # Collect ALL splicing annotations (variant may overlap multiple exon splice sites)
        splicing_details = []
        for k in range(tx.exoncount):
            exon_start = tx.exonstarts[k]
            exon_end = tx.exonends[k]

            # Check splicing at acceptor site (before exon start)
            if k > 0 or tx.exoncount == 1:
                if self._has_overlap(start, end, exon_start - self.splicing_threshold, exon_start - 1):
                    detail = self._make_splicing_annotation(tx, k, start, end, ref, alt)
                    if detail:
                        splicing_details.append(detail)

            # Check splicing at donor site (after exon end)
            if k < tx.exoncount - 1 or tx.exoncount == 1:
                if self._has_overlap(start, end, exon_end + 1, exon_end + self.splicing_threshold):
                    detail = self._make_splicing_annotation(tx, k, start, end, ref, alt)
                    if detail:
                        splicing_details.append(detail)

        if splicing_details:
            # Sort by exon number (ANNOVAR puts lower exon numbers first)
            def _exon_sort_key(s):
                m = re.search(r'exon(\d+)', s)
                return int(m.group(1)) if m else 0
            splicing_details.sort(key=_exon_sort_key)

            # ANNOVAR also checks if the variant is in UTR territory and appends UTR annotation
            # after each splicing detail for that transcript/exon
            if not is_ncrna:
                utr_type = None
                if tx.strand == '+':
                    if start > tx.cdsend:
                        utr_type = 'UTR3'
                    elif end < tx.cdsstart:
                        utr_type = 'UTR5'
                else:
                    if end < tx.cdsstart:
                        utr_type = 'UTR3'
                    elif start > tx.cdsend:
                        utr_type = 'UTR5'
                if utr_type:
                    # Append UTR annotation after each splicing detail that has cDNA info
                    # (not when the detail already IS a UTR annotation)
                    augmented = []
                    for detail in splicing_details:
                        augmented.append(detail)
                        if ':c.' in detail:
                            m = re.match(r'([^:]+:exon\d+)', detail)
                            if m:
                                augmented.append(f"{m.group(1)}:{utr_type}")
                    splicing_details = augmented

            return ('splicing', ';'.join(splicing_details))

        # If we're inside the transcript but not in any exon and not near splice sites
        if start >= tx.txstart and end <= tx.txend:
            return ('intronic', '')

        return ('intergenic', '')

    def _check_exonic_splicing(self, tx, start: int, end: int, ref: str, alt: str) -> str:
        """Check if an exonic variant is also near a splice site of another exon.

        ANNOVAR iterates exons in order (low-to-high for + strand, high-to-low for
        - strand) and breaks when it finds the exonic match. So splicing is only
        checked for exons visited BEFORE the variant's exon in that order:
        - For + strand: only exons with k < variant_exon
        - For - strand: only exons with k > variant_exon
        """
        # Find which exon the variant is in
        variant_exon = -1
        for k in range(tx.exoncount):
            if start >= tx.exonstarts[k] and end <= tx.exonends[k]:
                variant_exon = k
                break
            if start <= tx.exonends[k] and end >= tx.exonstarts[k]:
                variant_exon = k
                break

        if variant_exon < 0:
            return ''

        # Only check exons that ANNOVAR would visit before the variant's exon
        if tx.strand == '+':
            exon_range = range(0, variant_exon)
        else:
            exon_range = range(tx.exoncount - 1, variant_exon, -1)

        for k in exon_range:
            exon_start = tx.exonstarts[k]
            exon_end = tx.exonends[k]

            # Check acceptor site (before exon start)
            if k > 0 or tx.exoncount == 1:
                if self._has_overlap(start, end, exon_start - self.splicing_threshold, exon_start - 1):
                    detail = self._make_splicing_annotation(tx, k, start, end, ref, alt)
                    if detail:
                        return detail

            # Check donor site (after exon end)
            if k < tx.exoncount - 1 or tx.exoncount == 1:
                if self._has_overlap(start, end, exon_end + 1, exon_end + self.splicing_threshold):
                    detail = self._make_splicing_annotation(tx, k, start, end, ref, alt)
                    if detail:
                        return detail

        return ''

    def _has_overlap(self, start1, end1, start2, end2):
        """Check if two ranges overlap."""
        return start1 <= end2 and end1 >= start2

    @staticmethod
    def _parse_dist(detail: str) -> int:
        """Extract numeric distance from 'dist=NNN' string."""
        if detail and detail.startswith('dist='):
            try:
                return int(detail[5:])
            except ValueError:
                pass
        return float('inf')

    def _make_utr_annotation(self, tx: Transcript, start: int, ref: str, alt: str, utr_type: str) -> str:
        """Generate UTR annotation like NM_000506:c.*97G>A (UTR3) or NM_005514:c.-18G>A (UTR5).

        Computes the distance from the CDS boundary along the spliced mRNA.
        """
        if start != start:  # indels - skip detailed annotation
            return ''

        if utr_type == 'UTR3':
            if tx.strand == '+':
                # Count exonic bases from cdsend+1 to variant
                dist = 0
                for k in range(tx.exoncount):
                    es, ee = tx.exonstarts[k], tx.exonends[k]
                    if ee <= tx.cdsend:
                        continue
                    region_start = max(es, tx.cdsend + 1)
                    if start >= region_start and start <= ee:
                        dist += start - region_start + 1
                        break
                    elif start > ee:
                        dist += ee - region_start + 1
                if dist > 0:
                    return f"{tx.name}:c.*{dist}{ref}>{alt}"
            else:
                # - strand: UTR3 is genomically before cdsstart
                dist = 0
                for k in range(tx.exoncount - 1, -1, -1):
                    es, ee = tx.exonstarts[k], tx.exonends[k]
                    if es >= tx.cdsstart:
                        continue
                    region_end = min(ee, tx.cdsstart - 1)
                    if start >= es and start <= region_end:
                        dist += region_end - start + 1
                        break
                    elif start < es:
                        dist += region_end - es + 1
                if dist > 0:
                    return f"{tx.name}:c.*{dist}{revcom(ref)}>{revcom(alt)}"
        elif utr_type == 'UTR5':
            if tx.strand == '+':
                # Count exonic bases from variant to cdsstart-1
                dist = 0
                for k in range(tx.exoncount - 1, -1, -1):
                    es, ee = tx.exonstarts[k], tx.exonends[k]
                    if es >= tx.cdsstart:
                        continue
                    region_end = min(ee, tx.cdsstart - 1)
                    if start >= es and start <= region_end:
                        dist += region_end - start + 1
                        break
                    elif start < es:
                        dist += region_end - es + 1
                if dist > 0:
                    return f"{tx.name}:c.-{dist}{ref}>{alt}"
            else:
                # - strand: UTR5 is genomically after cdsend
                dist = 0
                for k in range(tx.exoncount):
                    es, ee = tx.exonstarts[k], tx.exonends[k]
                    if ee <= tx.cdsend:
                        continue
                    region_start = max(es, tx.cdsend + 1)
                    if start >= region_start and start <= ee:
                        dist += start - region_start + 1
                        break
                    elif start > ee:
                        dist += ee - region_start + 1
                if dist > 0:
                    return f"{tx.name}:c.-{dist}{revcom(ref)}>{revcom(alt)}"
        return ''

    def _make_splicing_annotation(self, tx: Transcript, exon_idx: int, start: int, end: int, ref: str, alt: str) -> str:
        """Generate splicing annotation like NM_004972:exon12:c.1514-88G>A

        Follows ANNOVAR's exact algorithm from annotate_variation.pl:
        - For + strand: lenexon accumulates (exonend-exonstart+1) per CDS exon,
          reset at the exon containing cdsstart.
        - For acceptor site: lenexon -= (exonend[k]-exonstart[k]) (NOT +1)
        - For donor site: lenexon is used as-is
        - For - strand: iterate from last to first exon, similar logic with revcom.
        """
        exon_start = tx.exonstarts[exon_idx]
        exon_end = tx.exonends[exon_idx]

        # For indels, use r.spl notation with correct exon display number
        if start != end:
            display_exon = (tx.exoncount - exon_idx) if tx.strand == '-' else (exon_idx + 1)
            if end >= exon_start - self.splicing_threshold and end < exon_start:
                return f"{tx.name}:exon{display_exon}:r.spl"
            elif start > exon_end and start <= exon_end + self.splicing_threshold:
                return f"{tx.name}:exon{display_exon}:r.spl"
            return ''

        # Non-coding transcripts: return simple r.spl without cDNA details
        is_ncrna = (tx.cdsstart >= tx.cdsend)

        if tx.strand == '+':
            return self._splicing_plus(tx, exon_idx, start, ref, alt, is_ncrna)
        else:
            return self._splicing_minus(tx, exon_idx, start, ref, alt, is_ncrna)

    def _splicing_plus(self, tx: Transcript, exon_idx: int, start: int, ref: str, alt: str, is_ncrna: bool = False) -> str:
        """Splicing annotation for + strand, matching ANNOVAR's exact lenexon logic."""
        if is_ncrna:
            return ''  # Non-coding transcripts don't get cDNA splicing annotations

        lenexon = 0
        for k in range(exon_idx + 1):
            es = tx.exonstarts[k]
            ee = tx.exonends[k]
            lenexon += (ee - es + 1)
            if tx.cdsstart >= es and tx.cdsstart <= ee:
                # CDS starts in this exon - reset lenexon
                lenexon = ee - tx.cdsstart + 1

        exon_start = tx.exonstarts[exon_idx]
        exon_end = tx.exonends[exon_idx]

        if start >= exon_start - self.splicing_threshold and start < exon_start:
            # Acceptor site (before exon start)
            dist = exon_start - start
            if start >= tx.cdsstart:
                # ANNOVAR: lenexon -= (exonend[k]-exonstart[k]) -- note: NOT +1
                lenexon -= (exon_end - exon_start)
                return f"{tx.name}:exon{exon_idx+1}:c.{lenexon}-{dist}{ref}>{alt}"
            else:
                return f"{tx.name}:exon{exon_idx+1}:UTR5"
        elif start > exon_end and start <= exon_end + self.splicing_threshold:
            # Donor site (after exon end)
            dist = start - exon_end
            if start >= tx.cdsstart:
                return f"{tx.name}:exon{exon_idx+1}:c.{lenexon}+{dist}{ref}>{alt}"
            else:
                return f"{tx.name}:exon{exon_idx+1}:UTR5"
        return ''

    def _splicing_minus(self, tx: Transcript, exon_idx: int, start: int, ref: str, alt: str, is_ncrna: bool = False) -> str:
        """Splicing annotation for - strand, matching ANNOVAR's logic.

        ANNOVAR line 1078-1084: For - strand (iterating from last exon to first):
        - start < exonstart[k] (before exon in genomic) → c.$lenexon+dist (donor)
        - start > exonend[k] (after exon in genomic) → c.$lenexon-dist (acceptor), lenexon adjusted
        """
        if is_ncrna:
            return ''  # Non-coding transcripts don't get cDNA splicing annotations

        n = tx.exoncount

        lenexon = 0
        # Iterate from last exon backwards to exon_idx (matching ANNOVAR's - strand loop)
        for k in range(n - 1, exon_idx - 1, -1):
            es = tx.exonstarts[k]
            ee = tx.exonends[k]
            lenexon += (ee - es + 1)
            if tx.cdsend >= es and tx.cdsend <= ee:
                lenexon = tx.cdsend - es + 1

        exon_start = tx.exonstarts[exon_idx]
        exon_end = tx.exonends[exon_idx]
        display_exon = n - exon_idx  # ANNOVAR: @exonstart-1-$k+1

        if start >= exon_start - self.splicing_threshold and start < exon_start:
            # Before exon start (genomic) = donor site on - strand → c.X+Y
            dist = exon_start - start
            if start <= tx.cdsend:
                return f"{tx.name}:exon{display_exon}:c.{lenexon}+{dist}{revcom(ref)}>{revcom(alt)}"
            else:
                return f"{tx.name}:exon{display_exon}:UTR5"
        elif start > exon_end and start <= exon_end + self.splicing_threshold:
            # After exon end (genomic) = acceptor site on - strand → c.X-Y
            dist = start - exon_end
            if start <= tx.cdsend:
                lenexon -= (exon_end - exon_start)
                return f"{tx.name}:exon{display_exon}:c.{lenexon}-{dist}{revcom(ref)}>{revcom(alt)}"
            else:
                return f"{tx.name}:exon{display_exon}:UTR5"
        return ''

    def _get_mrna_key(self, tx: Transcript) -> str:
        """Build composite mRNA key matching ANNOVAR's multimap convention.
        Key format: name#chr(no prefix)#txstart(0-based)
        """
        chrom = tx.chrom.replace('chr', '')
        txstart_0based = tx.txstart - 1  # convert back to 0-based
        composite = f"{tx.name}#{chrom}#{txstart_0based}"
        if composite in self.mrna_seqs:
            return composite
        # Fall back to plain name for transcripts without multimap
        return tx.name

    def _compute_exonic_annotation(self, tx: Transcript, start: int, end: int, ref: str, alt: str) -> Tuple[str, str]:
        """Compute exonic function and AA change annotation.
        Returns (exonic_func, aa_change_string).
        Uses ANNOVAR's exact algorithm for refvarstart/refvarend/refcdsstart.
        """
        mrna_seq = self.mrna_seqs.get(self._get_mrna_key(tx))
        if not mrna_seq:
            return ('.', '.')

        # Find which exon the variant is in
        exon_pos = 0
        for k in range(tx.exoncount):
            if start <= tx.exonends[k] and end >= tx.exonstarts[k]:
                exon_pos = k + 1
                break
        if exon_pos == 0:
            return ('.', '.')

        # Calculate refvarstart, refvarend, refcdsstart using ANNOVAR's exact algorithm
        if tx.strand == '+':
            refvarstart, refvarend, refcdsstart = self._calc_mrna_positions_plus(tx, start, end)
        else:
            refvarstart, refvarend, refcdsstart = self._calc_mrna_positions_minus(tx, start, end)
            exon_pos = tx.exoncount - (exon_pos - 1)  # reverse exon numbering

        if refvarstart == 0 or refvarend == 0 or refcdsstart == 0:
            return ('.', '.')

        cds_start = refvarstart - refcdsstart  # 0-based in CDS
        cds_end = refvarend - refcdsstart

        if cds_start < 0:
            return ('.', '.')

        # Reverse complement ref/obs for minus strand (ANNOVAR line 1612-1614)
        if tx.strand == '-':
            obs = revcom(alt) if alt != '-' else '-'
            qref = revcom(ref) if ref != '-' else '-'
        else:
            obs = alt
            qref = ref

        # Wildtype codon at variant start
        fs = cds_start % 3
        wt_codon_start = refvarstart - fs - 1  # 0-based in mRNA
        if wt_codon_start < 0 or wt_codon_start + 3 > len(mrna_seq):
            return ('.', '.')
        wtnt3 = mrna_seq[wt_codon_start:wt_codon_start + 3]
        if len(wtnt3) < 3:
            return ('.', '.')
        wtaa = translate_dna(wtnt3)
        varpos = cds_start // 3 + 1

        if refvarstart == refvarend:
            # Single position in mRNA
            if qref == '-':
                # Insertion
                return self._annotate_insertion(tx, mrna_seq, exon_pos, refvarstart, refcdsstart,
                                                cds_start, fs, wtnt3, wtaa, varpos, obs)
            elif obs == '-':
                # Single nucleotide deletion
                return self._annotate_single_deletion(tx, mrna_seq, exon_pos, refvarstart, refcdsstart,
                                                      cds_start, fs, wtnt3, wtaa, varpos)
            elif len(obs) > 1:
                # Block substitution (single pos to multi-base)
                canno = f"c.{cds_start+1}delins{obs}"
                ref_len = refvarend - refvarstart + 1
                if (ref_len - len(obs)) % 3 == 0:
                    exonic_func = 'nonframeshift substitution'
                else:
                    exonic_func = 'frameshift substitution'
                panno = self._compute_block_sub_panno(mrna_seq, refvarstart, refvarend, refcdsstart, obs, fs, cds_start)
                return (exonic_func, f"{tx.name2}:{tx.name}:exon{exon_pos}:{canno}{panno}")
            else:
                # SNV
                return self._annotate_snv(tx, exon_pos, cds_start, fs, wtnt3, wtaa, varpos, obs, qref)
        else:
            # Multiple positions in mRNA
            if obs == '-':
                # Multi-nucleotide deletion
                return self._annotate_deletion(tx, mrna_seq, exon_pos, refvarstart, refvarend,
                                               refcdsstart, cds_start, cds_end, fs, wtnt3, wtaa, varpos)
            else:
                # Block substitution (multi-base)
                canno = f"c.{cds_start+1}_{cds_end+1}delins{obs}"
                ref_len = refvarend - refvarstart + 1
                if (ref_len - len(obs)) % 3 == 0:
                    exonic_func = 'nonframeshift substitution'
                else:
                    exonic_func = 'frameshift substitution'
                panno = self._compute_block_sub_panno(mrna_seq, refvarstart, refvarend, refcdsstart, obs, fs, cds_start)
                return (exonic_func, f"{tx.name2}:{tx.name}:exon{exon_pos}:{canno}{panno}")

    def _annotate_snv(self, tx, exon_pos, cds_start, fs, wtnt3, wtaa, varpos, obs, qref):
        """Annotate a single nucleotide variant."""
        codon_offset = fs
        var_codon = list(wtnt3)
        var_codon[codon_offset] = obs
        var_codon = ''.join(var_codon)
        var_aa = translate_dna(var_codon)

        if wtaa == var_aa:
            exonic_func = 'synonymous SNV'
        elif var_aa == 'X':
            exonic_func = 'stopgain'
        elif wtaa == 'X':
            exonic_func = 'stoploss'
        else:
            exonic_func = 'nonsynonymous SNV'

        canno = f"c.{wtnt3[codon_offset]}{cds_start+1}{obs}"
        panno = f"p.{wtaa}{varpos}{var_aa}"
        return (exonic_func, f"{tx.name2}:{tx.name}:exon{exon_pos}:{canno}:{panno}")

    def _annotate_insertion(self, tx, mrna_seq, exon_pos, refvarstart, refcdsstart,
                            cds_start, fs, wtnt3, wtaa, varpos, obs):
        """Annotate an insertion variant."""
        ins_len = len(obs)
        if ins_len % 3 == 0:
            exonic_func = 'nonframeshift insertion'
        else:
            exonic_func = 'frameshift insertion'

        # cDNA annotation: check for dup, handle strand-specific position
        if tx.strand == '+':
            # Check for duplication
            if obs == mrna_seq[refvarstart - 1:refvarstart - 1 + len(obs)] if refvarstart - 1 + len(obs) <= len(mrna_seq) else False:
                canno = f"c.{cds_start+1}dup{obs}"
            elif refvarstart < len(mrna_seq) and obs == mrna_seq[refvarstart:refvarstart + len(obs)] if refvarstart + len(obs) <= len(mrna_seq) else False:
                canno = f"c.{cds_start+2}dup{obs}"
            else:
                canno = f"c.{cds_start+1}_{cds_start+2}ins{obs}"
        else:
            # Minus strand: "after current site" becomes "before current site"
            if refvarstart >= 2 and obs == mrna_seq[refvarstart - 2:refvarstart - 2 + len(obs)] if refvarstart - 2 + len(obs) <= len(mrna_seq) else False:
                canno = f"c.{cds_start}dup{obs}"
            elif obs == mrna_seq[refvarstart - 1:refvarstart - 1 + len(obs)] if refvarstart - 1 + len(obs) <= len(mrna_seq) else False:
                canno = f"c.{cds_start+1}dup{obs}"
            else:
                canno = f"c.{cds_start}_{cds_start+1}ins{obs}"

        # Protein annotation
        # Build variant codon for translation
        if tx.strand == '+':
            if fs == 1:
                varnt3 = wtnt3[0] + wtnt3[1] + obs + wtnt3[2]
            elif fs == 2:
                varnt3 = wtnt3[0] + wtnt3[1] + wtnt3[2] + obs
            else:
                varnt3 = wtnt3[0] + obs + wtnt3[1] + wtnt3[2]
        else:
            if fs == 1:
                varnt3 = wtnt3[0] + obs + wtnt3[1] + wtnt3[2]
            elif fs == 2:
                varnt3 = wtnt3[0] + wtnt3[1] + obs + wtnt3[2]
            else:
                varnt3 = obs + wtnt3[0] + wtnt3[1] + wtnt3[2]

        varaa = translate_dna(varnt3)

        if ins_len % 3 == 0:
            # Nonframeshift insertion
            if wtaa == 'X':
                # Mutation on stop codon
                if varaa and 'X' in varaa:
                    varaa = varaa[:varaa.index('X')+1]
                    panno = f"p.X{varpos}delins{varaa}"
                else:
                    exonic_func = 'stoploss'
                    panno = f"p.X{varpos}delins{varaa}"
            elif varaa and 'X' in varaa:
                exonic_func = 'stopgain'
                varaa = varaa[:varaa.index('X')+1]
                panno = f"p.{wtaa}{varpos}delins{varaa}"
            else:
                # Normal nonframeshift insertion: compare wt and var protein windows
                # to find the proper p.{prevAA}{prevPos}_{curAA}{curPos}ins{insertedAAs} format
                panno = self._compute_nfs_ins_annotation(mrna_seq, refvarstart, refcdsstart, obs, tx.strand)
        else:
            # Frameshift insertion
            if wtaa == 'X':
                # Mutation on stop codon
                if varaa and 'X' in varaa:
                    varaa = varaa[:varaa.index('X')+1]
                    panno = f"p.X{varpos}delins{varaa}"
                else:
                    exonic_func = 'stoploss'
                    panno = f"p.X{varpos}delins{varaa}"
            elif varaa and 'X' in varaa:
                exonic_func = 'stopgain'
                varaa = varaa[:varaa.index('X')+1]
                # Get wtaa_after
                wt_after_start = (refvarstart - 1) + (3 - fs) if fs > 0 else refvarstart - 1 + 3
                if wt_after_start + 3 <= len(mrna_seq):
                    wtaa_after = translate_dna(mrna_seq[wt_after_start:wt_after_start + 3])
                else:
                    wtaa_after = '?'
                panno = f"p.{wtaa}{varpos}_{wtaa_after}{varpos+1}delins{varaa}"
            else:
                # Compute detailed frameshift: p.{aa}{pos}{mutAA}fs*{dist}
                panno = self._compute_fs_annotation(mrna_seq, refvarstart, refcdsstart, varnt3, fs)

        return (exonic_func, f"{tx.name2}:{tx.name}:exon{exon_pos}:{canno}:{panno}")

    def _annotate_single_deletion(self, tx, mrna_seq, exon_pos, refvarstart, refcdsstart,
                                   cds_start, fs, wtnt3, wtaa, varpos):
        """Annotate a single nucleotide deletion."""
        # Get the next codon for frameshift computation
        wtnt3_after_start = refvarstart - fs + 2
        wtnt3_after = mrna_seq[wtnt3_after_start:wtnt3_after_start + 1] if wtnt3_after_start < len(mrna_seq) else ''

        if fs == 1:
            deletent = wtnt3[1]
            varnt3 = wtnt3[0] + wtnt3[2] + wtnt3_after
        elif fs == 2:
            deletent = wtnt3[2]
            varnt3 = wtnt3[0] + wtnt3[1] + wtnt3_after
        else:
            deletent = wtnt3[0]
            varnt3 = wtnt3[1] + wtnt3[2] + wtnt3_after

        varaa = translate_dna(varnt3)
        canno = f"c.{cds_start+1}del{deletent}"

        if wtaa == 'X':
            # Mutation on stop codon
            if varaa == 'X':
                exonic_func = 'synonymous SNV'  # stop -> stop (nfsdel in ANNOVAR)
                panno = f"p.X{varpos}X"
            else:
                exonic_func = 'stoploss'
                panno = f"p.X{varpos}{varaa}"
        elif varaa == 'X':
            exonic_func = 'stopgain'
            panno = f"p.{wtaa}{varpos}X"
        else:
            exonic_func = 'frameshift deletion'
            # Compute detailed: p.{wtaa}{pos}{mutAA}fs*{dist}
            panno = self._compute_fs_annotation(mrna_seq, refvarstart, refcdsstart, varnt3, fs, is_del=True, refvarend=refvarstart)

        return (exonic_func, f"{tx.name2}:{tx.name}:exon{exon_pos}:{canno}:{panno}")

    def _annotate_deletion(self, tx, mrna_seq, exon_pos, refvarstart, refvarend,
                           refcdsstart, cds_start, cds_end, fs, wtnt3, wtaa, varpos):
        """Annotate a multi-nucleotide deletion."""
        del_len = refvarend - refvarstart + 1
        canno = f"c.{cds_start+1}_{cds_end+1}del"

        if del_len % 3 == 0:
            exonic_func = 'nonframeshift deletion'
            num_del_aa = del_len // 3
            # Compare wt and var protein to find actual deleted AA(s)
            cds_start_0 = refcdsstart - 1  # 0-based mRNA position of CDS start
            # Get a window around the deletion, codon-aligned
            codon_pos_start = cds_start // 3  # 0-based codon number
            window_start_codon = codon_pos_start
            window_end_codon = (cds_end // 3) + 2  # extra codons after
            wt_start = cds_start_0 + window_start_codon * 3
            wt_end = min(len(mrna_seq), cds_start_0 + window_end_codon * 3)
            wt_window = mrna_seq[wt_start:wt_end]
            # Build variant window: remove the deleted bases
            del_offset = cds_start - window_start_codon * 3  # offset within window
            var_window = wt_window[:del_offset] + wt_window[del_offset + del_len:]
            wt_prot = translate_dna(wt_window)
            var_prot = translate_dna(var_window)
            # Find first difference
            first_diff = 0
            for i in range(min(len(wt_prot), len(var_prot))):
                if wt_prot[i] != var_prot[i]:
                    first_diff = i
                    break
            else:
                first_diff = min(len(wt_prot), len(var_prot))
            del_start_pos = window_start_codon + first_diff + 1  # 1-based AA position
            del_start_aa = wt_prot[first_diff] if first_diff < len(wt_prot) else '?'
            if num_del_aa == 1:
                panno = f"p.{del_start_aa}{del_start_pos}del"
            else:
                del_end_pos = del_start_pos + num_del_aa - 1
                del_end_aa = wt_prot[first_diff + num_del_aa - 1] if first_diff + num_del_aa - 1 < len(wt_prot) else '?'
                panno = f"p.{del_start_aa}{del_start_pos}_{del_end_aa}{del_end_pos}del"
        else:
            # Build varnt3 for frameshift deletion
            end_fs = cds_end % 3
            wt_codon_start = refvarstart - fs - 1  # 0-based in mRNA
            # Get the remaining portion after deletion
            after_start = refvarend  # 0-based position right after deleted region
            after_nt = mrna_seq[after_start:after_start + (2 - end_fs)] if after_start + (2 - end_fs) <= len(mrna_seq) else ''
            # Build the variant codon: prefix (before deletion in first codon) + suffix (after deletion)
            prefix = wtnt3[:fs] if fs > 0 else ''
            varnt3 = prefix + after_nt

            # Build full first variant codon (varnt3 may be incomplete)
            wt_codon_start_0 = refvarstart - fs - 1
            full_var_prefix = mrna_seq[wt_codon_start_0:wt_codon_start_0 + fs] if fs > 0 else ''
            full_var_rest = mrna_seq[refvarend:]  # everything after deletion
            full_var_seq = full_var_prefix + full_var_rest
            varaa = translate_dna(full_var_seq[:3]) if len(full_var_seq) >= 3 else ''

            if wtaa == 'X':
                if varaa == 'X':
                    exonic_func = 'synonymous SNV'
                    panno = f"p.X{varpos}X"
                else:
                    exonic_func = 'stoploss'
                    panno = f"p.X{varpos}{varaa}"
            elif varaa == 'X':
                exonic_func = 'stopgain'
                panno = f"p.{wtaa}{varpos}*"
            else:
                exonic_func = 'frameshift deletion'
                panno = self._compute_fs_annotation(mrna_seq, refvarstart, refcdsstart, varnt3, fs, is_del=True, refvarend=refvarend)

        return (exonic_func, f"{tx.name2}:{tx.name}:exon{exon_pos}:{canno}:{panno}")

    def _compute_nfs_ins_annotation(self, mrna_seq, refvarstart, refcdsstart, obs, strand):
        """Compute nonframeshift insertion protein annotation.

        Builds wildtype and variant CDS from CDS start, translates both,
        compares to find the inserted amino acids, and formats as
        p.{prevAA}{prevPos}_{nextAA}{nextPos}ins{insertedAAs}.
        """
        cds_start_0 = refcdsstart - 1  # 0-based mRNA position of CDS start
        cds_pos = refvarstart - refcdsstart  # 0-based position within CDS
        ins_len = len(obs)
        num_new_aa = ins_len // 3

        # Build wildtype and variant CDS windows (a generous window around the insertion)
        # Start a few codons before the insertion point
        codon_pos = cds_pos // 3  # 0-based codon number
        window_start_codon = max(0, codon_pos - 2)
        window_start = cds_start_0 + window_start_codon * 3  # 0-based in mRNA
        window_end = min(len(mrna_seq), window_start + (num_new_aa + 10) * 3)

        wt_window = mrna_seq[window_start:window_end]

        # Build variant window by inserting obs into the mRNA at the right position
        # For + strand: insertion is AFTER refvarstart
        # For - strand: insertion is BEFORE refvarstart (because "after" in genomic = "before" in cDNA)
        ins_pos_in_window = (cds_start_0 + cds_pos) - window_start  # position in window
        if strand == '-':
            # Insert before the current position
            var_window = wt_window[:ins_pos_in_window] + obs + wt_window[ins_pos_in_window:]
        else:
            # Insert after the current position
            var_window = wt_window[:ins_pos_in_window + 1] + obs + wt_window[ins_pos_in_window + 1:]

        wt_protein = translate_dna(wt_window)
        var_protein = translate_dna(var_window)

        if not wt_protein or not var_protein:
            varpos = cds_pos // 3 + 1
            return f"p.?{varpos}?"

        # Compare wt and var proteins to find inserted AAs
        # The var protein should be num_new_aa longer than wt protein
        # Find the first position that differs
        first_diff = None
        for i in range(len(wt_protein)):
            if i >= len(var_protein) or wt_protein[i] != var_protein[i]:
                first_diff = i
                break

        if first_diff is None:
            # No difference found in the compared region - synonymous-like
            first_diff = len(wt_protein)

        # The inserted AAs are var_protein[first_diff:first_diff+num_new_aa]
        inserted_aas = var_protein[first_diff:first_diff + num_new_aa]

        # The flanking AAs in the wildtype
        prev_pos_in_wt = first_diff - 1  # position in wt_protein (0-based)
        next_pos_in_wt = first_diff  # position in wt_protein (0-based)

        if prev_pos_in_wt >= 0 and next_pos_in_wt < len(wt_protein):
            prev_aa = wt_protein[prev_pos_in_wt]
            next_aa = wt_protein[next_pos_in_wt]
            # Convert to 1-based AA positions
            prev_aa_pos = window_start_codon + prev_pos_in_wt + 1
            next_aa_pos = window_start_codon + next_pos_in_wt + 1
            return f"p.{prev_aa}{prev_aa_pos}_{next_aa}{next_aa_pos}ins{inserted_aas}"
        else:
            varpos = cds_pos // 3 + 1
            return f"p.{wt_protein[0] if wt_protein else '?'}{varpos}delins{var_protein[:num_new_aa+1] if var_protein else '?'}"

    def _compute_fs_annotation(self, mrna_seq, refvarstart, refcdsstart, varnt3, fs,
                               is_del=False, refvarend=None):
        """Compute detailed frameshift annotation: p.{wtAA}{pos}{mutAA}fs*{dist}.

        Builds the full mutant mRNA by splicing the variant into the wildtype,
        translates both from the affected codon, finds the first different AA
        and the distance to the next stop codon.
        """
        wt_codon_start = refvarstart - fs - 1  # 0-based position in mRNA
        varpos = (refvarstart - refcdsstart) // 3 + 1

        # Build wildtype CDS from affected codon to end of mRNA
        max_len = min(len(mrna_seq), wt_codon_start + 3000)
        wt_seq = mrna_seq[wt_codon_start:max_len]
        wt_protein = translate_dna(wt_seq)

        # Build mutant mRNA by splicing
        # prefix: wt bases before variant in the affected codon
        # varnt3 = prefix + inserted/modified bases + suffix from wt
        if is_del:
            # For deletion: the mutant mRNA skips the deleted bases
            # varnt3 was built from prefix (fs bases) + bases from after deletion
            # The continuation is the rest of the mRNA after the deletion
            if refvarend is None:
                refvarend = refvarstart  # single deletion
            # Mutant mRNA from this codon = prefix + (mRNA from refvarend onward)
            prefix = mrna_seq[wt_codon_start:wt_codon_start + fs] if fs > 0 else ''
            var_full = prefix + mrna_seq[refvarend:max_len]
        else:
            # For insertion: varnt3 = modified codon with insertion
            # The continuation is the rest of the mRNA after the original codon
            var_full = varnt3 + mrna_seq[wt_codon_start + 3:max_len]

        var_protein = translate_dna(var_full)

        if not var_protein:
            return f"p.{wt_protein[0] if wt_protein else '?'}{varpos}fs"

        # Find first mismatch position
        first_diff = 0
        min_len = min(len(wt_protein), len(var_protein))
        for i in range(min_len):
            if wt_protein[i] != var_protein[i]:
                first_diff = i
                break
        else:
            first_diff = min_len

        if first_diff >= len(wt_protein) or first_diff >= len(var_protein):
            return f"p.{wt_protein[0] if wt_protein else '?'}{varpos}fs"

        diff_wt_aa = wt_protein[first_diff]
        diff_var_aa = var_protein[first_diff]
        diff_pos = varpos + first_diff

        # Find distance to next stop codon in variant protein (from first_diff onward)
        stop_dist = None
        for i in range(first_diff, len(var_protein)):
            if var_protein[i] == 'X':
                stop_dist = i - first_diff + 1
                break

        if stop_dist is not None:
            return f"p.{diff_wt_aa}{diff_pos}{diff_var_aa}fs*{stop_dist}"
        else:
            return f"p.{diff_wt_aa}{diff_pos}{diff_var_aa}fs"

    def _compute_block_sub_panno(self, mrna_seq, refvarstart, refvarend, refcdsstart, obs, fs, cds_start):
        """Compute protein annotation for block substitutions.

        Builds wildtype and variant CDS windows, translates both, and compares
        to find the protein change. Returns the protein annotation string
        (e.g., ':p.E201E' or ':p.R106_G107delinsLR').
        """
        ref_len = refvarend - refvarstart + 1
        is_nfs = (ref_len - len(obs)) % 3 == 0
        cds_start_0 = refcdsstart - 1  # 0-based mRNA pos of CDS start

        if not is_nfs:
            return ''

        # Build wildtype and variant CDS windows, aligned to codon boundaries
        cds_pos = refvarstart - refcdsstart  # 0-based CDS position of variant start
        codon_start = (cds_pos // 3) * 3  # align to codon boundary (0-based in CDS)
        cds_end_pos = cds_pos + ref_len  # 0-based CDS end (exclusive)
        codon_end = ((cds_end_pos + 2) // 3) * 3  # align to next codon boundary

        wt_start = cds_start_0 + codon_start  # 0-based in mRNA
        wt_end = cds_start_0 + codon_end

        if wt_start < 0 or wt_end > len(mrna_seq):
            return ''

        wt_seq = mrna_seq[wt_start:wt_end]
        # Build variant: replace the affected bases
        offset_in_window = cds_pos - codon_start  # offset from window start to variant
        var_seq = wt_seq[:offset_in_window] + obs + wt_seq[offset_in_window + ref_len:]

        wt_aa = translate_dna(wt_seq)
        var_aa = translate_dna(var_seq)

        if not wt_aa or not var_aa:
            return ''

        varpos = codon_start // 3 + 1  # 1-based AA position

        if wt_aa == var_aa:
            # Synonymous
            return f":p.{wt_aa[0]}{varpos}{var_aa[0]}"
        elif len(wt_aa) == 1 and len(var_aa) == 1:
            return f":p.{wt_aa}{varpos}{var_aa}"
        else:
            # Find the range of affected amino acids
            first_diff = 0
            for i in range(min(len(wt_aa), len(var_aa))):
                if wt_aa[i] != var_aa[i]:
                    first_diff = i
                    break
            last_diff = min(len(wt_aa), len(var_aa)) - 1
            for i in range(min(len(wt_aa), len(var_aa)) - 1, first_diff - 1, -1):
                if wt_aa[i] != var_aa[i]:
                    last_diff = i
                    break

            start_pos = varpos + first_diff
            end_pos = varpos + last_diff
            wt_range = wt_aa[first_diff:last_diff + 1]
            var_range = var_aa[first_diff:last_diff + 1]

            if start_pos == end_pos:
                return f":p.{wt_range}{start_pos}{var_range}"
            else:
                return f":p.{wt_range[0]}{start_pos}_{wt_range[-1]}{end_pos}delins{var_range}"

    def _get_aa_at(self, mrna_seq: str, refcdsstart: int, cds_pos: int) -> str:
        """Get amino acid at a given 0-based CDS position."""
        codon_num = cds_pos // 3
        codon_start = (refcdsstart - 1) + codon_num * 3  # 0-based in mRNA
        if codon_start + 3 > len(mrna_seq) or codon_start < 0:
            return '?'
        codon = mrna_seq[codon_start:codon_start + 3]
        return translate_dna(codon) if len(codon) == 3 else '?'

    def _calc_mrna_positions_plus(self, tx: Transcript, start: int, end: int) -> Tuple[int, int, int]:
        """Calculate mRNA positions for + strand, matching ANNOVAR's exact algorithm.
        Returns (refvarstart, refvarend, refcdsstart) - all 1-based in mRNA.
        """
        lenintron = 0
        refvarstart = 0
        refvarend = 0
        refcdsstart = 0

        for k in range(tx.exoncount):
            es = tx.exonstarts[k]
            ee = tx.exonends[k]
            if k > 0:
                lenintron += es - tx.exonends[k-1] - 1

            if tx.cdsstart >= es:
                refcdsstart = tx.cdsstart - tx.txstart - lenintron + 1

            if start >= es and start <= ee:
                refvarstart = start - tx.txstart - lenintron + 1
            if end >= es and end <= ee:
                refvarend = end - tx.txstart - lenintron + 1

        return (refvarstart, refvarend, refcdsstart)

    def _calc_mrna_positions_minus(self, tx: Transcript, start: int, end: int) -> Tuple[int, int, int]:
        """Calculate mRNA positions for - strand, matching ANNOVAR's exact algorithm.
        Returns (refvarstart, refvarend, refcdsstart) - all 1-based in mRNA.
        """
        lenintron = 0
        refvarstart = 0
        refvarend = 0
        refcdsstart = 0

        for k in range(tx.exoncount - 1, -1, -1):
            es = tx.exonstarts[k]
            ee = tx.exonends[k]
            if k < tx.exoncount - 1:
                lenintron += tx.exonstarts[k+1] - ee - 1

            if tx.cdsend >= es and tx.cdsend <= ee:
                refcdsstart = tx.txend - tx.cdsend - lenintron + 1

            if end >= es and end <= ee:
                refvarstart = tx.txend - end - lenintron + 1
            if start >= es and start <= ee:
                refvarend = tx.txend - start - lenintron + 1

        return (refvarstart, refvarend, refcdsstart)

    def _sort_exonic_annotation(self, aa_change: str) -> str:
        """Sort AA change entries by exon number then transcript name.
        Matches ANNOVAR's sortExonicAnnotation function.
        """
        parts = [p for p in aa_change.split(',') if p]
        sortable = []
        for part in parts:
            fields = part.split(':')
            if len(fields) >= 3:
                exon_str = fields[2]
                exon_num = 0
                m = re.match(r'exon(\d+)', exon_str)
                if m:
                    exon_num = int(m.group(1))
                tx_name = fields[1] if len(fields) > 1 else ''
                sortable.append((exon_num, tx_name, part))
            else:
                sortable.append((999999, '', part))
        sortable.sort(key=lambda x: (x[0], x[1]))
        return ','.join(s[2] for s in sortable)

    def _format_result(self, func, gene, gene_detail, exonic_func, aa_change) -> Dict[str, str]:
        """Format the result dictionary matching the column names."""
        suffix = self.db_name  # 'refGene' or 'ensGene'
        return {
            f'Func.{suffix}': func,
            f'Gene.{suffix}': gene,
            f'GeneDetail.{suffix}': gene_detail,
            f'ExonicFunc.{suffix}': exonic_func,
            f'AAChange.{suffix}': aa_change,
        }

    def reopen_tabix(self):
        """Re-open tabix file handle (required after fork to avoid sharing file descriptors)."""
        config = DATABASE_CONFIG[self.db_name]
        gene_file = os.path.join(HUMANDB_DIR, config['file'])
        if self.tabix:
            self.tabix.close()
        if os.path.exists(gene_file) and os.path.exists(gene_file + '.tbi'):
            self.tabix = pysam.TabixFile(gene_file)

    def close(self):
        if self.tabix:
            self.tabix.close()
