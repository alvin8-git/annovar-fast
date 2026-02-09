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

from config import HUMANDB_TBI_DIR, HUMANDB_DIR, NA_STRING, DATABASE_CONFIG

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
        gene_file = os.path.join(HUMANDB_TBI_DIR, config['file'])
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
        exonic_func = '.'
        aa_change = '.'
        exonic_found = False
        splicing_found = False
        intronic_found = False

        for tx in transcripts:
            func, detail = self._classify_variant(tx, start, end, ref, alt)

            if func == 'exonic':
                exonic_found = True
                gene_names.add(tx.name2)
                # Compute AA change
                ef, ac = self._compute_exonic_annotation(tx, start, end, ref, alt)
                if ef != '.':
                    exonic_func = ef
                if ac != '.':
                    if aa_change == '.':
                        aa_change = ac
                    else:
                        aa_change += ',' + ac
            elif func == 'splicing':
                splicing_found = True
                gene_names.add(tx.name2)
                if detail:
                    gene_detail_parts.append(detail)
            elif func == 'intronic':
                intronic_found = True
                gene_names.add(tx.name2)
            elif func in ('UTR5', 'UTR3'):
                gene_names.add(tx.name2)
                if best_func is None or best_func in ('intronic',):
                    best_func = func
            elif func == 'ncRNA_exonic':
                gene_names.add(tx.name2)
                if best_func is None:
                    best_func = func

        # Determine final function based on priority
        if exonic_found:
            final_func = 'exonic'
        elif splicing_found:
            final_func = 'splicing'
        elif best_func:
            final_func = best_func
        elif intronic_found:
            final_func = 'intronic'
        else:
            final_func = 'intergenic'

        gene_name = ','.join(sorted(gene_names)) if gene_names else '.'
        gene_detail = ';'.join(gene_detail_parts) if gene_detail_parts else '.'

        if final_func == 'splicing':
            exonic_func = '.'
            aa_change = '.'
        elif final_func == 'intronic':
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
        # Check if variant is within transcript boundaries
        if end < tx.txstart or start > tx.txend:
            return ('intergenic', '')

        # Check if this is a non-coding transcript (cdsstart == cdsend)
        is_ncrna = (tx.cdsstart >= tx.cdsend)

        # Iterate through exons to classify
        for k in range(tx.exoncount):
            exon_start = tx.exonstarts[k]
            exon_end = tx.exonends[k]

            # Check if variant overlaps this exon
            if start >= exon_start and end <= exon_end:
                # Variant is within this exon
                if is_ncrna:
                    return ('ncRNA_exonic', '')

                # Check if in CDS or UTR
                if start >= tx.cdsstart and end <= tx.cdsend:
                    return ('exonic', '')
                elif end < tx.cdsstart:
                    if tx.strand == '+':
                        return ('UTR5', '')
                    else:
                        return ('UTR3', '')
                elif start > tx.cdsend:
                    if tx.strand == '+':
                        return ('UTR3', '')
                    else:
                        return ('UTR5', '')
                else:
                    # Spans CDS boundary - treat as exonic
                    return ('exonic', '')

            # Check if variant partially overlaps exon (for indels)
            if start <= exon_end and end >= exon_start:
                if is_ncrna:
                    return ('ncRNA_exonic', '')
                # Check CDS overlap
                if start <= tx.cdsend and end >= tx.cdsstart:
                    return ('exonic', '')
                elif end < tx.cdsstart:
                    return ('UTR5' if tx.strand == '+' else 'UTR3', '')
                elif start > tx.cdsend:
                    return ('UTR3' if tx.strand == '+' else 'UTR5', '')
                else:
                    return ('exonic', '')

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
            return ('splicing', ';'.join(splicing_details))

        # If we're inside the transcript but not in any exon and not near splice sites
        if start >= tx.txstart and end <= tx.txend:
            return ('intronic', '')

        return ('intergenic', '')

    def _has_overlap(self, start1, end1, start2, end2):
        """Check if two ranges overlap."""
        return start1 <= end2 and end1 >= start2

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
                return (exonic_func, f"{tx.name2}:{tx.name}:exon{exon_pos}:{canno}")
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
                return (exonic_func, f"{tx.name2}:{tx.name}:exon{exon_pos}:{canno}")

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
            if varaa and '*' in varaa:
                varaa = varaa[:varaa.index('*')] + 'X'
            panno = f"p.{wtaa}{varpos}delins{varaa}"
        else:
            # Frameshift insertion
            panno = f"p.{wtaa}{varpos}fs"

        return (exonic_func, f"{tx.name2}:{tx.name}:exon{exon_pos}:{canno}:{panno}")

    def _annotate_single_deletion(self, tx, mrna_seq, exon_pos, refvarstart, refcdsstart,
                                   cds_start, fs, wtnt3, wtaa, varpos):
        """Annotate a single nucleotide deletion."""
        # Get the next codon for frameshift computation
        wt_codon_start = refvarstart - fs - 1
        wtnt3_after_start = refvarstart - fs + 2
        wtnt3_after = mrna_seq[wtnt3_after_start:wtnt3_after_start + 3] if wtnt3_after_start + 3 <= len(mrna_seq) else ''

        if fs == 1:
            deletent = wtnt3[1]
            varnt3 = wtnt3[0] + wtnt3[2] + (wtnt3_after[0] if wtnt3_after else '')
        elif fs == 2:
            deletent = wtnt3[2]
            varnt3 = wtnt3[0] + wtnt3[1] + (wtnt3_after[0] if wtnt3_after else '')
        else:
            deletent = wtnt3[0]
            varnt3 = wtnt3[1] + wtnt3[2] + (wtnt3_after[0] if wtnt3_after else '')

        canno = f"c.{cds_start+1}del"
        # Single nt deletion is always frameshift
        exonic_func = 'frameshift deletion'
        panno = f"p.{wtaa}{varpos}fs"
        return (exonic_func, f"{tx.name2}:{tx.name}:exon{exon_pos}:{canno}:{panno}")

    def _annotate_deletion(self, tx, mrna_seq, exon_pos, refvarstart, refvarend,
                           refcdsstart, cds_start, cds_end, fs, wtnt3, wtaa, varpos):
        """Annotate a multi-nucleotide deletion."""
        del_len = refvarend - refvarstart + 1
        canno = f"c.{cds_start+1}_{cds_end+1}del"

        if del_len % 3 == 0:
            exonic_func = 'nonframeshift deletion'
            if fs > 0:
                adj_varpos = (cds_start - fs + 3) // 3 + 1
            else:
                adj_varpos = varpos
            varposend = cds_end // 3 + 1
            start_aa = self._get_aa_at(mrna_seq, refcdsstart, (adj_varpos - 1) * 3)
            end_aa = self._get_aa_at(mrna_seq, refcdsstart, (varposend - 1) * 3)
            panno = f"p.{start_aa}{adj_varpos}_{end_aa}{varposend}del"
        else:
            exonic_func = 'frameshift deletion'
            panno = f"p.{wtaa}{varpos}fs"

        return (exonic_func, f"{tx.name2}:{tx.name}:exon{exon_pos}:{canno}:{panno}")

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

    def close(self):
        if self.tabix:
            self.tabix.close()
