#!/usr/bin/env python3
"""
VCF processor using cyvcf2 for C-speed parsing.
Converts VCF variants to ANNOVAR format matching convert2annovar.pl --format vcf4old --includeinfo --withzyg
"""

import os
import sys
import logging
from collections import defaultdict
from typing import List, Dict, Any, Iterator, Tuple

try:
    import cyvcf2
except ImportError:
    print("Error: cyvcf2 not installed. Please install with: pip install cyvcf2")
    sys.exit(1)

logger = logging.getLogger(__name__)


class Variant:
    """Represents a genetic variant in ANNOVAR format"""
    __slots__ = ['chrom', 'start', 'end', 'ref', 'alt', 'zygosity', 'vcf_fields']

    def __init__(self, chrom, start, end, ref, alt, zygosity, vcf_fields):
        self.chrom = chrom
        self.start = start
        self.end = end
        self.ref = ref
        self.alt = alt
        self.zygosity = zygosity
        self.vcf_fields = vcf_fields  # list: [CHROM, POS, ID, REF, ALT, QUAL, FILTER, INFO, FORMAT, SAMPLE1, SAMPLE2, ...]

    def __str__(self):
        return f"{self.chrom}:{self.start}:{self.ref}>{self.alt}"


class VCFProcessor:
    """VCF processor matching convert2annovar.pl --format vcf4old --includeinfo --withzyg"""

    def __init__(self, vcf_path: str):
        self.vcf_path = vcf_path
        self.vcf = None
        self.sample_names = []
        self.num_samples = 0
        self._raw_lines = []
        self._load_vcf()

    def _load_vcf(self):
        """Load VCF file with cyvcf2 and also read raw lines for exact field preservation."""
        self.vcf = cyvcf2.VCF(self.vcf_path)
        self.sample_names = self.vcf.samples
        self.num_samples = len(self.sample_names)

        # Read raw data lines (non-header) for exact field preservation
        import gzip
        opener = gzip.open if self.vcf_path.endswith('.gz') else open
        with opener(self.vcf_path, 'rt') as f:
            for line in f:
                if not line.startswith('#'):
                    self._raw_lines.append(line.rstrip('\n\r'))

        logger.info(f"Loaded VCF: {self.vcf_path} with {self.num_samples} samples")

    def get_variants(self) -> Iterator[Variant]:
        """Iterate over all variants, converting to ANNOVAR format"""
        if self.vcf is None:
            raise RuntimeError("VCF not loaded")

        for line_idx, variant in enumerate(self.vcf):
            alts = variant.ALT
            if not alts:
                continue

            raw_line = self._raw_lines[line_idx] if line_idx < len(self._raw_lines) else None

            # Only process first ALT allele (matches convert2annovar.pl --format vcf4old)
            alt = alts[0]
            if alt == '<M>' or alt == '*':
                continue
            yield self._convert_to_annovar_format(variant, alt, raw_line)

    def _convert_to_annovar_format(self, cyvcf_variant, alt: str, raw_line: str) -> Variant:
        """Convert cyvcf2 variant to ANNOVAR format.

        Matches convert2annovar.pl --format vcf4old logic:
        - SNP: start=pos, end=pos, ref=REF, alt=ALT
        - Deletion: strip common prefix, start=pos+len(prefix), end=start+len(deleted)-1, ref=deleted, alt="-"
        - Insertion: strip common prefix, start=pos+len(prefix)-1, end=start, ref="-", alt=inserted
        """
        chrom = cyvcf_variant.CHROM
        pos = int(cyvcf_variant.POS)
        ref = cyvcf_variant.REF

        # Compute zygosity from GT
        zygosity = self._compute_zygosity(cyvcf_variant, alt)

        # Build VCF fields from raw line for exact formatting preservation
        vcf_fields = self._build_vcf_fields(raw_line)

        # Convert to ANNOVAR coordinate format
        if len(ref) == 1 and len(alt) == 1:
            # SNP
            annovar_start = pos
            annovar_end = pos
            annovar_ref = ref
            annovar_alt = alt
        elif len(ref) > len(alt):
            # Deletion (or complex where ref is longer)
            prefix_len = 0
            min_len = min(len(ref), len(alt))
            while prefix_len < min_len and ref[prefix_len] == alt[prefix_len]:
                prefix_len += 1

            if prefix_len > 0 and prefix_len == len(alt):
                # Pure deletion
                deleted = ref[prefix_len:]
                annovar_start = pos + prefix_len
                annovar_end = annovar_start + len(deleted) - 1
                annovar_ref = deleted
                annovar_alt = '-'
            else:
                # Block substitution
                annovar_start = pos
                annovar_end = pos + len(ref) - 1
                annovar_ref = ref
                annovar_alt = alt
        elif len(alt) > len(ref):
            # Insertion (or complex where alt is longer)
            prefix_len = 0
            min_len = min(len(ref), len(alt))
            while prefix_len < min_len and ref[prefix_len] == alt[prefix_len]:
                prefix_len += 1

            if prefix_len > 0 and prefix_len == len(ref):
                # Pure insertion
                inserted = alt[prefix_len:]
                annovar_start = pos + prefix_len - 1
                annovar_end = annovar_start
                annovar_ref = '-'
                annovar_alt = inserted
            else:
                # Block substitution
                annovar_start = pos
                annovar_end = pos + len(ref) - 1
                annovar_ref = ref
                annovar_alt = alt
        else:
            # Same length, multi-base substitution
            annovar_start = pos
            annovar_end = pos + len(ref) - 1
            annovar_ref = ref
            annovar_alt = alt

        return Variant(
            chrom=chrom,
            start=annovar_start,
            end=annovar_end,
            ref=annovar_ref,
            alt=annovar_alt,
            zygosity=zygosity,
            vcf_fields=vcf_fields,
        )

    def _compute_zygosity(self, cyvcf_variant, alt: str) -> str:
        """Compute zygosity: 'hom' if 1/1, 'het' if 0/1, 'unknown' if missing."""
        try:
            gt = cyvcf_variant.genotypes[0]
            a1, a2 = gt[0], gt[1]
            if a1 < 0 or a2 < 0:
                return 'unknown'
            if a1 == a2 and a1 > 0:
                return 'hom'
            elif a1 != a2:
                return 'het'
            else:
                return 'hom'  # 0/0 shouldn't happen for called variants
        except Exception:
            return 'unknown'

    def get_variants_by_chrom(self) -> Tuple[Dict[str, List[Tuple[int, 'Variant']]], int]:
        """Group all variants by chromosome, preserving original order.

        Returns:
            (chrom_dict, total_count) where chrom_dict maps chrom -> list of (original_index, Variant)
        """
        by_chrom = defaultdict(list)
        idx = 0
        for variant in self.get_variants():
            by_chrom[variant.chrom].append((idx, variant))
            idx += 1
        return dict(by_chrom), idx

    def _build_vcf_fields(self, raw_line: str) -> list:
        """Build VCF fields for Otherinfo columns from the raw VCF line.

        Returns [CHROM, POS, ID, REF, ALT, QUAL, FILTER, INFO, FORMAT, SAMPLE1, SAMPLE2, ...]
        matching what convert2annovar.pl --includeinfo produces.
        Uses raw line to preserve exact formatting (e.g., QUAL as 68.00 not 68).
        """
        if raw_line:
            fields = raw_line.split('\t')
            # Return all fields: fixed VCF cols + all sample cols
            return fields
        return ['.'] * 10
