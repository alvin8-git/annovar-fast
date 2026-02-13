#!/usr/bin/env python3
"""
ANNOVAR-fast: High-performance variant annotation using tabix-indexed databases.
Drop-in replacement for convert2annovar.pl + table_annovar.pl pipeline.
"""

import os
import sys
import argparse
import logging
import time
from pathlib import Path

import pysam
# Suppress benign htslib warnings (e.g., "Coordinate <= 0", "index older than data")
# Level 1 = errors only; warnings are noise from database records with 0-based coords
pysam.set_verbosity(1)

# Setup logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)


def parse_arguments():
    parser = argparse.ArgumentParser(
        description='ANNOVAR-fast: High-performance variant annotation using tabix-indexed databases'
    )
    parser.add_argument('vcf_file', help='Input VCF file to annotate')
    parser.add_argument('-o', '--output', help='Output file name', default=None)
    parser.add_argument('-v', '--verbose', action='store_true', help='Enable verbose logging')
    return parser.parse_args()


def main():
    args = parse_arguments()

    if args.verbose:
        logging.getLogger().setLevel(logging.DEBUG)

    vcf_file = args.vcf_file
    if not os.path.exists(vcf_file):
        logger.error(f"Input VCF file does not exist: {vcf_file}")
        sys.exit(1)

    # Determine output path
    if args.output:
        output_path = args.output
    else:
        sample_name = Path(vcf_file).stem
        output_path = f"{sample_name}.annovar.txt"

    logger.info(f"Input: {vcf_file}")
    logger.info(f"Output: {output_path}")

    from annotation_engine import AnnotationEngine

    engine = AnnotationEngine()
    try:
        engine.process_vcf(vcf_file, output_path)
    finally:
        engine.close()

    # Compare with reference if available
    reference_path = vcf_file.replace('.vcf', '.annovar.txt')
    if os.path.exists(reference_path) and output_path != reference_path:
        _compare(output_path, reference_path)

    logger.info("Done!")


def _compare(output_path, reference_path):
    """Compare output with reference."""
    with open(output_path) as f1, open(reference_path) as f2:
        out_lines = f1.readlines()
        ref_lines = f2.readlines()

    if len(out_lines) != len(ref_lines):
        logger.warning(f"Line count mismatch: {len(out_lines)} vs {len(ref_lines)}")
        return

    mismatches = 0
    for i, (ol, rl) in enumerate(zip(out_lines, ref_lines)):
        if ol != rl:
            mismatches += 1
            if mismatches <= 3:
                # Show first differing column
                ocols = ol.strip().split('\t')
                rcols = rl.strip().split('\t')
                for j, (oc, rc) in enumerate(zip(ocols, rcols)):
                    if oc != rc:
                        col_name = f"col{j}"
                        if j < len(['Chr', 'Start', 'End', 'Ref', 'Alt']):
                            pass
                        logger.warning(f"Line {i+1} col {j}: got '{oc[:80]}' expected '{rc[:80]}'")
                        break

    if mismatches == 0:
        logger.info("Output matches reference perfectly!")
    else:
        logger.warning(f"{mismatches} lines differ from reference")


if __name__ == "__main__":
    main()
