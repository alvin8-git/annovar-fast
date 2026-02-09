#!/usr/bin/env python3
"""
Annotation engine for processing variants with all databases.
Orchestrates gene annotation and database queries, writes output.
"""

import logging
import time
from typing import Dict, List, Any

from vcf_processor import Variant, VCFProcessor
from database_manager import DatabaseManager
from gene_annotator import GeneAnnotator
from config import DATABASE_ORDER, DATABASE_CONFIG, OUTPUT_COLUMNS, NA_STRING

logger = logging.getLogger(__name__)


# Pre-build the mapping from OUTPUT_COLUMNS index to (db_name, col_index_within_db).
# This handles duplicate column names (e.g., gnomad211_exome and gnomad30_genome both have 'AF').
def _build_column_map():
    """Build a list of (db_name, col_idx) for each annotation column in OUTPUT_COLUMNS."""
    # OUTPUT_COLUMNS[0:5] = Chr,Start,End,Ref,Alt
    # OUTPUT_COLUMNS[-11:] = Otherinfo1..Otherinfo11
    # Everything in between maps to database columns, in DATABASE_ORDER order.
    col_map = []
    col_offset = 5  # skip first 5 variant columns

    for db_name in DATABASE_ORDER:
        config = DATABASE_CONFIG[db_name]
        columns = config['columns']
        for i, col in enumerate(columns):
            col_map.append((db_name, i))

    return col_map

_COLUMN_MAP = _build_column_map()


class AnnotationEngine:
    """Annotation engine that processes variants sequentially."""

    def __init__(self):
        self.db_manager = DatabaseManager()
        self.gene_annotators = {}
        for db_name in ('refGene', 'ensGene'):
            if db_name in DATABASE_CONFIG:
                self.gene_annotators[db_name] = GeneAnnotator(db_name)
        logger.info("Annotation engine initialized")

    def process_vcf(self, vcf_path: str, output_path: str):
        """Process entire VCF file and write annotations."""
        start_time = time.time()
        vcf_processor = VCFProcessor(vcf_path)

        # Build header: annotation columns + dynamic Otherinfo columns
        # Otherinfo1 = zygosity, Otherinfo2-10 = VCF fixed fields (CHROM..FORMAT),
        # Otherinfo11+ = one per sample
        num_samples = vcf_processor.num_samples
        num_otherinfo = 1 + 9 + num_samples  # zygosity + 9 VCF fixed fields + samples
        annotation_cols = [c for c in OUTPUT_COLUMNS if not c.startswith('Otherinfo')]
        otherinfo_cols = [f'Otherinfo{i}' for i in range(1, num_otherinfo + 1)]
        header = annotation_cols + otherinfo_cols

        with open(output_path, 'w') as outfile:
            outfile.write('\t'.join(header) + '\n')
            variant_count = 0
            for variant in vcf_processor.get_variants():
                line = self._annotate_and_format(variant)
                outfile.write(line + '\n')
                variant_count += 1

        elapsed = time.time() - start_time
        logger.info(f"Completed {variant_count} variants in {elapsed:.1f}s")

    def _annotate_and_format(self, variant: Variant) -> str:
        """Annotate a variant and format as output line."""
        # Query each database and store results keyed by db_name
        db_results = {}  # db_name -> list of values (matching columns order)

        for db_name in DATABASE_ORDER:
            config = DATABASE_CONFIG[db_name]
            db_type = config['type']
            columns = config['columns']

            if db_type == 'gene':
                if db_name in self.gene_annotators:
                    result = self.gene_annotators[db_name].annotate(
                        variant.chrom, variant.start, variant.end,
                        variant.ref, variant.alt
                    )
                else:
                    result = {col: NA_STRING for col in columns}
                # Convert dict to ordered list
                db_results[db_name] = [result.get(col, NA_STRING) for col in columns]
            else:
                result = self.db_manager.query_variant(
                    db_name, variant.chrom, variant.start, variant.end,
                    variant.ref, variant.alt
                )
                db_results[db_name] = [result.get(col, NA_STRING) for col in columns]

        # Build output fields
        fields = [
            variant.chrom,
            str(variant.start),
            str(variant.end),
            variant.ref,
            variant.alt,
        ]

        # Add annotation columns in order using the column map
        for db_name, col_idx in _COLUMN_MAP:
            vals = db_results.get(db_name)
            if vals and col_idx < len(vals):
                fields.append(vals[col_idx])
            else:
                fields.append(NA_STRING)

        # Add Otherinfo columns
        fields.append(variant.zygosity)
        fields.extend(variant.vcf_fields)

        return '\t'.join(fields)

    def close(self):
        self.db_manager.close()
        for ga in self.gene_annotators.values():
            ga.close()
