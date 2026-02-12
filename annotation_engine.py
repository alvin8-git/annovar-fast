#!/usr/bin/env python3
"""
Annotation engine for processing variants with all databases.
Orchestrates gene annotation and database queries, writes output.
Supports parallel annotation via fork-based chromosome chunking.
"""

import logging
import multiprocessing
import time
from typing import Dict, List, Any, Tuple

from vcf_processor import Variant, VCFProcessor
from database_manager import DatabaseManager
from gene_annotator import GeneAnnotator
from config import DATABASE_ORDER, DATABASE_CONFIG, OUTPUT_COLUMNS, NA_STRING, MAX_WORKERS

logger = logging.getLogger(__name__)


# Pre-build the mapping from OUTPUT_COLUMNS index to (db_name, col_index_within_db).
# This handles the mapping from database columns to output positions.
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


# -------------------------------------------------------
# Module-level worker state and functions for multiprocessing
# -------------------------------------------------------

# Set by the parent process before forking so workers inherit it via COW.
# Workers then re-open tabix file handles in _init_worker.
_worker_engine = None


def _init_worker():
    """Initialize worker process after fork: re-open tabix file handles.

    The _worker_engine global is inherited from the parent via fork (COW).
    We only need to re-open tabix file descriptors which can't be shared.
    """
    _worker_engine.db_manager.reopen_all_tabix()
    for ga in _worker_engine.gene_annotators.values():
        ga.reopen_tabix()


def _annotate_chunk(chunk):
    """Annotate a list of (index, variant) tuples. Returns list of (index, line).

    Runs in a forked worker process with _worker_engine set by _init_worker.
    """
    results = []
    for idx, variant in chunk:
        line = _worker_engine._annotate_and_format(variant)
        results.append((idx, line))
    return results


class AnnotationEngine:
    """Annotation engine that processes variants, with optional parallelism."""

    def __init__(self):
        self.db_manager = DatabaseManager()
        self.gene_annotators = {}
        for db_name in ('refGene', 'ensGene'):
            if db_name in DATABASE_CONFIG:
                self.gene_annotators[db_name] = GeneAnnotator(db_name)
        logger.info("Annotation engine initialized")

    def process_vcf(self, vcf_path: str, output_path: str, parallel: bool = True):
        """Process entire VCF file and write annotations.

        Args:
            vcf_path: Path to input VCF file.
            output_path: Path to output annotation file.
            parallel: Use fork-based parallelism (default True).
        """
        start_time = time.time()
        vcf_processor = VCFProcessor(vcf_path)

        # Build header: annotation columns + dynamic Otherinfo columns
        num_samples = vcf_processor.num_samples
        num_otherinfo = 1 + 9 + num_samples
        annotation_cols = [c for c in OUTPUT_COLUMNS if not c.startswith('Otherinfo')]
        otherinfo_cols = [f'Otherinfo{i}' for i in range(1, num_otherinfo + 1)]
        header = annotation_cols + otherinfo_cols

        if parallel:
            variant_count = self._process_parallel(vcf_processor, output_path, header)
        else:
            variant_count = self._process_sequential(vcf_processor, output_path, header)

        elapsed = time.time() - start_time
        logger.info(f"Completed {variant_count} variants in {elapsed:.1f}s")

    def _process_sequential(self, vcf_processor, output_path, header):
        """Original sequential processing."""
        with open(output_path, 'w') as outfile:
            outfile.write('\t'.join(header) + '\n')
            variant_count = 0
            for variant in vcf_processor.get_variants():
                line = self._annotate_and_format(variant)
                outfile.write(line + '\n')
                variant_count += 1
        return variant_count

    def _process_parallel(self, vcf_processor, output_path, header):
        """Fork-based parallel processing with equal-size chunking.

        Strategy:
        - Read all variants, assign sequential indices
        - Split into MAX_WORKERS equal-size chunks (not by chromosome)
        - Fork workers that inherit mRNA data via copy-on-write
        - Each worker re-opens tabix handles and annotates its chunk
        - Merge results in original order
        """
        # Collect all variants with indices
        all_variants = []
        for variant in vcf_processor.get_variants():
            all_variants.append((len(all_variants), variant))

        total_count = len(all_variants)
        if total_count == 0:
            with open(output_path, 'w') as outfile:
                outfile.write('\t'.join(header) + '\n')
            return 0

        num_workers = min(MAX_WORKERS, total_count)
        # For very small inputs, skip parallelism overhead
        if total_count < 100:
            logger.info(f"Small input ({total_count} variants), using sequential processing")
            with open(output_path, 'w') as outfile:
                outfile.write('\t'.join(header) + '\n')
                for idx, variant in all_variants:
                    outfile.write(self._annotate_and_format(variant) + '\n')
            return total_count

        # Split into equal-size chunks
        chunk_size = (total_count + num_workers - 1) // num_workers
        chunks = [all_variants[i:i + chunk_size] for i in range(0, total_count, chunk_size)]
        num_workers = len(chunks)  # may be fewer if total < MAX_WORKERS * chunk_size

        logger.info(f"Parallel annotation: {total_count} variants, {num_workers} workers, ~{chunk_size} variants/worker")

        # Set the module-level engine so forked workers inherit it via COW
        global _worker_engine
        _worker_engine = self

        # Use fork context to inherit mRNA data via copy-on-write
        ctx = multiprocessing.get_context('fork')
        with ctx.Pool(processes=num_workers, initializer=_init_worker) as pool:
            chunk_results = pool.map(_annotate_chunk, chunks)

        # Flatten and sort by original index
        all_results = []
        for chunk_result in chunk_results:
            all_results.extend(chunk_result)
        all_results.sort(key=lambda x: x[0])

        # Write output
        with open(output_path, 'w') as outfile:
            outfile.write('\t'.join(header) + '\n')
            for idx, line in all_results:
                outfile.write(line + '\n')

        return total_count

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
