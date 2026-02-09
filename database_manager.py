#!/usr/bin/env python3
"""
Database manager for ANNOVAR-fast using pysam/tabix for random access.
Handles all database formats: 1000g, generic filter, snp141, region, gene.
"""

import os
import re
import logging
from typing import Dict, List, Optional, Any

import pysam

from config import DATABASE_CONFIG, HUMANDB_TBI_DIR, NA_STRING

logger = logging.getLogger(__name__)


def _decode_annovar_escapes(val: str) -> str:
    """Decode ANNOVAR-style escape sequences like \\x2c -> comma."""
    if '\\x' in val:
        val = val.replace('\\x2c', ',').replace('\\x3b', ';').replace('\\x3d', '=')
    return val


def _encode_obs_for_query(start, end, ref, alt):
    """Encode the 'obs' field for ANNOVAR-style matching in filter databases.

    ANNOVAR's matching key is (chr, start, encoded_obs).
    For generic databases:
      - insertion (start==end, ref=='-'): obs = '0' + alt_seq
      - deletion (alt=='-'): obs = str(end - start + 1)
      - block substitution: obs = str(end - start + 1) + alt_seq
      - SNP: obs = alt_seq
    """
    if start == end and ref == '-':
        # Insertion
        return '0' + alt.upper()
    elif alt == '-':
        # Deletion
        return str(end - start + 1)
    elif start != end or (start == end and len(alt) > 1):
        # Block substitution
        return str(end - start + 1) + alt.upper()
    else:
        # SNP
        return alt.upper()


class TabixDatabase:
    """Interface for a single tabix-indexed ANNOVAR database"""

    def __init__(self, name: str, config: Dict[str, Any]):
        self.name = name
        self.config = config
        self.file_path = os.path.join(HUMANDB_TBI_DIR, config['file'])
        self.db_type = config['type']
        self.operation = config['operation']
        self.columns = config['columns']
        self.tabix = None
        self._contigs = None

        if os.path.exists(self.file_path) and os.path.exists(self.file_path + '.tbi'):
            try:
                self.tabix = pysam.TabixFile(self.file_path)
                self._contigs = set(self.tabix.contigs)
            except Exception as e:
                logger.warning(f"Failed to open {name}: {e}")

    def is_available(self) -> bool:
        return self.tabix is not None

    def close(self):
        if self.tabix:
            self.tabix.close()

    def _get_query_chrom(self, chrom: str) -> Optional[str]:
        """Get the chromosome name as it exists in this database's contigs."""
        if self._contigs is None:
            return None
        if chrom in self._contigs:
            return chrom
        # Try with/without chr prefix
        if chrom.startswith('chr'):
            bare = chrom[3:]
            if bare in self._contigs:
                return bare
        else:
            with_chr = 'chr' + chrom
            if with_chr in self._contigs:
                return with_chr
        return None

    def _fetch(self, chrom: str, start: int, end: int) -> List[str]:
        """Fetch records from tabix. Returns list of tab-separated line strings."""
        qchrom = self._get_query_chrom(chrom)
        if qchrom is None:
            return []
        try:
            return list(self.tabix.fetch(qchrom, max(0, start - 1), end))
        except ValueError:
            return []

    # -------------------------------------------------------
    # Filter database queries
    # -------------------------------------------------------

    def query_filter(self, chrom: str, start: int, end: int, ref: str, alt: str) -> Dict[str, str]:
        """Query a filter-type database. Returns {column_name: value}."""
        fmt = self.config.get('format', 'generic')
        if fmt == '1000g':
            return self._query_1000g(chrom, start, end, ref, alt)
        elif fmt == 'snp':
            return self._query_snp(chrom, start, end, ref, alt)
        else:
            return self._query_generic_filter(chrom, start, end, ref, alt)

    def _query_1000g(self, chrom: str, start: int, end: int, ref: str, alt: str) -> Dict[str, str]:
        """Query 1000 Genomes database.

        Format: chr, pos, ref, obs, AF, rsid
        The obs field is already in ANNOVAR's encoded format (pre-processed):
          SNP: alt base; Deletion: length; Insertion: 0+seq; Block sub: length+seq
        """
        na = {col: NA_STRING for col in self.columns}
        records = self._fetch(chrom, start, end + 1)
        if not records:
            return na

        query_obs = _encode_obs_for_query(start, end, ref, alt)

        for record in records:
            fields = record.split('\t')
            if len(fields) < 5:
                continue

            try:
                db_start = int(fields[1])
            except ValueError:
                continue
            db_obs = fields[3]
            db_af = fields[4]

            if db_start == start and db_obs == query_obs:
                return {self.columns[0]: db_af}

        return na

    def _query_generic_filter(self, chrom: str, start: int, end: int, ref: str, alt: str) -> Dict[str, str]:
        """Query generic ANNOVAR filter database.

        Format: chr, start, end, ref, obs, data_col1, data_col2, ...
        The database ref/obs are encoded the same way as query.
        """
        na = {col: NA_STRING for col in self.columns}
        records = self._fetch(chrom, start, end)
        if not records:
            return na

        query_obs = _encode_obs_for_query(start, end, ref, alt)

        for record in records:
            if record.startswith('#'):
                continue
            fields = record.split('\t')
            if len(fields) < 5:
                continue

            db_chr = fields[0]
            try:
                db_start = int(fields[1])
                db_end = int(fields[2])
            except ValueError:
                continue
            db_ref = fields[3].upper()
            db_obs_raw = fields[4].upper()

            # Encode db obs the same way as query
            db_obs = _encode_obs_for_query(db_start, db_end, db_ref, db_obs_raw)

            if db_start == start and db_obs == query_obs:
                # Match found - extract data columns
                data_start = self.config.get('data_start_col', 5)
                data_fields = fields[data_start:] if len(fields) > data_start else []
                result = {}
                for i, col in enumerate(self.columns):
                    if i < len(data_fields):
                        val = _decode_annovar_escapes(data_fields[i].strip())
                        result[col] = val if val != '' else NA_STRING
                    else:
                        result[col] = NA_STRING
                return result

        return na

    def _query_snp(self, chrom: str, start: int, end: int, ref: str, alt: str) -> Dict[str, str]:
        """Query snp141 database (UCSC snp format) using ANNOVAR's exact matching logic.

        Format: bin, chr, start(0-based), end, name(rsID), score, strand, refNCBI, refUCSC, observed, class, ...
        ANNOVAR matches using class-specific allele encoding and exact (chr, start+1, obs) matching.
        Only processes records with class: single, deletion, in-del, insertion.
        """
        na = {col: NA_STRING for col in self.columns}
        records = self._fetch(chrom, start, end)
        if not records:
            return na

        # Encode query obs
        query_obs = _encode_obs_for_query(start, end, ref, alt)

        for record in records:
            fields = record.split('\t')
            if len(fields) < 12:
                continue

            try:
                db_start = int(fields[2])  # 0-based
            except ValueError:
                continue

            rsid = fields[4]
            ucscallele = fields[8]  # refUCSC
            twoallele = fields[9]   # observed alleles
            snp_class = fields[11]  # class field

            # Only process recognized classes
            if snp_class not in ('single', 'deletion', 'in-del', 'insertion'):
                continue

            alleles = twoallele.split('/')
            db_start_1based = db_start + 1

            if snp_class == 'single':
                # Find non-reference alleles
                for i, allele in enumerate(alleles):
                    if ucscallele == allele:
                        for j, obs_allele in enumerate(alleles):
                            if j != i and db_start_1based == start and obs_allele == query_obs:
                                return {self.columns[0]: rsid}
            elif snp_class == 'insertion':
                ins_start = db_start_1based - 1  # ANNOVAR: $start--
                if len(alleles) > 1:
                    db_obs = '0' + alleles[1]
                    if ins_start == start and db_obs == query_obs:
                        return {self.columns[0]: rsid}
            elif snp_class == 'deletion':
                db_obs = str(len(ucscallele))
                if db_start_1based == start and db_obs == query_obs:
                    return {self.columns[0]: rsid}
            elif snp_class == 'in-del':
                if len(alleles) > 1:
                    db_obs = str(len(ucscallele)) + alleles[1]
                    if db_start_1based == start and db_obs == query_obs:
                        return {self.columns[0]: rsid}

        return na

    # -------------------------------------------------------
    # Region database queries
    # -------------------------------------------------------

    def query_region(self, chrom: str, start: int, end: int, ref: str, alt: str) -> Dict[str, str]:
        """Query a region-type database. Returns {column_name: value}."""
        na = {col: NA_STRING for col in self.columns}
        records = self._fetch(chrom, start, end)
        if not records:
            return na

        pos_cols = self.config.get('positionCols', (1, 2, 3))
        score_cols = self.config.get('scoreCols', ())
        output_cols = self.config.get('colsToOutput', (4,))
        is_phastcons = self.name.startswith('phastConsElements')
        is_cytoband = (self.name == 'cytoBand')

        found_names = []
        best_score = 0
        best_name = None

        for record in records:
            fields = record.split('\t')
            try:
                chr_idx, start_idx, end_idx = pos_cols
                db_chr = fields[chr_idx]
                db_start = int(fields[start_idx])
                db_end = int(fields[end_idx])
            except (ValueError, IndexError):
                continue

            # ANNOVAR position adjustment: $start == $end or $start++
            if db_start != db_end:
                db_start += 1

            # Check overlap: query [start, end] overlaps db [db_start, db_end]
            if start > db_end or end < db_start:
                continue

            # Get score
            score = 0
            normscore = 0
            if score_cols:
                try:
                    if len(score_cols) == 2:
                        score_val = fields[score_cols[0]]
                        normscore_val = fields[score_cols[1]]
                        # phastCons special: lod= prefix
                        if is_phastcons:
                            score_val = score_val.replace('lod=', '')
                        score = float(score_val) if score_val else 0
                        normscore = float(normscore_val) if normscore_val else 0
                    elif len(score_cols) == 1:
                        normscore_val = fields[score_cols[0]]
                        normscore = float(normscore_val) if normscore_val else 0
                        score = normscore
                except (ValueError, IndexError):
                    pass

            # Get name
            try:
                name = ':'.join(fields[c] for c in output_cols)
            except IndexError:
                name = ''

            if is_cytoband:
                # cytoBand: concatenate chr+band, strip 'chr' prefix from first part
                name = name.replace(':', '')
                # The name is already chr:band from colsToOutput=(0,3)
                # e.g., "chr9:p24.1" -> strip chr -> "9p24.1"
                name = name.replace('chr', '')

            found_names.append(name)

            if score_cols:
                if score > best_score:
                    best_score = score
                    best_name = name
                elif score == best_score and best_name and name != best_name:
                    pass  # Keep existing best_name

        if not found_names:
            return na

        if is_cytoband:
            # cytoBand: special formatting
            unique_names = list(dict.fromkeys(found_names))  # preserve order, deduplicate
            if len(unique_names) >= 2:
                # Remove leading digits from last name for range format
                last = re.sub(r'^\d+', '', unique_names[-1])
                result_name = unique_names[0] + '-' + last
            else:
                result_name = unique_names[0]
        elif score_cols and not is_phastcons and len(score_cols) >= 2:
            # Databases with 2 score columns (score + normscore): format as Score=X;Name=Y
            # ANNOVAR stores normscore in regiondb; with only 1 scoreCols, normscore is undef
            # so Score= only appears when there are 2 score columns
            unique_names = list(dict.fromkeys(found_names))
            result_name = ','.join(unique_names)
            parts = []
            if best_score:
                parts.append(f"Score={int(best_score)}")
            if result_name:
                parts.append(f"Name={result_name}")
            result_name = ';'.join(parts)
        elif is_phastcons:
            # phastCons: Score=normscore;Name=lod=score
            unique_names = list(dict.fromkeys(found_names))
            result_name = ','.join(unique_names)
            parts = []
            if best_score:
                parts.append(f"Score={int(normscore)}")
            if result_name:
                parts.append(f"Name={result_name}")
            result_name = ';'.join(parts)
        else:
            # Databases without scores: just Name=...
            unique_names = list(dict.fromkeys(found_names))
            result_name = ','.join(unique_names)
            if result_name:
                result_name = f"Name={result_name}"

        return {self.columns[0]: result_name if result_name else NA_STRING}


class DatabaseManager:
    """Manages all ANNOVAR databases"""

    def __init__(self):
        self.databases: Dict[str, TabixDatabase] = {}
        self._load_databases()

    def _load_databases(self):
        """Load all available databases"""
        logger.info("Loading ANNOVAR databases...")
        for db_name, config in DATABASE_CONFIG.items():
            if config['type'] == 'gene':
                continue  # Gene databases are handled separately
            db = TabixDatabase(db_name, config)
            if db.is_available():
                self.databases[db_name] = db
                logger.debug(f"Loaded database: {db_name}")
            else:
                logger.warning(f"Database not available: {db_name}")
        logger.info(f"Loaded {len(self.databases)} databases")

    def query_variant(self, db_name: str, chrom: str, start: int, end: int,
                      ref: str, alt: str) -> Dict[str, str]:
        """Query a specific database for a variant."""
        if db_name not in self.databases:
            config = DATABASE_CONFIG.get(db_name, {})
            columns = config.get('columns', [db_name])
            return {col: NA_STRING for col in columns}

        db = self.databases[db_name]
        if db.db_type == 'filter':
            return db.query_filter(chrom, start, end, ref, alt)
        elif db.db_type == 'region':
            return db.query_region(chrom, start, end, ref, alt)
        else:
            return {col: NA_STRING for col in db.columns}

    def close(self):
        for db in self.databases.values():
            db.close()
