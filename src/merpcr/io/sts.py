"""
STS file loading functionality.
"""

import logging
import os
import time
from typing import Dict, List

from ..core.models import STSRecord

logger = logging.getLogger(__name__)


class STSLoader:
    """Class for loading STS files."""

    def __init__(self, wordsize: int, margin: int, default_pcr_size: int):
        """Initialize STS loader with parameters."""
        self.wordsize = wordsize
        self.margin = margin
        self.default_pcr_size = default_pcr_size

    def load_file(self, filename: str) -> tuple[List[STSRecord], Dict[int, List[STSRecord]], int]:
        """
        Load STS records from a tab-delimited file.

        Args:
            filename: Path to the STS file

        Returns:
            Tuple of (sts_records, sts_table, max_pcr_size)
        """
        start_time = time.time()
        file_size = os.path.getsize(filename)

        if file_size == 0:
            logger.error(f"STS file '{filename}' is empty")
            return [], {}, 0

        logger.info(f"Reading STS file: {filename}")

        sts_records = []
        sts_table = {}
        max_pcr_size = 0

        bad_primers_short = 0
        bad_primers_ambig = 0
        bad_pcr_size = 0

        with open(filename, "r") as file:
            lines = file.readlines()
            line_no = 0

            for line in lines:
                line_no += 1
                line = line.strip()

                # Skip comments and blank lines
                if not line or line.startswith("#"):
                    continue

                # Parse tab-delimited fields
                fields = line.split("\t")
                if len(fields) < 4:
                    logger.error(
                        f"Bad STS file format at line {line_no}. Expected at least 4 fields."
                    )
                    return [], {}, 0

                sts_id = fields[0]
                primer1 = fields[1].upper()
                primer2 = fields[2].upper()

                # Parse PCR size
                pcr_size = self._parse_pcr_size(fields[3])
                alias = fields[4] if len(fields) > 4 else ""

                # Check if primer length and PCR size are valid
                if len(primer1) < self.wordsize or len(primer2) < self.wordsize:
                    bad_primers_short += 1
                    continue

                if len(primer1) + len(primer2) > pcr_size:
                    bad_pcr_size += 1
                    pcr_size = len(primer1) + len(primer2)

                # Keep track of the maximum PCR size
                if pcr_size > max_pcr_size:
                    max_pcr_size = pcr_size

                # Create STS records for both directions
                self._create_sts_records(
                    sts_records,
                    sts_table,
                    sts_id,
                    primer1,
                    primer2,
                    pcr_size,
                    alias,
                    line_no,
                    bad_primers_ambig,
                )

        # Report statistics
        if bad_primers_short > 0:
            logger.warning(
                f"{bad_primers_short} STSs have primer shorter than word size ({self.wordsize}): not included in search"
            )

        if bad_primers_ambig > 0:
            logger.warning(
                f"{bad_primers_ambig} primers have ambiguities which prevent computation of a hash value: not included in search"
            )

        if bad_pcr_size > 0:
            logger.warning(
                f"{bad_pcr_size} STSs have a primer length sum greater than the pcr size: expected pcr size adjusted"
            )

        logger.info(
            f"Loaded {len(sts_records)} STS records in {time.time() - start_time:.2f} seconds"
        )
        return sts_records, sts_table, max_pcr_size

    def _parse_pcr_size(self, pcr_size_str: str) -> int:
        """Parse PCR size from string, handling ranges."""
        if "-" in pcr_size_str:
            try:
                size_range = pcr_size_str.split("-")
                if len(size_range) == 2 and size_range[0] and size_range[1]:
                    low = int(size_range[0])
                    high = int(size_range[1])
                    return (low + high) // 2
                else:
                    return self.default_pcr_size
            except ValueError:
                return self.default_pcr_size
        else:
            try:
                pcr_size = int(pcr_size_str)
                return pcr_size if pcr_size > 0 else self.default_pcr_size
            except ValueError:
                return self.default_pcr_size

    def _create_sts_records(
        self,
        sts_records: List[STSRecord],
        sts_table: Dict[int, List[STSRecord]],
        sts_id: str,
        primer1: str,
        primer2: str,
        pcr_size: int,
        alias: str,
        line_no: int,
        bad_primers_ambig: int,
    ):
        """Create STS records for both forward and reverse directions."""
        from ..core.utils import hash_value, reverse_complement

        # Create base STS record
        sts = STSRecord(
            id=sts_id,
            primer1=primer1,
            primer2=primer2,
            pcr_size=pcr_size,
            alias=alias,
            offset=line_no,
            direct="+",
        )

        # Forward primer hash
        hash_offset1, hash_value1 = hash_value(primer1, self.wordsize)
        if hash_offset1 >= 0:
            sts_for = STSRecord(**{**vars(sts), "hash_offset": hash_offset1, "direct": "+"})
            self._insert_sts(sts_records, sts_table, sts_for, hash_value1)
        else:
            bad_primers_ambig += 1

        # Reverse direction: search for primer2 (forward) followed by primer1_rc
        rev_primer1 = reverse_complement(primer1)
        hash_offset2, hash_value2 = hash_value(primer2, self.wordsize)
        if hash_offset2 >= 0:
            sts_rev = STSRecord(**{**vars(sts), "hash_offset": hash_offset2, "direct": "-"})
            sts_rev.primer1 = primer2
            sts_rev.primer2 = rev_primer1
            self._insert_sts(sts_records, sts_table, sts_rev, hash_value2)
        else:
            bad_primers_ambig += 1

    def _insert_sts(
        self,
        sts_records: List[STSRecord],
        sts_table: Dict[int, List[STSRecord]],
        sts: STSRecord,
        hash_value: int,
    ):
        """Insert an STS record into the hash table."""
        if hash_value not in sts_table:
            sts_table[hash_value] = []
        sts_table[hash_value].append(sts)
        sts_records.append(sts)
