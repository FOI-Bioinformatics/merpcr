"""
FASTA file loading functionality.
"""

import logging
import os
import time
from typing import List

from ..core.models import FASTARecord

logger = logging.getLogger(__name__)


class FASTALoader:
    """Class for loading FASTA files."""

    @staticmethod
    def load_file(filename: str) -> List[FASTARecord]:
        """
        Load sequences from a FASTA file.

        Args:
            filename: Path to the FASTA file

        Returns:
            List of FASTARecord objects
        """
        start_time = time.time()
        file_size = os.path.getsize(filename)

        if file_size == 0:
            logger.error(f"FASTA file '{filename}' is empty")
            return []

        logger.info(f"Reading FASTA file: {filename}")

        fasta_records = []
        current_defline = None
        current_sequence = []

        with open(filename, "r") as file:
            for line in file:
                line = line.strip()

                if not line:
                    continue

                if line.startswith(">"):
                    # If we were already working on a sequence, save it
                    if current_defline is not None:
                        seq = "".join(current_sequence)
                        fasta_records.append(FASTARecord(defline=current_defline, sequence=seq))

                    # Start a new sequence
                    current_defline = line
                    current_sequence = []
                else:
                    # Add to current sequence, keeping only valid nucleotide characters
                    filtered_line = "".join(c for c in line if c.upper() in "ACGTBDHKMNRSVWXY")
                    current_sequence.append(filtered_line)

        # Don't forget the last sequence
        if current_defline is not None:
            seq = "".join(current_sequence)
            fasta_records.append(FASTARecord(defline=current_defline, sequence=seq))

        logger.info(
            f"Loaded {len(fasta_records)} sequences in {time.time() - start_time:.2f} seconds"
        )
        return fasta_records
