"""
STS file loading functionality.

The parsing logic lives in MerPCR.load_sts_file so that there is a single
implementation of the me-PCR rules. This module keeps a small loader class
for callers that only need the parsed records and hash table.
"""

import logging
from typing import Dict, List

from ..core.models import STSRecord

logger = logging.getLogger(__name__)


class STSLoader:
    """Class for loading STS files."""

    def __init__(
        self,
        wordsize: int,
        margin: int,
        default_pcr_size: int,
        iupac_mode: int = 0,
        max_sts_line_length: int = 1022,
        one_sided: int = 0,
    ):
        """Initialize STS loader with parameters."""
        self.wordsize = wordsize
        self.margin = margin
        self.default_pcr_size = default_pcr_size
        self.iupac_mode = iupac_mode
        self.max_sts_line_length = max_sts_line_length
        self.one_sided = one_sided

    def load_file(self, filename: str) -> tuple[List[STSRecord], Dict[int, List[STSRecord]], int]:
        """
        Load STS records from a tab-delimited file.

        Args:
            filename: Path to the STS file

        Returns:
            Tuple of (sts_records, sts_table, max_pcr_size). On failure the
            list and table are empty and max_pcr_size is 0.
        """
        from ..core.engine import MerPCR

        engine = MerPCR(
            wordsize=self.wordsize,
            margin=self.margin,
            iupac_mode=self.iupac_mode,
            default_pcr_size=self.default_pcr_size,
            max_sts_line_length=self.max_sts_line_length,
            one_sided=self.one_sided,
        )
        if not engine.load_sts_file(filename):
            return [], {}, 0
        return engine.sts_records, engine.sts_table, engine.max_pcr_size

    def _parse_pcr_size(self, pcr_size_str: str) -> int:
        """Parse PCR size from string, handling ranges (see MerPCR._parse_pcr_size)."""
        from ..core.engine import MerPCR

        engine = MerPCR(
            wordsize=self.wordsize, margin=self.margin, default_pcr_size=self.default_pcr_size
        )
        return engine._parse_pcr_size(pcr_size_str)[0]
