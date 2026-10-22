"""
Data models for merPCR.
"""

from dataclasses import dataclass, field
from enum import Enum
from typing import List, Optional


class SeqType(Enum):
    """Sequence type enumeration."""

    AMINO_ACID = 1
    NUCLEOTIDE = 2


@dataclass
class STSRecord:
    """Class representing an STS record."""

    id: str
    primer1: str
    primer2: str
    pcr_size: int
    # Text following the fourth tab of the STS line. None means the line had no
    # alias field at all, which me-PCR reports differently from an empty alias.
    alias: Optional[str] = None
    offset: int = 0  # Line number in the original file
    hash_offset: int = 0  # Offset of hash word within primer
    direct: str = "+"  # Orientation: '+' for forward, '-' for reverse
    ambig_primer: int = 0  # Flags indicating which primers have ambiguities
    margin: int = 0  # Per-STS margin (may be adjusted for PCR size ranges)


@dataclass
class FASTARecord:
    """Class representing a FASTA sequence record."""

    defline: str
    sequence: str
    label: str = ""

    def __post_init__(self):
        """Extract the label from the defline if not provided."""
        if not self.label:
            # Follow me-PCR: skip the leading '>' and blanks, then take the text up to
            # the first space or tab. A defline with no text gives an empty label.
            defline = self.defline.lstrip(" \t").lstrip(">").lstrip(" \t")
            words = defline.split(None, 1)
            self.label = words[0] if words else ""


@dataclass
class STSHit:
    """An STS hit within one sequence (0-based inclusive positions)."""

    pos1: int
    pos2: int
    sts: STSRecord


@dataclass
class SearchHit:
    """
    A reported hit, as returned by MerPCR.find_hits.

    Positions are 1-based and inclusive, as in the text output. format()
    returns the corresponding output line without the trailing newline.

    mismatches1 and mismatches2 are the numbers of mismatching positions
    between each primer and the sequence, counted with the same rule the
    search used (exact comparison, or IUPAC expansion for primers that
    contain ambiguity codes in IUPAC mode). seq1 and seq2 are the matched
    sequence windows on the reported strand: seq1 under the first primer of
    the record and seq2 under the reverse complement of the second primer.
    """

    seq_label: str
    start: int
    end: int
    sts: STSRecord
    mismatches1: int = 0
    mismatches2: int = 0
    seq1: str = ""
    seq2: str = ""

    @property
    def mismatches(self) -> int:
        """Total number of mismatches in both primers."""
        return self.mismatches1 + self.mismatches2

    @property
    def sts_id(self) -> str:
        return self.sts.id

    @property
    def alias(self) -> Optional[str]:
        return self.sts.alias

    @property
    def strand(self) -> str:
        return self.sts.direct

    @property
    def size(self) -> int:
        """Product size in bp."""
        return self.end - self.start + 1

    def format(self) -> str:
        """Return the me-PCR output line for this hit."""
        alias_field = "" if self.sts.alias is None else f"\t{self.sts.alias}"
        return f"{self.seq_label}\t{self.start}..{self.end}\t{self.sts.id}{alias_field}\t({self.sts.direct})"


@dataclass
class ThreadData:
    """Class for thread-specific search data."""

    thread_id: int
    sequence: str
    offset: int
    length: int
    hits: List[STSHit] = field(default_factory=list)
