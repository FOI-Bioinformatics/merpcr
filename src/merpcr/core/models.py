"""
Data models for merPCR.
"""

from dataclasses import dataclass, field
from enum import Enum
from typing import List


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
    alias: str = ""
    offset: int = 0  # Line number in the original file
    hash_offset: int = 0  # Offset of hash word within primer
    direct: str = "+"  # Orientation: '+' for forward, '-' for reverse
    ambig_primer: int = 0  # Flags indicating which primers have ambiguities


@dataclass
class FASTARecord:
    """Class representing a FASTA sequence record."""

    defline: str
    sequence: str
    label: str = ""

    def __post_init__(self):
        """Extract the label from the defline if not provided."""
        if not self.label:
            if ">" in self.defline:
                defline = self.defline.strip()[1:]  # Remove '>' character
            else:
                defline = self.defline.strip()

            # Extract label as the first word in the defline
            self.label = defline.split()[0]


@dataclass
class STSHit:
    """Class representing an STS hit in a sequence."""

    pos1: int
    pos2: int
    sts: STSRecord


@dataclass
class ThreadData:
    """Class for thread-specific search data."""

    thread_id: int
    sequence: str
    offset: int
    length: int
    hits: List[STSHit] = field(default_factory=list)
