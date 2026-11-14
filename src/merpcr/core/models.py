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
    # Optional internal probe, as written in the STS file (see MerPCR probe_column).
    # Several probes are written as a comma-separated list.
    probe: Optional[str] = None

    @property
    def probes(self) -> List[str]:
        """The probes of the marker, in the order written; empty without a probe."""
        return self.probe.split(",") if self.probe else []


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
class ProbeSite:
    """
    The site of one probe inside an amplicon.

    start and end are 1-based inclusive positions in the sequence, strand is
    "+" when the probe as written matches the given sequence and "-" when its
    reverse complement does, mismatches is the number of mismatching
    positions and seq is the sequence window under the probe.
    """

    start: int
    end: int
    strand: str
    mismatches: int
    seq: str = ""


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

    The probe fields are only filled for markers that carry a probe (see the
    probe_column option of MerPCR). probe_start and probe_end are 1-based
    positions of the probe site in the sequence, probe_strand is "+" when the
    probe as written matches the given sequence and "-" when its reverse
    complement does, probe_mismatches is the number of mismatching positions
    and probe_seq is the sequence window under the probe. They are None (or
    empty) when the marker has no probe or no probe site was found.

    For a marker with several probes, probe_sites holds one entry per probe
    in the order written, a ProbeSite or None when that probe has no site,
    and the single-probe fields above describe the first probe.
    """

    seq_label: str
    start: int
    end: int
    sts: STSRecord
    mismatches1: int = 0
    mismatches2: int = 0
    seq1: str = ""
    seq2: str = ""
    probe_start: Optional[int] = None
    probe_end: Optional[int] = None
    probe_strand: Optional[str] = None
    probe_mismatches: Optional[int] = None
    probe_seq: str = ""
    probe_sites: List[Optional[ProbeSite]] = field(default_factory=list)

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

    def _sites(self) -> List[Optional[ProbeSite]]:
        """One entry per probe; falls back to the single-probe fields."""
        if self.probe_sites:
            return self.probe_sites
        if (
            self.probe_start is not None
            and self.probe_end is not None
            and self.probe_strand is not None
            and self.probe_mismatches is not None
        ):
            return [
                ProbeSite(
                    self.probe_start,
                    self.probe_end,
                    self.probe_strand,
                    self.probe_mismatches,
                    self.probe_seq,
                )
            ]
        return [None] * len(self.sts.probes)

    @property
    def probe_found(self) -> bool:
        """True when the marker has probes and a site was found for each of them."""
        sites = self._sites()
        return bool(sites) and all(site is not None for site in sites)

    def format(self) -> str:
        """
        Return the output line for this hit.

        For a marker without a probe this is the me-PCR output line. For a
        marker with probes, tab-separated fields are appended for each probe
        in the order written: the probe site and its mismatch count, or the
        single field "probe=none" when that probe has no site.
        """
        alias_field = "" if self.sts.alias is None else f"\t{self.sts.alias}"
        line = (
            f"{self.seq_label}\t{self.start}..{self.end}\t{self.sts.id}{alias_field}"
            f"\t({self.sts.direct})"
        )
        if not self.sts.probe:
            return line
        for site in self._sites():
            if site is None:
                line += "\tprobe=none"
            else:
                line += f"\tprobe={site.start}..{site.end}({site.strand})\tmm={site.mismatches}"
        return line


@dataclass
class ThreadData:
    """Class for thread-specific search data."""

    thread_id: int
    sequence: str
    offset: int
    length: int
    hits: List[STSHit] = field(default_factory=list)
