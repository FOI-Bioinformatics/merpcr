#!/usr/bin/env python3
"""
merPCR - Modern Electronic Rapid PCR

Python implementation of the me-PCR (Multithreaded Electronic PCR) program originally
developed by Gregory Schuler at NCBI and enhanced by Kevin Murphy at Children's
Hospital of Philadelphia.

This tool searches large sequences for Sequence-Tagged Site (STS) markers, which
are defined as two short subsequences (primers) separated by an approximate distance.
"""

import argparse
import concurrent.futures
import logging
import os
import re
import sys
import time
from dataclasses import dataclass, field
from enum import Enum
from pathlib import Path
from typing import Dict, List, Tuple, Optional, Set, Any, Iterator

__version__ = "1.0.0"

# Global constants
AMBIG = 100  # Ambiguous base code
SEQTYPE_NT = 2  # Nucleotide sequence type

# Default parameters
DEFAULT_MARGIN = 50  # Default margin
DEFAULT_WORDSIZE = 11  # Default word size
DEFAULT_MISMATCHES = 0  # Default number of mismatches allowed
DEFAULT_THREE_PRIME_MATCH = 1  # Default number of 3' bases which must match
DEFAULT_IUPAC_MODE = 0  # Default IUPAC mode (do not honor ambiguity symbols)
DEFAULT_THREADS = 1  # Default number of threads
DEFAULT_PCR_SIZE = 240  # Default PCR size (if not specified in STS file)
MIN_FILESIZE_FOR_THREADING = 100000  # Minimum file size to enable threading

# Min/Max values for parameters
MIN_WORDSIZE = 3
MAX_WORDSIZE = 16
MIN_MISMATCHES = 0
MAX_MISMATCHES = 10
MIN_MARGIN = 0
MAX_MARGIN = 10000
MIN_THREE_PRIME_MATCH = 0
MIN_PCR_SIZE = 1
MAX_PCR_SIZE = 10000

# Setup logging
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s - %(levelname)s - %(message)s"
)
logger = logging.getLogger("merPCR")


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
    offset: int = 0  # File offset for STS in the original file
    hash_offset: int = 0  # Offset of hash word within primer
    direct: str = '+'  # Orientation: '+' for forward, '-' for reverse
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
            if '>' in self.defline:
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


class MerPCR:
    """Main merPCR class that handles all the e-PCR functionality."""

    def __init__(self,
                 wordsize: int = DEFAULT_WORDSIZE,
                 margin: int = DEFAULT_MARGIN,
                 mismatches: int = DEFAULT_MISMATCHES,
                 three_prime_match: int = DEFAULT_THREE_PRIME_MATCH,
                 iupac_mode: int = DEFAULT_IUPAC_MODE,
                 default_pcr_size: int = DEFAULT_PCR_SIZE,
                 threads: int = DEFAULT_THREADS):
        """Initialize the MerPCR class with search parameters."""
        self.wordsize = wordsize
        self.margin = margin
        self.mismatches = mismatches
        self.three_prime_match = three_prime_match
        self.iupac_mode = iupac_mode
        self.default_pcr_size = default_pcr_size
        self.threads = threads

        # Storage for loaded data
        self.sts_records = []
        self.sts_table = {}  # Hash table for STS lookup
        self.max_pcr_size = 0
        self.total_hits = 0

        # Initialize lookup tables
        self._init_lookup_tables()

        # Validate parameters
        self._validate_parameters()

    def _validate_parameters(self):
        """Validate input parameters."""
        if not (MIN_WORDSIZE <= self.wordsize <= MAX_WORDSIZE):
            raise ValueError(f"Word size must be between {MIN_WORDSIZE} and {MAX_WORDSIZE}")

        if not (MIN_MISMATCHES <= self.mismatches <= MAX_MISMATCHES):
            raise ValueError(f"Number of mismatches must be between {MIN_MISMATCHES} and {MAX_MISMATCHES}")

        if not (MIN_MARGIN <= self.margin <= MAX_MARGIN):
            raise ValueError(f"Margin must be between {MIN_MARGIN} and {MAX_MARGIN}")

        if self.three_prime_match < MIN_THREE_PRIME_MATCH:
            raise ValueError(f"Three prime match must be at least {MIN_THREE_PRIME_MATCH}")

        if not (MIN_PCR_SIZE <= self.default_pcr_size <= MAX_PCR_SIZE):
            raise ValueError(f"Default PCR size must be between {MIN_PCR_SIZE} and {MAX_PCR_SIZE}")

    def _init_lookup_tables(self):
        """Initialize lookup tables for sequence processing."""
        # Nucleotide coding table: A=0, C=1, G=2, T=3, others=AMBIG
        self.scode = [AMBIG] * 128
        self.scode[ord('A')] = 0
        self.scode[ord('C')] = 1
        self.scode[ord('G')] = 2
        self.scode[ord('T')] = 3

        # Complement table for reverse complement
        self.compl = {}
        compl_pairs = {
            'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A',
            'B': 'V', 'D': 'H', 'H': 'D', 'K': 'M',
            'M': 'K', 'N': 'N', 'R': 'Y', 'S': 'S',
            'V': 'B', 'W': 'W', 'X': 'X', 'Y': 'R'
        }
        for k, v in compl_pairs.items():
            self.compl[k] = v

        # IUPAC ambiguity table
        self.iupac_mapping = {
            'A': "A",
            'C': "C",
            'G': "G",
            'T': "TU",
            'U': "TU",
            'R': "AGR",
            'Y': "CTUY",
            'M': "ACM",
            'K': "GTUK",
            'S': "CGS",
            'W': "ATUW",
            'B': "CGTUYKSB",
            'D': "AGTURKWD",
            'H': "ACTUYMWH",
            'V': "ACGRMSV",
            'N': "ACGTURYMKSWBDHVN"
        }

        # Initialize IUPAC match matrix if IUPAC mode is enabled
        if self.iupac_mode:
            self.iupac_match_matrix = [[False for _ in range(256)] for _ in range(256)]
            for base, matches in self.iupac_mapping.items():
                for match in matches:
                    self.iupac_match_matrix[ord(match)][ord(base)] = True

        # Ambiguity detection lookup
        self.ambig = {}
        for base in "BDHKMNRSVWXY":
            self.ambig[base] = True

    def load_sts_file(self, filename: str) -> bool:
        """
        Load STS records from a tab-delimited file.

        Args:
            filename: Path to the STS file

        Returns:
            True if successful, False otherwise
        """
        start_time = time.time()
        file_size = os.path.getsize(filename)

        if file_size == 0:
            logger.error(f"STS file '{filename}' is empty")
            return False

        logger.info(f"Reading STS file: {filename}")

        self.sts_records = []
        self.sts_table = {}
        self.max_pcr_size = 0

        bad_primers_short = 0
        bad_primers_ambig = 0
        bad_pcr_size = 0

        with open(filename, 'r') as file:
            line_no = 0

            for line in file:
                line_no += 1
                line = line.strip()

                # Skip comments and blank lines
                if not line or line.startswith('#'):
                    continue

                # Parse tab-delimited fields
                fields = line.split('\t')
                if len(fields) < 4:
                    logger.error(f"Bad STS file format at line {line_no}. Expected at least 4 fields.")
                    return False

                sts_id = fields[0]
                primer1 = fields[1].upper()
                primer2 = fields[2].upper()

                # Parse PCR size, which may be a range like "100-150"
                pcr_size_str = fields[3]
                if '-' in pcr_size_str:
                    # Handle range format
                    try:
                        size_range = pcr_size_str.split('-')
                        if len(size_range) == 2 and size_range[0] and size_range[1]:
                            low = int(size_range[0])
                            high = int(size_range[1])
                            pcr_size = (low + high) // 2
                            # Adjust margin to cover the range plus the default margin
                            extra_margin = (high - low) // 2 + 1
                            effective_margin = self.margin + extra_margin
                        else:
                            pcr_size = self.default_pcr_size
                            effective_margin = self.margin
                    except ValueError:
                        pcr_size = self.default_pcr_size
                        effective_margin = self.margin
                else:
                    try:
                        pcr_size = int(pcr_size_str)
                        if pcr_size == 0:
                            pcr_size = self.default_pcr_size
                        effective_margin = self.margin
                    except ValueError:
                        pcr_size = self.default_pcr_size
                        effective_margin = self.margin

                # Get alias if available (optional field)
                alias = fields[4] if len(fields) > 4 else ""

                # Check if primer length and PCR size are valid
                if len(primer1) < self.wordsize or len(primer2) < self.wordsize:
                    bad_primers_short += 1
                    continue

                if len(primer1) + len(primer2) > pcr_size:
                    bad_pcr_size += 1
                    pcr_size = len(primer1) + len(primer2)

                # Keep track of the maximum PCR size
                if pcr_size > self.max_pcr_size:
                    self.max_pcr_size = pcr_size

                # Create STS record
                sts = STSRecord(
                    id=sts_id,
                    primer1=primer1,
                    primer2=primer2,
                    pcr_size=pcr_size,
                    alias=alias,
                    offset=file.tell(),  # Current position in file
                    direct='+'
                )

                # Forward primer hash
                hash_offset1, hash_value1 = self._hash_value(primer1)
                if hash_offset1 >= 0:
                    sts_for = STSRecord(**{**vars(sts), 'hash_offset': hash_offset1, 'direct': '+'})
                    self._insert_sts(sts_for, hash_value1)
                else:
                    bad_primers_ambig += 1

                # Reverse primer (reverse complement of primer2)
                rev_primer2 = self._reverse_complement(primer2)
                hash_offset2, hash_value2 = self._hash_value(primer2)
                if hash_offset2 >= 0:
                    sts_rev = STSRecord(**{**vars(sts), 'hash_offset': hash_offset2, 'direct': '-'})
                    # For reverse direction, primers are swapped and complemented
                    sts_rev.primer1 = primer2
                    sts_rev.primer2 = self._reverse_complement(primer1)
                    self._insert_sts(sts_rev, hash_value2)
                else:
                    bad_primers_ambig += 1

        # Report statistics
        if bad_primers_short > 0:
            logger.warning(f"{bad_primers_short} STSs have primer shorter than word size ({self.wordsize}): not included in search")

        if bad_primers_ambig > 0:
            logger.warning(f"{bad_primers_ambig} primers have ambiguities which prevent computation of a hash value: not included in search")

        if bad_pcr_size > 0:
            logger.warning(f"{bad_pcr_size} STSs have a primer length sum greater than the pcr size: expected pcr size adjusted")

        logger.info(f"Loaded {len(self.sts_records)} STS records in {time.time() - start_time:.2f} seconds")
        return True

    def _insert_sts(self, sts: STSRecord, hash_value: int):
        """Insert an STS record into the hash table."""
        if hash_value not in self.sts_table:
            self.sts_table[hash_value] = []
        self.sts_table[hash_value].append(sts)
        self.sts_records.append(sts)

    def _hash_value(self, primer: str) -> Tuple[int, int]:
        """
        Compute a hash value for the specified primer.

        Args:
            primer: The primer sequence

        Returns:
            Tuple of (offset, hash_value). If no valid hash can be computed,
            offset will be -1.
        """
        primer_len = len(primer)
        if primer_len < self.wordsize:
            return -1, 0

        # Try all possible positions for the hash word
        for offset in range(primer_len - self.wordsize + 1):
            hash_value = 0
            valid_hash = True

            for i in range(self.wordsize):
                base = primer[offset + i]
                code = self.scode[ord(base)]
                if code == AMBIG:
                    valid_hash = False
                    break

                hash_value = (hash_value << 2) | code

            if valid_hash:
                return offset, hash_value

        return -1, 0

    def _reverse_complement(self, sequence: str) -> str:
        """Return the reverse complement of a DNA sequence."""
        return ''.join(self.compl.get(base, 'N') for base in reversed(sequence))

    def load_fasta_file(self, filename: str) -> List[FASTARecord]:
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
        with open(filename, 'r') as file:
            content = file.read()

        # Parse FASTA content
        entries = re.split(r'(?=>)', content)
        for entry in entries:
            if not entry.strip():
                continue

            lines = entry.strip().split('\n')
            if not lines[0].startswith('>'):
                logger.warning(f"Invalid FASTA format, skipping entry: {lines[0][:50]}...")
                continue

            defline = lines[0]
            sequence = ''.join(lines[1:])

            # Filter out non-ACGT characters and convert to uppercase
            filtered_seq = ''.join(c.upper() for c in sequence if c.upper() in 'ACGTBDHKMNRSVWXY')

            fasta_records.append(FASTARecord(defline=defline, sequence=filtered_seq))

        logger.info(f"Loaded {len(fasta_records)} sequences in {time.time() - start_time:.2f} seconds")
        return fasta_records

    def search(self, fasta_records: List[FASTARecord], output_file: str = None) -> int:
        """
        Search for STS markers in the provided FASTA sequences.

        Args:
            fasta_records: List of FASTARecord objects to search in
            output_file: Optional file to write results to (stdout if None)

        Returns:
            Number of hits found
        """
        total_hits = 0

        output = open(output_file, 'w') if output_file else sys.stdout

        for record in fasta_records:
            seq_label = record.label
            sequence = record.sequence
            seq_len = len(sequence)

            logger.info(f"Processing sequence: {seq_label} ({seq_len} bp)")

            # Determine number of threads to use
            num_threads = self.threads
            if seq_len < MIN_FILESIZE_FOR_THREADING:
                logger.info("Sequence too small for threading, using single thread.")
                num_threads = 1

            # Calculate overlap size
            overlap = self.max_pcr_size + self.margin - 1

            # Ensure we don't use too many threads
            while num_threads > 1 and (num_threads + 1) * overlap > seq_len:
                num_threads -= 1
                logger.info(f"Reduced threads to {num_threads} due to sequence size limitations")

            # Calculate chunk size for each thread
            chunk_size = int((seq_len - (num_threads + 1) * overlap) / num_threads) + 2 * overlap

            # Prepare thread data
            thread_data = []
            offset = 0

            for i in range(num_threads):
                length = chunk_size if i < num_threads - 1 else seq_len - offset
                thread_data.append(ThreadData(
                    thread_id=i,
                    sequence=sequence[offset:offset+length],
                    offset=offset,
                    length=length
                ))
                offset += (length - overlap)

            # Execute search in parallel
            if num_threads > 1:
                with concurrent.futures.ProcessPoolExecutor(max_workers=num_threads) as executor:
                    futures = [executor.submit(self._process_thread, data) for data in thread_data]
                    for future in concurrent.futures.as_completed(futures):
                        result = future.result()
                        thread_data[result.thread_id] = result
            else:
                # Single thread version
                thread_data[0] = self._process_thread(thread_data[0])

            # Process results, filtering redundant hits
            hits = []
            for i, data in enumerate(thread_data):
                for hit in data.hits:
                    # Skip redundant hits (those that appear in the overlap region)
                    if data.offset > 0 and hit.pos2 < overlap and i > 0:
                        continue
                    hits.append(hit)

            # Sort hits by position
            hits.sort(key=lambda h: h.pos1)

            # Report hits
            for hit in hits:
                sts = hit.sts
                pos1 = hit.pos1 + 1  # Convert to 1-based
                pos2 = hit.pos2 + 1  # Convert to 1-based

                output_line = f"{seq_label}\t{pos1}..{pos2}\t{sts.id}\t({sts.direct})"
                print(output_line, file=output)
                total_hits += 1

        if output_file:
            output.close()

        logger.info(f"Total hits found: {total_hits}")
        self.total_hits = total_hits
        return total_hits

    def _process_thread(self, thread_data: ThreadData) -> ThreadData:
        """
        Process a single thread's worth of sequence data.

        Args:
            thread_data: ThreadData object containing sequence and position info

        Returns:
            Updated ThreadData with hits found
        """
        sequence = thread_data.sequence
        seq_len = len(sequence)

        if seq_len <= self.wordsize:
            return thread_data

        # Slide a window of wordsize along the sequence
        p = 0
        h = 0
        N = 0  # Number of ambiguous bases in current window

        # Initialize the hash with the first wordsize bases
        for i in range(self.wordsize):
            h = (h << 2)
            if p + i >= seq_len:
                break

            code = self.scode[ord(sequence[p + i])]
            if code == AMBIG:
                N = self.wordsize
            else:
                if N > 0:
                    N -= 1
                h |= code

        # Slide the window along the sequence
        for pos in range(seq_len - self.wordsize + 1):
            # Check if current window has a valid hash (no ambiguities)
            if N == 0 and h in self.sts_table:
                for sts in self.sts_table[h]:
                    # Verify sequence match at hash position
                    k = pos - sts.hash_offset
                    if k >= 0 and k + len(sts.primer1) <= seq_len:
                        # Try to match this STS
                        self._match_sts(sequence, seq_len, k, sts, thread_data)

            # Shift the window by one position
            if pos + self.wordsize < seq_len:
                # Remove leftmost base from hash
                h = (h << 2) & ((1 << (2 * self.wordsize)) - 1)

                # Add new rightmost base to hash
                code = self.scode[ord(sequence[pos + self.wordsize])]
                if code == AMBIG:
                    N = self.wordsize
                else:
                    if N > 0:
                        N -= 1
                    h |= code

        return thread_data

    def _match_sts(self, sequence: str, seq_len: int, k: int, sts: STSRecord, thread_data: ThreadData) -> int:
        """
        Try to match an STS at position k in the sequence.

        Args:
            sequence: The sequence to search in
            seq_len: Length of the sequence
            k: Position to start matching
            sts: STS record to match
            thread_data: ThreadData to store hits

        Returns:
            Number of hits found
        """
        primer1 = sts.primer1
        len_p1 = len(primer1)

        # Check if the whole primer matches (not just the hash word)
        if k + len_p1 <= seq_len and self._compare_seqs(sequence[k:k+len_p1], primer1, '+'):
            primer2 = sts.primer2
            len_p2 = len(primer2)
            exp_size = sts.pcr_size

            # Calculate margins for searching the second primer
            if exp_size > seq_len - k:
                # Can't fit the expected size, try to adjust
                if seq_len - k < len_p1 + len_p2:
                    return 0  # Not enough room for both primers

                exp_size = seq_len - k
                hi_margin = 0
            else:
                hi_margin = self.margin
                # Make sure hi_margin doesn't extend beyond sequence end
                if hi_margin + exp_size > seq_len - k:
                    hi_margin = seq_len - k - exp_size

            lo_margin = self.margin
            if lo_margin > exp_size - len_p1 - len_p2:
                lo_margin = exp_size - len_p1 - len_p2

            # Try the expected position first
            p2_pos = k + exp_size - len_p2
            if p2_pos >= 0 and p2_pos + len_p2 <= seq_len and self._compare_seqs(sequence[p2_pos:p2_pos+len_p2], primer2, '-'):
                thread_data.hits.append(STSHit(pos1=k, pos2=k+exp_size-1, sts=sts))
                count = 1
            else:
                count = 0

            # Try positions within the margin
            for i in range(1, self.margin + 1):
                if i <= lo_margin:
                    p2_pos = k + exp_size - len_p2 - i
                    if p2_pos >= 0 and p2_pos + len_p2 <= seq_len and self._compare_seqs(sequence[p2_pos:p2_pos+len_p2], primer2, '-'):
                        thread_data.hits.append(STSHit(pos1=k, pos2=k+exp_size-i-1, sts=sts))
                        count += 1

                if i <= hi_margin:
                    p2_pos = k + exp_size - len_p2 + i
                    if p2_pos >= 0 and p2_pos + len_p2 <= seq_len and self._compare_seqs(sequence[p2_pos:p2_pos+len_p2], primer2, '-'):
                        thread_data.hits.append(STSHit(pos1=k, pos2=k+exp_size+i-1, sts=sts))
                        count += 1

            return count

        return 0

    def _compare_seqs(self, seq1: str, seq2: str, strand: str) -> bool:
        """
        Compare two sequences allowing for mismatches.

        Args:
            seq1: First sequence
            seq2: Second sequence
            strand: '+' for forward strand, '-' for reverse strand

        Returns:
            True if sequences match within mismatch tolerance, False otherwise
        """
        if len(seq1) != len(seq2):
            return False

        mismatches = 0
        seq_len = len(seq1)

        for i in range(seq_len):
            # Check if the position is in the 3' protected region
            in_3prime_protected = ((strand == '+' and i >= seq_len - self.three_prime_match) or
                                  (strand == '-' and i < self.three_prime_match))

            # Use IUPAC comparison if enabled
            if self.iupac_mode:
                match = self.iupac_match_matrix[ord(seq1[i])][ord(seq2[i])]
            else:
                match = seq1[i] == seq2[i]

            if not match:
                # No mismatches allowed in 3' protected region
                if in_3prime_protected:
                    return False

                mismatches += 1
                if mismatches > self.mismatches:
                    return False

        return True


def main():
    """Main function to run the merPCR program."""
    parser = argparse.ArgumentParser(
        description="merPCR - Modern Electronic Rapid PCR",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )

    parser.add_argument('sts_file', type=str, help="STS file (tab-delimited)")
    parser.add_argument('fasta_file', type=str, help="FASTA sequence file")

    parser.add_argument('-M', '--margin', type=int, default=DEFAULT_MARGIN,
                        help=f"Margin (default: {DEFAULT_MARGIN})")

    parser.add_argument('-N', '--mismatches', type=int, default=DEFAULT_MISMATCHES,
                        help=f"Number of mismatches allowed (default: {DEFAULT_MISMATCHES})")

    parser.add_argument('-W', '--wordsize', type=int, default=DEFAULT_WORDSIZE,
                        help=f"Word size (default: {DEFAULT_WORDSIZE})")

    parser.add_argument('-T', '--threads', type=int, default=DEFAULT_THREADS,
                        help=f"Number of threads (default: {DEFAULT_THREADS})")

    parser.add_argument('-X', '--three-prime-match', type=int, default=DEFAULT_THREE_PRIME_MATCH,
                        help=f"Number of 3'-ward bases in which to disallow mismatches (default: {DEFAULT_THREE_PRIME_MATCH})")

    parser.add_argument('-O', '--output', type=str, default=None,
                        help="Output file name (default: stdout)")

    parser.add_argument('-Q', '--quiet', type=int, choices=[0, 1], default=1,
                        help="Quiet flag (0=verbose, 1=quiet)")

    parser.add_argument('-Z', '--default-pcr-size', type=int, default=DEFAULT_PCR_SIZE,
                        help=f"Default PCR size (default: {DEFAULT_PCR_SIZE})")

    parser.add_argument('-I', '--iupac', type=int, choices=[0, 1], default=DEFAULT_IUPAC_MODE,
                        help="IUPAC flag (0=don't honor IUPAC ambiguity symbols, 1=honor IUPAC symbols)")

    parser.add_argument('-v', '--version', action='version',
                        version=f"merPCR version {__version__}")

    parser.add_argument('--debug', action='store_true',
                        help="Enable debug logging")

    args = parser.parse_args()

    # Set up logging based on arguments
    if args.debug:
        logger.setLevel(logging.DEBUG)
    elif args.quiet == 0:
        logger.setLevel(logging.INFO)
    else:
        logger.setLevel(logging.WARNING)

    try:
        # Initialize the merPCR instance with command line arguments
        mer_pcr = MerPCR(
            wordsize=args.wordsize,
            margin=args.margin,
            mismatches=args.mismatches,
            three_prime_match=args.three_prime_match,
            iupac_mode=args.iupac,
            default_pcr_size=args.default_pcr_size,
            threads=args.threads
        )

        # Load STS file
        if not mer_pcr.load_sts_file(args.sts_file):
            logger.error(f"Failed to load STS file: {args.sts_file}")
            return 1

        # Load FASTA file
        fasta_records = mer_pcr.load_fasta_file(args.fasta_file)
        if not fasta_records:
            logger.error(f"Failed to load FASTA file: {args.fasta_file}")
            return 1

        # Run the search
        hit_count = mer_pcr.search(fasta_records, args.output)

        logger.info(f"Search complete: {hit_count} hits found")
        return 0

    except Exception as e:
        logger.error(f"Error: {str(e)}")
        if args.debug:
            import traceback
            traceback.print_exc()
        return 1


if __name__ == "__main__":
    sys.exit(main())