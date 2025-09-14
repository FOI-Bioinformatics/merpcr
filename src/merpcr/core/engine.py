"""
merPCR core search engine.

This module contains the main MerPCR class that handles STS searching functionality.
"""

import concurrent.futures
import logging
import os
import sys
import time
from typing import Dict, List, Optional

from ..io.fasta import FASTALoader
from .models import FASTARecord, STSHit, STSRecord, ThreadData

# Constants
AMBIG = 100
MIN_FILESIZE_FOR_THREADING = 100000

# Default parameters
DEFAULT_MARGIN = 50
DEFAULT_WORDSIZE = 11
DEFAULT_MISMATCHES = 0
DEFAULT_THREE_PRIME_MATCH = 1
DEFAULT_IUPAC_MODE = 0
DEFAULT_THREADS = 1
DEFAULT_PCR_SIZE = 240

# Parameter bounds
MIN_WORDSIZE = 3
MAX_WORDSIZE = 16
MIN_MISMATCHES = 0
MAX_MISMATCHES = 10
MIN_MARGIN = 0
MAX_MARGIN = 10000
MIN_THREE_PRIME_MATCH = 0
MIN_PCR_SIZE = 1
MAX_PCR_SIZE = 10000

logger = logging.getLogger(__name__)


class MerPCR:
    """Main merPCR class that handles all the e-PCR functionality."""

    def __init__(
        self,
        wordsize: int = DEFAULT_WORDSIZE,
        margin: int = DEFAULT_MARGIN,
        mismatches: int = DEFAULT_MISMATCHES,
        three_prime_match: int = DEFAULT_THREE_PRIME_MATCH,
        iupac_mode: int = DEFAULT_IUPAC_MODE,
        default_pcr_size: int = DEFAULT_PCR_SIZE,
        threads: int = DEFAULT_THREADS,
        max_sts_line_length: int = 1022,
    ):
        """Initialize the MerPCR class with search parameters."""
        self.wordsize = wordsize
        self.margin = margin
        self.mismatches = mismatches
        self.three_prime_match = three_prime_match
        self.iupac_mode = iupac_mode
        self.default_pcr_size = default_pcr_size
        self.threads = threads
        self.max_sts_line_length = max_sts_line_length

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
            raise ValueError(
                f"Number of mismatches must be between {MIN_MISMATCHES} and {MAX_MISMATCHES}"
            )

        if not (MIN_MARGIN <= self.margin <= MAX_MARGIN):
            raise ValueError(f"Margin must be between {MIN_MARGIN} and {MAX_MARGIN}")

        if self.three_prime_match < MIN_THREE_PRIME_MATCH:
            raise ValueError(f"Three prime match must be at least {MIN_THREE_PRIME_MATCH}")

        if not (MIN_PCR_SIZE <= self.default_pcr_size <= MAX_PCR_SIZE):
            raise ValueError(f"Default PCR size must be between {MIN_PCR_SIZE} and {MAX_PCR_SIZE}")

    def _init_lookup_tables(self):
        """Initialize lookup tables for sequence processing."""
        # Nucleotide coding table: A=0, C=1, G=2, T=3, others=AMBIG
        self.scode = [AMBIG] * 256

        # Add both uppercase and lowercase codes
        self.scode[ord("A")] = self.scode[ord("a")] = 0
        self.scode[ord("C")] = self.scode[ord("c")] = 1
        self.scode[ord("G")] = self.scode[ord("g")] = 2
        self.scode[ord("T")] = self.scode[ord("t")] = 3
        self.scode[ord("U")] = self.scode[ord("u")] = 3  # Treat U (RNA) as T

        # Complement table for reverse complement
        self.compl = {}
        compl_pairs = {
            "A": "T",
            "C": "G",
            "G": "C",
            "T": "A",
            "U": "A",
            "B": "V",
            "D": "H",
            "H": "D",
            "K": "M",
            "M": "K",
            "N": "N",
            "R": "Y",
            "S": "S",
            "V": "B",
            "W": "W",
            "X": "X",
            "Y": "R",
        }
        # Add lowercase complements too
        for k, v in compl_pairs.items():
            self.compl[k] = v
            self.compl[k.lower()] = v.lower()

        # IUPAC ambiguity table
        self.iupac_mapping = {
            "A": "A",
            "C": "C",
            "G": "G",
            "T": "TU",
            "U": "TU",
            "R": "AGR",
            "Y": "CTUY",
            "M": "ACM",
            "K": "GTUK",
            "S": "CGS",
            "W": "ATUW",
            "B": "CGTUYKSB",
            "D": "AGTURKWD",
            "H": "ACTUYMWH",
            "V": "ACGRMSV",
            "N": "ACGTURYMKSWBDHVN",
            # Also add lowercase entries
            "a": "A",
            "c": "C",
            "g": "G",
            "t": "TU",
            "u": "TU",
            "r": "AGR",
            "y": "CTUY",
            "m": "ACM",
            "k": "GTUK",
            "s": "CGS",
            "w": "ATUW",
            "b": "CGTUYKSB",
            "d": "AGTURKWD",
            "h": "ACTUYMWH",
            "v": "ACGRMSV",
            "n": "ACGTURYMKSWBDHVN",
        }

        # Initialize IUPAC match matrix if IUPAC mode is enabled
        if self.iupac_mode:
            self.iupac_match_matrix = [[False for _ in range(256)] for _ in range(256)]
            for base1, matches1 in self.iupac_mapping.items():
                for base2, matches2 in self.iupac_mapping.items():
                    # Two bases match if their possible interpretations overlap
                    set1 = set(matches1.upper())
                    set2 = set(matches2.upper())
                    if set1.intersection(set2):
                        self.iupac_match_matrix[ord(base1)][ord(base2)] = True
                        # Also add the uppercase-lowercase combinations
                        self.iupac_match_matrix[ord(base1.upper())][ord(base2.lower())] = True
                        self.iupac_match_matrix[ord(base1.lower())][ord(base2.upper())] = True

        # Ambiguity detection lookup
        self.ambig = {}
        for base in "BDHKMNRSVWXYbdhkmnrsvwxy":
            self.ambig[base] = True

    def load_sts_file(self, filename: str) -> bool:
        """Load STS records from a tab-delimited file."""
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
                    return False

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
                if pcr_size > self.max_pcr_size:
                    self.max_pcr_size = pcr_size

                # Create STS record for forward direction
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
                hash_offset1, hash_value1 = self._hash_value(primer1)
                if hash_offset1 >= 0:
                    sts_for = STSRecord(**{**vars(sts), "hash_offset": hash_offset1, "direct": "+"})
                    self._insert_sts(sts_for, hash_value1)
                else:
                    bad_primers_ambig += 1

                # Reverse direction: search for primer2 (forward) followed by primer1_rc
                rev_primer1 = self._reverse_complement(primer1)
                hash_offset2, hash_value2 = self._hash_value(primer2)
                if hash_offset2 >= 0:
                    sts_rev = STSRecord(**{**vars(sts), "hash_offset": hash_offset2, "direct": "-"})
                    sts_rev.primer1 = primer2
                    sts_rev.primer2 = rev_primer1
                    self._insert_sts(sts_rev, hash_value2)
                else:
                    bad_primers_ambig += 1

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
            f"Loaded {len(self.sts_records)} STS records in {time.time() - start_time:.2f} seconds"
        )
        return True

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

    def _insert_sts(self, sts: STSRecord, hash_value: int):
        """Insert an STS record into the hash table."""
        if hash_value not in self.sts_table:
            self.sts_table[hash_value] = []
        self.sts_table[hash_value].append(sts)
        self.sts_records.append(sts)

    def _hash_value(self, primer: str) -> tuple[int, int]:
        """Compute a hash value for the specified primer."""
        primer = primer.upper()
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
        return "".join(self.compl.get(base, "N") for base in reversed(sequence))

    def load_fasta_file(self, filename: str) -> List[FASTARecord]:
        """Load sequences from a FASTA file."""
        return FASTALoader.load_file(filename)

    def search(self, fasta_records: List[FASTARecord], output_file: str = None) -> int:
        """Search for STS markers in the provided FASTA sequences."""
        total_hits = 0
        if output_file and output_file.lower() != "stdout":
            output = open(output_file, "w")
        else:
            output = sys.stdout

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
                thread_data.append(
                    ThreadData(
                        thread_id=i,
                        sequence=sequence[offset : offset + length],
                        offset=offset,
                        length=length,
                    )
                )
                offset += length - overlap

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

                output_line = f"{seq_label}\t{pos1}..{pos2}\t{sts.id}\t{sts.alias}\t({sts.direct})"
                print(output_line, file=output)
                total_hits += 1

        if output_file and output_file.lower() != "stdout":
            output.close()

        logger.info(f"Total hits found: {total_hits}")
        self.total_hits = total_hits
        return total_hits

    def _process_thread(self, thread_data: ThreadData) -> ThreadData:
        """Process a single thread's worth of sequence data."""
        sequence = thread_data.sequence.upper()
        seq_len = len(sequence)

        if seq_len <= self.wordsize:
            return thread_data

        # Slide a window of wordsize along the sequence
        p = 0
        h = 0
        N = 0  # Number of ambiguous bases in current window

        # Initialize the hash with the first wordsize bases
        for i in range(self.wordsize):
            h = h << 2
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

    def _match_sts(
        self, sequence: str, seq_len: int, k: int, sts: STSRecord, thread_data: ThreadData
    ) -> int:
        """Try to match an STS at position k in the sequence."""
        primer1 = sts.primer1
        len_p1 = len(primer1)

        # Check if the whole primer matches (not just the hash word)
        if k + len_p1 <= seq_len and self._compare_seqs(sequence[k : k + len_p1], primer1, "+"):
            primer2 = sts.primer2
            len_p2 = len(primer2)
            exp_size = sts.pcr_size

            # Calculate actual available sequence length from the end of primer1
            avail_length = seq_len - (k + len_p1)

            # Check if we can fit the second primer within available sequence
            if avail_length < len_p2:
                return 0  # Not enough room for second primer

            # Calculate margins for searching
            actual_size = avail_length + len_p1  # Total available size including primer1

            # For small sequences, adjust the expected size to not exceed available sequence
            if exp_size > actual_size:
                exp_size = actual_size
                hi_margin = 0
            else:
                hi_margin = min(self.margin, seq_len - k - exp_size)

            # Ensure lo_margin doesn't push second primer before the end of first primer
            lo_margin = min(self.margin, exp_size - len_p1 - len_p2)
            if lo_margin < 0:
                lo_margin = 0

            # Try the expected position first
            p2_pos = k + exp_size - len_p2

            # Ensure p2_pos is valid
            if k + len_p1 <= p2_pos and p2_pos + len_p2 <= seq_len:
                if self._compare_seqs(sequence[p2_pos : p2_pos + len_p2], primer2, "-"):
                    actual_product_size = (p2_pos + len_p2) - k
                    thread_data.hits.append(
                        STSHit(
                            pos1=k + thread_data.offset,
                            pos2=k + actual_product_size - 1 + thread_data.offset,
                            sts=sts,
                        )
                    )
                    count = 1
                else:
                    count = 0
            else:
                count = 0

            # Try positions within the margin
            for i in range(1, self.margin + 1):
                # Try lower margin
                if i <= lo_margin:
                    p2_pos = k + exp_size - len_p2 - i
                    # Ensure second primer is after first primer
                    if k + len_p1 <= p2_pos and p2_pos + len_p2 <= seq_len:
                        if self._compare_seqs(sequence[p2_pos : p2_pos + len_p2], primer2, "-"):
                            actual_product_size = (p2_pos + len_p2) - k
                            thread_data.hits.append(
                                STSHit(
                                    pos1=k + thread_data.offset,
                                    pos2=k + actual_product_size - 1 + thread_data.offset,
                                    sts=sts,
                                )
                            )
                            count += 1

                # Try higher margin
                if i <= hi_margin:
                    p2_pos = k + exp_size - len_p2 + i
                    if p2_pos + len_p2 <= seq_len:
                        if self._compare_seqs(sequence[p2_pos : p2_pos + len_p2], primer2, "-"):
                            actual_product_size = (p2_pos + len_p2) - k
                            thread_data.hits.append(
                                STSHit(
                                    pos1=k + thread_data.offset,
                                    pos2=k + actual_product_size - 1 + thread_data.offset,
                                    sts=sts,
                                )
                            )
                            count += 1

            return count

        return 0

    def _compare_seqs(self, seq1: str, seq2: str, strand: str) -> bool:
        """Compare two sequences allowing for mismatches."""
        if len(seq1) != len(seq2):
            return False

        mismatches = 0
        seq_len = len(seq1)

        for i in range(seq_len):
            # Check if the position is in the 3' protected region
            in_3prime_protected = (strand == "+" and i >= seq_len - self.three_prime_match) or (
                strand == "-" and i < self.three_prime_match
            )

            # Use IUPAC comparison if enabled
            if self.iupac_mode:
                # Get the characters at this position
                c1 = seq1[i].upper()
                c2 = seq2[i].upper()

                # Check if either character is in the ambiguity set
                if c1 in self.iupac_mapping and c2 in self.iupac_mapping:
                    # Get the expanded character sets
                    set1 = set(self.iupac_mapping[c1])
                    set2 = set(self.iupac_mapping[c2])

                    # If the sets have any common elements, it's a match
                    match = bool(set1.intersection(set2))
                else:
                    # If one of the characters isn't recognized, it's not a match
                    match = c1 == c2
            else:
                match = seq1[i].upper() == seq2[i].upper()

            if not match:
                # No mismatches allowed in 3' protected region
                if in_3prime_protected:
                    return False

                mismatches += 1
                if mismatches > self.mismatches:
                    return False

        return True
