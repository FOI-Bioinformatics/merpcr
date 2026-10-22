"""
Utility functions for merPCR.
"""

from typing import Dict, Tuple

# Global constants
AMBIG = 100  # Ambiguous base code

# Initialize lookup tables
_scode = [AMBIG] * 256
_scode[ord("A")] = _scode[ord("a")] = 0
_scode[ord("C")] = _scode[ord("c")] = 1
_scode[ord("G")] = _scode[ord("g")] = 2
_scode[ord("T")] = _scode[ord("t")] = 3
# As in me-PCR, only A, C, G and T can be hashed; U counts as ambiguous.

# Complement table
_compl = {
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
# Add lowercase complements
for k, v in list(_compl.items()):
    _compl[k.lower()] = v.lower()


def reverse_complement(sequence: str) -> str:
    """Return the reverse complement of a DNA sequence."""
    return "".join(_compl.get(base, "N") for base in reversed(sequence))


def hash_value(primer: str, wordsize: int) -> Tuple[int, int]:
    """
    Compute a hash value for the specified primer.

    This follows me-PCR (stsmatch.cpp, HashValue): the word is taken from the
    3' end of the primer and the search moves backwards when an ambiguous base
    prevents a hash from being formed.

    Args:
        primer: The primer sequence
        wordsize: Word size for hashing

    Returns:
        Tuple of (offset, hash_value). If no valid hash can be computed,
        offset will be -1.
    """
    primer = primer.upper()
    primer_len = len(primer)
    if primer_len < wordsize:
        return -1, 0

    offset = primer_len - wordsize
    while offset >= 0:
        hash_val = 0
        i = 0
        while i < wordsize:
            char_code = ord(primer[offset + i])
            code = _scode[char_code] if char_code < 256 else AMBIG
            if code == AMBIG:
                # Jump back past the ambiguous base, as me-PCR does
                offset -= wordsize - i - 1
                break
            hash_val = (hash_val << 2) | code
            i += 1

        if i == wordsize:
            return offset, hash_val
        offset -= 1

    return -1, 0


def init_iupac_tables(iupac_mode: bool = False) -> Dict:
    """Initialize IUPAC lookup tables if needed."""
    if not iupac_mode:
        return {}

    iupac_mapping = {
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
    }

    # Add lowercase entries
    for k, v in list(iupac_mapping.items()):
        iupac_mapping[k.lower()] = v

    return iupac_mapping
