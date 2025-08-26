"""
Utility functions for merPCR.
"""

from typing import Tuple, Dict

# Global constants
AMBIG = 100  # Ambiguous base code

# Initialize lookup tables
_scode = [AMBIG] * 256
_scode[ord('A')] = _scode[ord('a')] = 0
_scode[ord('C')] = _scode[ord('c')] = 1
_scode[ord('G')] = _scode[ord('g')] = 2
_scode[ord('T')] = _scode[ord('t')] = 3
_scode[ord('U')] = _scode[ord('u')] = 3

# Complement table
_compl = {
    'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A', 'U': 'A',
    'B': 'V', 'D': 'H', 'H': 'D', 'K': 'M',
    'M': 'K', 'N': 'N', 'R': 'Y', 'S': 'S',
    'V': 'B', 'W': 'W', 'X': 'X', 'Y': 'R'
}
# Add lowercase complements
for k, v in list(_compl.items()):
    _compl[k.lower()] = v.lower()


def reverse_complement(sequence: str) -> str:
    """Return the reverse complement of a DNA sequence."""
    return ''.join(_compl.get(base, 'N') for base in reversed(sequence))


def hash_value(primer: str, wordsize: int) -> Tuple[int, int]:
    """
    Compute a hash value for the specified primer.

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

    # Try all possible positions for the hash word
    for offset in range(primer_len - wordsize + 1):
        hash_val = 0
        valid_hash = True

        for i in range(wordsize):
            base = primer[offset + i]
            code = _scode[ord(base)]
            if code == AMBIG:
                valid_hash = False
                break

            hash_val = (hash_val << 2) | code

        if valid_hash:
            return offset, hash_val

    return -1, 0


def init_iupac_tables(iupac_mode: bool = False) -> Dict:
    """Initialize IUPAC lookup tables if needed."""
    if not iupac_mode:
        return {}
    
    iupac_mapping = {
        'A': "A", 'C': "C", 'G': "G", 'T': "TU", 'U': "TU",
        'R': "AGR", 'Y': "CTUY", 'M': "ACM", 'K': "GTUK", 'S': "CGS", 'W': "ATUW",
        'B': "CGTUYKSB", 'D': "AGTURKWD", 'H': "ACTUYMWH", 'V': "ACGRMSV",
        'N': "ACGTURYMKSWBDHVN"
    }
    
    # Add lowercase entries
    for k, v in list(iupac_mapping.items()):
        iupac_mapping[k.lower()] = v
    
    return iupac_mapping