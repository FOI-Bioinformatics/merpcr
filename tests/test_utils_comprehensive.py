"""
Comprehensive tests for utility functions with full coverage.
"""

import pytest

from merpcr.core.utils import (AMBIG, _compl, _scode, hash_value,
                               init_iupac_tables, reverse_complement)


class TestReverseComplement:
    """Comprehensive tests for reverse complement function."""

    def test_basic_complement(self):
        """Test basic DNA complement."""
        assert reverse_complement("A") == "T"
        assert reverse_complement("T") == "A"
        assert reverse_complement("C") == "G"
        assert reverse_complement("G") == "C"

    def test_case_handling(self):
        """Test case handling in reverse complement."""
        assert reverse_complement("ATCG") == "CGAT"
        assert reverse_complement("atcg") == "cgat"
        assert reverse_complement("AtCg") == "cGaT"

    def test_rna_support(self):
        """Test RNA (U) support."""
        assert reverse_complement("U") == "A"
        assert reverse_complement("u") == "a"
        assert reverse_complement("AUCG") == "CGAT"  # U -> A, not U -> U

    def test_ambiguous_bases(self):
        """Test IUPAC ambiguous base support."""
        assert reverse_complement("R") == "Y"  # A/G -> C/T
        assert reverse_complement("Y") == "R"  # C/T -> A/G
        assert reverse_complement("M") == "K"  # A/C -> T/G
        assert reverse_complement("K") == "M"  # G/T -> C/A
        assert reverse_complement("S") == "S"  # G/C -> C/G
        assert reverse_complement("W") == "W"  # A/T -> T/A
        assert reverse_complement("N") == "N"  # any -> any

    def test_complex_ambiguous_codes(self):
        """Test complex ambiguous codes."""
        assert reverse_complement("B") == "V"  # CGT -> GCA
        assert reverse_complement("D") == "H"  # AGT -> ACT
        assert reverse_complement("H") == "D"  # ACT -> AGT
        assert reverse_complement("V") == "B"  # ACG -> CGT

    def test_unknown_bases(self):
        """Test handling of unknown bases."""
        assert reverse_complement("Z") == "N"  # Unknown -> N
        assert reverse_complement("*") == "N"
        assert reverse_complement("1") == "N"

    def test_empty_string(self):
        """Test empty string handling."""
        assert reverse_complement("") == ""

    def test_long_sequence(self):
        """Test longer sequences."""
        seq = "ATCGATCGATCG"
        expected = "CGATCGATCGAT"
        assert reverse_complement(seq) == expected

    def test_palindromic_sequences(self):
        """Test palindromic sequences."""
        # These should be their own reverse complement
        palindromes = ["ATAT", "GCGC", "CGCG"]  # These are actually palindromic
        for seq in palindromes:
            assert reverse_complement(seq) == seq

    def test_mixed_case_complex(self):
        """Test complex mixed case sequences."""
        seq = "AtCgTaGc"
        expected = "gCtAcGaT"
        assert reverse_complement(seq) == expected


class TestHashValue:
    """Comprehensive tests for hash value computation."""

    def test_basic_hash_computation(self):
        """Test basic hash computation."""
        # Test with wordsize 3
        offset, hash_val = hash_value("ATCG", 3)
        assert offset == 0  # Should find valid hash at beginning
        assert hash_val > 0

    def test_wordsize_too_large(self):
        """Test when wordsize is larger than primer length."""
        offset, hash_val = hash_value("ATG", 5)
        assert offset == -1
        assert hash_val == 0

    def test_all_positions(self):
        """Test hash computation at different positions."""
        primer = "NATCGTT"  # N is ambiguous, should skip to position 1
        offset, hash_val = hash_value(primer, 3)
        assert offset == 1  # Should skip the N at position 0

    def test_no_valid_hash(self):
        """Test when no valid hash can be computed."""
        primer = "NNNNN"  # All ambiguous
        offset, hash_val = hash_value(primer, 3)
        assert offset == -1
        assert hash_val == 0

    def test_case_insensitive(self):
        """Test case insensitive hashing."""
        offset1, hash1 = hash_value("ATCG", 3)
        offset2, hash2 = hash_value("atcg", 3)
        assert offset1 == offset2
        assert hash1 == hash2

    def test_different_wordsizes(self):
        """Test different word sizes."""
        primer = "ATCGATCG"

        # Test various word sizes
        for wordsize in [3, 4, 5, 6, 7, 8]:
            offset, hash_val = hash_value(primer, wordsize)
            if wordsize <= len(primer):
                assert offset >= 0
                assert hash_val > 0
            else:
                assert offset == -1

    def test_hash_uniqueness(self):
        """Test that different sequences produce different hashes."""
        primers = ["ATCG", "GCTA", "TTTT", "AAAA", "CGCG"]
        hashes = []

        for primer in primers:
            offset, hash_val = hash_value(primer, 4)
            if offset >= 0:
                hashes.append(hash_val)

        # All valid hashes should be unique for these sequences
        assert len(hashes) == len(set(hashes))

    def test_ambiguous_bases_skipped(self):
        """Test that ambiguous bases are properly skipped."""
        test_cases = [
            ("NATCG", 3, 1),  # Skip N at start
            ("ATNCG", 3, -1),  # N in middle, no valid hash
            ("ATCGN", 3, 0),  # N at end, hash at start
            ("RATCG", 3, 1),  # Skip R (ambiguous)
        ]

        for primer, wordsize, expected_offset in test_cases:
            offset, hash_val = hash_value(primer, wordsize)
            if expected_offset == -1:
                assert offset == -1
            else:
                assert offset >= expected_offset

    def test_boundary_conditions(self):
        """Test boundary conditions."""
        # Minimum wordsize
        offset, hash_val = hash_value("ATG", 3)
        assert offset == 0
        assert hash_val > 0

        # Empty primer
        offset, hash_val = hash_value("", 3)
        assert offset == -1

        # Single base
        offset, hash_val = hash_value("A", 3)
        assert offset == -1

    def test_hash_bit_operations(self):
        """Test that hash computation uses correct bit operations."""
        # Test known values
        # A=0, T=3, C=1, G=2
        # ATCG with wordsize 4 should be: (0<<6)|(3<<4)|(1<<2)|(2<<0) = 0+48+4+2 = 54
        offset, hash_val = hash_value("ATCG", 4)
        assert offset == 0
        expected_hash = (0 << 6) | (3 << 4) | (1 << 2) | (2 << 0)
        assert hash_val == expected_hash


class TestIUPACTables:
    """Comprehensive tests for IUPAC table initialization."""

    def test_iupac_disabled(self):
        """Test when IUPAC mode is disabled."""
        tables = init_iupac_tables(False)
        assert tables == {}

    def test_iupac_enabled(self):
        """Test when IUPAC mode is enabled."""
        tables = init_iupac_tables(True)
        assert len(tables) > 0

        # Check standard bases
        assert tables["A"] == "A"
        assert tables["C"] == "C"
        assert tables["G"] == "G"
        assert tables["T"] == "TU"
        assert tables["U"] == "TU"

    def test_iupac_ambiguous_codes(self):
        """Test IUPAC ambiguous code mappings."""
        tables = init_iupac_tables(True)

        # Test ambiguous codes
        assert "A" in tables["R"]  # R = A or G
        assert "G" in tables["R"]
        assert "C" in tables["Y"]  # Y = C or T
        assert "T" in tables["Y"]
        assert "A" in tables["M"]  # M = A or C
        assert "C" in tables["M"]

    def test_iupac_case_handling(self):
        """Test case handling in IUPAC tables."""
        tables = init_iupac_tables(True)

        # Test that both upper and lower case are present
        assert "A" in tables
        assert "a" in tables
        assert tables["A"] == tables["a"]

        assert "R" in tables
        assert "r" in tables
        assert tables["R"] == tables["r"]

    def test_iupac_complex_codes(self):
        """Test complex IUPAC codes."""
        tables = init_iupac_tables(True)

        # Test that complex codes contain expected bases
        assert "C" in tables["B"]  # B = C, G, T
        assert "G" in tables["B"]
        assert "T" in tables["B"]

        assert "N" in tables["N"]  # N = any base
        assert len(tables["N"]) > 10  # Should contain many bases


class TestScodeLookupTable:
    """Test the internal _scode lookup table."""

    def test_scode_basic_mapping(self):
        """Test basic nucleotide mappings."""
        assert _scode[ord("A")] == 0
        assert _scode[ord("a")] == 0
        assert _scode[ord("C")] == 1
        assert _scode[ord("c")] == 1
        assert _scode[ord("G")] == 2
        assert _scode[ord("g")] == 2
        assert _scode[ord("T")] == 3
        assert _scode[ord("t")] == 3
        assert _scode[ord("U")] == 3
        assert _scode[ord("u")] == 3

    def test_scode_ambiguous_mapping(self):
        """Test that ambiguous bases map to AMBIG."""
        ambiguous_bases = ["R", "Y", "M", "K", "S", "W", "B", "D", "H", "V", "N"]
        for base in ambiguous_bases:
            assert _scode[ord(base)] == AMBIG
            assert _scode[ord(base.lower())] == AMBIG

    def test_scode_unknown_chars(self):
        """Test that unknown characters map to AMBIG."""
        unknown_chars = ["Z", "X", "*", "1", "@"]
        for char in unknown_chars:
            assert _scode[ord(char)] == AMBIG


class TestComplLookupTable:
    """Test the internal _compl lookup table."""

    def test_compl_basic_mapping(self):
        """Test basic complement mappings."""
        assert _compl["A"] == "T"
        assert _compl["T"] == "A"
        assert _compl["C"] == "G"
        assert _compl["G"] == "C"
        assert _compl["U"] == "A"

    def test_compl_case_mapping(self):
        """Test case-sensitive complement mappings."""
        assert _compl["a"] == "t"
        assert _compl["t"] == "a"
        assert _compl["c"] == "g"
        assert _compl["g"] == "c"
        assert _compl["u"] == "a"

    def test_compl_ambiguous_mapping(self):
        """Test ambiguous base complement mappings."""
        expected_complements = {
            "R": "Y",
            "Y": "R",
            "M": "K",
            "K": "M",
            "S": "S",
            "W": "W",
            "B": "V",
            "V": "B",
            "D": "H",
            "H": "D",
            "N": "N",
            "X": "X",
        }

        for base, complement in expected_complements.items():
            assert _compl[base] == complement
            assert _compl[base.lower()] == complement.lower()


class TestUtilsIntegration:
    """Integration tests combining multiple utility functions."""

    def test_hash_and_reverse_complement_integration(self):
        """Test that hash works correctly with reverse complements."""
        primer = "ATCGATCG"
        rev_comp = reverse_complement(primer)

        offset1, hash1 = hash_value(primer, 6)
        offset2, hash2 = hash_value(rev_comp, 6)

        # Both should produce valid hashes
        assert offset1 >= 0
        assert offset2 >= 0
        # But hashes should be different (unless palindromic)
        if primer != rev_comp:
            assert hash1 != hash2

    def test_iupac_with_hash_computation(self):
        """Test that IUPAC tables work with sequences that would be hashed."""
        # Even though hash_value doesn't use IUPAC tables directly,
        # test that IUPAC tables contain sequences that would be hashable
        tables = init_iupac_tables(True)

        # Test that standard bases in IUPAC tables can be hashed
        test_sequences = ["ATCG", "GCTA", "TACG", "CGAT"]  # Use mixed sequences
        for test_seq in test_sequences:
            # Should be able to compute hash for sequences containing these bases
            offset, hash_val = hash_value(test_seq, 3)
            assert offset >= 0
            assert hash_val > 0
