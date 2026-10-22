"""
Test hash computation equivalence between merPCR and me-PCR.

This test verifies that the hash value calculation algorithm in merPCR
exactly matches the me-PCR C++ implementation, which is critical for
producing identical search results.

Key algorithm from me-PCR (stsmatch.cpp lines 1204-1240):
- Searches BACKWARDS from (primer_len - wordsize) position
- Finds first valid hash without ambiguous bases
- Returns offset where valid hash was found
"""

import pytest

from merpcr.core.engine import MerPCR


class TestHashComputationEquivalence:
    """Test that hash computation matches me-PCR exactly."""

    def setup_method(self):
        """Set up test instance."""
        self.engine = MerPCR(wordsize=11)

    def test_hash_backward_search_simple(self):
        """Test that hash search works backwards from end of primer."""
        # Simple primer with no ambiguities - should use rightmost position
        primer = "ATCGATCGATCG"  # 12 bases
        wordsize = 11

        # me-PCR searches from position (12-11=1) backwards to 0
        # For a clean primer, it should find the hash at position 1 (rightmost)
        offset, hash_value = self.engine._hash_value(primer)

        # Expected: offset should be 1 (primer_len - wordsize)
        assert offset == 1, f"Expected offset=1 (rightmost valid position), got {offset}"

    def test_hash_backward_search_with_ambiguity_at_end(self):
        """Test hash search when ambiguity is at the end."""
        # Primer with N at the end - should search backwards
        primer = "ATCGATCGATCN"  # 12 bases, N at position 11
        wordsize = 11

        offset, hash_value = self.engine._hash_value(primer)

        # me-PCR would start at position 1 (12-11), find 'N' at position 11,
        # then jump back and try position 0, which should work
        assert offset == 0, f"Expected offset=0 (fallback position), got {offset}"

    def test_hash_backward_search_with_middle_ambiguity(self):
        """Test hash search with ambiguity in middle."""
        # Primer with N in middle
        primer = "ATCGATNATCGATCG"  # 15 bases, N at position 6
        wordsize = 11

        offset, hash_value = self.engine._hash_value(primer)

        # me-PCR starts at position 4 (15-11), checks positions 4-14
        # All possible 11-base windows include position 6 which has 'N':
        # - Position 0-10: "ATCGATNATCG" - has N at position 6 (invalid)
        # - Position 1-11: "TCGATNATCGA" - has N at position 6 (invalid)
        # - Position 2-12: "CGATNATCGAT" - has N at position 6 (invalid)
        # - Position 3-13: "GATNATCGATC" - has N at position 6 (invalid)
        # - Position 4-14: "ATNATCGATCG" - has N at position 6 (invalid)
        # Therefore, no valid hash exists for this primer
        assert offset == -1, f"Expected no valid hash (offset=-1), got {offset}"

    def test_hash_computation_values(self):
        """Test that actual hash values are computed correctly."""
        # Test the 2-bit encoding: A=0, C=1, G=2, T=3
        engine11 = MerPCR(wordsize=11)

        # Test simple sequence "AAAAAAAAAAA" (11 A's)
        primer = "AAAAAAAAAAA"
        offset, hash_value = engine11._hash_value(primer)
        # All A's = all 0's = hash value 0
        assert hash_value == 0, f"Expected hash=0 for all A's, got {hash_value}"

        # Test "TTTTTTTTTTT" (11 T's)
        primer = "TTTTTTTTTTT"
        offset, hash_value = engine11._hash_value(primer)
        # All T's = all 3's = 0b11111111111111111111 = (4^11 - 1) = 4194303
        expected = 4**11 - 1
        assert hash_value == expected, f"Expected hash={expected} for all T's, got {hash_value}"

        # Test "ACGT" repeated pattern (wordsize 11)
        # ACGTACGTACG = A(0) C(1) G(2) T(3) A(0) C(1) G(2) T(3) A(0) C(1) G(2)
        primer = "ACGTACGTACG"
        offset, hash_value = engine11._hash_value(primer)
        # Binary (22 bits for 11 bases): 00 01 10 11 00 01 10 11 00 01 10
        # = 0b0001101100011011000110 = 444102
        expected = 0b0001101100011011000110
        assert hash_value == expected, f"Expected hash={expected}, got {hash_value}"

    def test_hash_no_valid_position(self):
        """Test primer with no valid hash position."""
        # Primer shorter than wordsize
        primer = "ATCG"
        engine = MerPCR(wordsize=11)
        offset, hash_value = engine._hash_value(primer)
        assert offset == -1, "Expected offset=-1 for too-short primer"

        # Primer with too many N's
        primer = "NNNNNNNNNNNNN"
        offset, hash_value = engine._hash_value(primer)
        assert offset == -1, "Expected offset=-1 for all-N primer"

    def test_hash_backward_algorithm_detailed(self):
        """
        Detailed test of the backward search algorithm.

        From me-PCR stsmatch.cpp:1204-1240, the algorithm is:
        1. Start at offset = primer_len - m_wsize
        2. Try to compute hash from that position
        3. If ambiguous base found at position i within the word:
           - offset -= (m_wsize - i - 1)  [jump back]
           - offset--  [at end of loop]
        4. Continue until valid hash found or offset < 0
        """
        engine = MerPCR(wordsize=7)  # Smaller wordsize for easier testing

        # Test case: "ATCGATNGATCG" (12 bases, N at position 6)
        # wordsize=7, so we need positions 6-12 for rightmost hash
        primer = "ATCGATNGATCG"

        offset, hash_value = engine._hash_value(primer)

        # Starting position: 12 - 7 = 5
        # Try positions 5-11: A T C G A T N - has N at position 11 (within word)
        # N is at position i=6 within the word (0-indexed within the 7-base word)
        # Jump back: offset -= (7 - 6 - 1) = 0, then offset-- makes offset = 4

        # Try positions 4-10: C G A T N G A - has N at position 8 (position 4 in word)
        # Jump back: offset -= (7 - 4 - 1) = 2, then offset-- makes offset = 1

        # Try positions 1-7: T C G A T N G - has N at position 6 (position 5 in word)
        # Jump back: offset -= (7 - 5 - 1) = 1, then offset-- makes offset = -1

        # Try positions 0-6: A T C G A T N - has N at position 6 (position 6 in word)
        # Jump back: offset -= (7 - 6 - 1) = 0, then offset-- makes offset = -1

        # No valid hash found
        assert offset == -1, f"Expected no valid hash (offset=-1), got {offset}"

    def test_hash_matches_expected_me_pcr_behavior(self):
        """
        Test specific cases to match known me-PCR behavior.

        These test cases are designed to verify edge cases in the backward
        search algorithm.
        """
        engine = MerPCR(wordsize=11)

        # Case 1: Clean primer, should use rightmost 11 bases
        test_cases = [
            # (primer, expected_offset, description)
            ("ATCGATCGATCGATCG", 5, "16-base clean primer, offset should be 16-11=5"),
            ("ATCGATCGATC", 0, "11-base clean primer, offset should be 0"),
            ("ATCGATCGATCGATCGATCG", 9, "20-base clean primer, offset should be 20-11=9"),
        ]

        for primer, expected_offset, description in test_cases:
            offset, hash_value = engine._hash_value(primer)
            assert (
                offset == expected_offset
            ), f"{description}: Expected offset={expected_offset}, got {offset}"
            assert hash_value >= 0, f"{description}: Hash value should be non-negative"

    def test_hash_with_iupac_codes(self):
        """Test that IUPAC ambiguity codes are treated as ambiguous."""
        engine = MerPCR(wordsize=11, iupac_mode=1)

        # Test with various IUPAC codes (all should be treated as ambiguous for hashing)
        iupac_codes = "RYMKSWBDHVN"
        for code in iupac_codes:
            primer = f"ATCGATCGATC{code}"  # 12 bases with IUPAC code at end
            offset, hash_value = engine._hash_value(primer)
            # Should fall back to position 0 since position 1-11 has ambiguity
            assert (
                offset == 0
            ), f"IUPAC code {code} should force fallback to position 0, got {offset}"
