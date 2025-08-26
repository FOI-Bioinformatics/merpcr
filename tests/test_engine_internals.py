#!/usr/bin/env python3
"""
Tests for internal engine functionality.
"""

import unittest
import sys
import os
import pytest

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'src'))

from merpcr import MerPCR
from merpcr.core.models import STSRecord, ThreadData


@pytest.mark.unit
class TestHashFunctions(unittest.TestCase):
    """Test hash computation functions."""
    
    def setUp(self):
        """Set up test engine."""
        self.mer_pcr = MerPCR(wordsize=8)  # Use smaller wordsize for testing
    
    def test_hash_value_simple(self):
        """Test hash computation for simple sequences."""
        offset, hash_val = self.mer_pcr._hash_value("AAAAAAAA")
        self.assertEqual(offset, 0)
        self.assertEqual(hash_val, 0)  # All A's = 0000...
        
        offset, hash_val = self.mer_pcr._hash_value("TTTTTTTT")
        self.assertEqual(offset, 0)
        self.assertEqual(hash_val, 65535)  # All T's = 1111... in binary for 8 bases
    
    def test_hash_value_mixed(self):
        """Test hash computation for mixed sequences."""
        offset, hash_val = self.mer_pcr._hash_value("ATCGATCG")
        self.assertEqual(offset, 0)
        # A=00, T=11, C=01, G=10 -> 00110110 00110110 = 0x3636
        
        offset, hash_val = self.mer_pcr._hash_value("GCTAGCTA")
        self.assertEqual(offset, 0)
        # G=10, C=01, T=11, A=00 -> 10011100 10011100
    
    def test_hash_value_with_ambiguities(self):
        """Test hash computation fails with ambiguous bases."""
        offset, hash_val = self.mer_pcr._hash_value("ATCGATNG")
        self.assertEqual(offset, -1)  # Should fail due to N
        
        offset, hash_val = self.mer_pcr._hash_value("ATCGATCR")
        self.assertEqual(offset, -1)  # Should fail due to R
    
    def test_hash_value_longer_sequence(self):
        """Test hash computation with longer sequences."""
        # Should find hash in the first valid position
        offset, hash_val = self.mer_pcr._hash_value("ATCGATCGATCG")
        self.assertEqual(offset, 0)  # Uses first 8 bases
        
        # Should skip ambiguous start and find hash later
        offset, hash_val = self.mer_pcr._hash_value("NNNATCGATCGATCG")
        self.assertEqual(offset, 3)  # Skips the N's
    
    def test_hash_value_too_short(self):
        """Test hash computation with sequences too short."""
        offset, hash_val = self.mer_pcr._hash_value("ATCG")  # Only 4 bases, need 8
        self.assertEqual(offset, -1)


@pytest.mark.unit
class TestSequenceComparison(unittest.TestCase):
    """Test sequence comparison functions."""
    
    def setUp(self):
        """Set up test engine."""
        self.mer_pcr = MerPCR(mismatches=1, three_prime_match=2)
    
    def test_exact_match(self):
        """Test exact sequence matches."""
        self.assertTrue(self.mer_pcr._compare_seqs("ATCG", "ATCG", '+'))
        self.assertTrue(self.mer_pcr._compare_seqs("GGCCTTAA", "GGCCTTAA", '-'))
    
    def test_mismatches_allowed(self):
        """Test mismatches within tolerance."""
        # One mismatch allowed, not in 3' region
        self.assertTrue(self.mer_pcr._compare_seqs("ATCGATCG", "TTCGATCG", '+'))  # Mismatch at pos 0
    
    def test_mismatches_3prime_protected(self):
        """Test 3' protection prevents mismatches."""
        # Mismatch in 3' protected region (last 2 bases for forward)
        self.assertFalse(self.mer_pcr._compare_seqs("ATCGATCG", "ATCGATCT", '+'))  # Mismatch at pos 7 (last)
        self.assertFalse(self.mer_pcr._compare_seqs("ATCGATCG", "ATCGATAG", '+'))  # Mismatch at pos 6 (second last)
        
        # Mismatch in 3' protected region (first 2 bases for reverse)  
        self.assertFalse(self.mer_pcr._compare_seqs("ATCGATCG", "TTCGATCG", '-'))  # Mismatch at pos 0 (first)
        self.assertFalse(self.mer_pcr._compare_seqs("ATCGATCG", "AGCGATCG", '-'))  # Mismatch at pos 1 (second)
    
    def test_too_many_mismatches(self):
        """Test rejection when too many mismatches."""
        self.mer_pcr.mismatches = 1
        # Two mismatches, only 1 allowed
        self.assertFalse(self.mer_pcr._compare_seqs("ATCGATCG", "TTCGATCT", '+'))
    
    def test_length_mismatch(self):
        """Test rejection of different length sequences."""
        self.assertFalse(self.mer_pcr._compare_seqs("ATCG", "ATCGA", '+'))
        self.assertFalse(self.mer_pcr._compare_seqs("ATCGATCG", "ATCG", '+'))
    
    def test_case_insensitive(self):
        """Test case insensitive comparison."""
        self.assertTrue(self.mer_pcr._compare_seqs("atcg", "ATCG", '+'))
        self.assertTrue(self.mer_pcr._compare_seqs("AtCg", "aTcG", '+'))


@pytest.mark.unit
class TestIUPACSupport(unittest.TestCase):
    """Test IUPAC ambiguity code support."""
    
    def setUp(self):
        """Set up test engine with IUPAC mode."""
        self.mer_pcr = MerPCR(iupac_mode=1, mismatches=0)
    
    def test_iupac_matches(self):
        """Test IUPAC ambiguity matches."""
        # N matches any base
        self.assertTrue(self.mer_pcr._compare_seqs("ATCG", "NTCG", '+'))
        self.assertTrue(self.mer_pcr._compare_seqs("ATCG", "ANCG", '+'))
        self.assertTrue(self.mer_pcr._compare_seqs("ATCG", "ATNG", '+'))
        self.assertTrue(self.mer_pcr._compare_seqs("ATCG", "ATCN", '+'))
        
        # R (A or G) matches A and G
        self.assertTrue(self.mer_pcr._compare_seqs("ATCG", "RTCG", '+'))  # R matches A
        self.assertTrue(self.mer_pcr._compare_seqs("GTCG", "RTCG", '+'))  # R matches G
        self.assertFalse(self.mer_pcr._compare_seqs("CTCG", "RTCG", '+'))  # R doesn't match C
        
        # Y (C or T) matches C and T  
        self.assertTrue(self.mer_pcr._compare_seqs("ATCG", "ATYG", '+'))  # Y matches C
        self.assertTrue(self.mer_pcr._compare_seqs("ATTG", "ATYG", '+'))  # Y matches T
    
    def test_iupac_no_match(self):
        """Test IUPAC codes that don't match."""
        # W (A or T) doesn't match C or G
        self.assertFalse(self.mer_pcr._compare_seqs("ACCG", "AWCG", '+'))
        self.assertFalse(self.mer_pcr._compare_seqs("AGCG", "AWCG", '+'))


@pytest.mark.unit
class TestReverseComplement(unittest.TestCase):
    """Test reverse complement functionality."""
    
    def setUp(self):
        """Set up test engine."""
        self.mer_pcr = MerPCR()
    
    def test_simple_reverse_complement(self):
        """Test reverse complement of simple sequences."""
        self.assertEqual(self.mer_pcr._reverse_complement("A"), "T")
        self.assertEqual(self.mer_pcr._reverse_complement("T"), "A")
        self.assertEqual(self.mer_pcr._reverse_complement("C"), "G")
        self.assertEqual(self.mer_pcr._reverse_complement("G"), "C")
    
    def test_complex_reverse_complement(self):
        """Test reverse complement of complex sequences."""
        self.assertEqual(self.mer_pcr._reverse_complement("ATCG"), "CGAT")
        self.assertEqual(self.mer_pcr._reverse_complement("GGCCTTAA"), "TTAAGGCC")
        self.assertEqual(self.mer_pcr._reverse_complement("ATCGATCGATCG"), "CGATCGATCGAT")
    
    def test_palindromic_sequences(self):
        """Test reverse complement of palindromic sequences."""
        self.assertEqual(self.mer_pcr._reverse_complement("ACGT"), "ACGT")
        self.assertEqual(self.mer_pcr._reverse_complement("GATC"), "GATC")
        self.assertEqual(self.mer_pcr._reverse_complement("AATT"), "AATT")
    
    def test_ambiguous_bases(self):
        """Test reverse complement with ambiguous bases."""
        self.assertEqual(self.mer_pcr._reverse_complement("ATCGN"), "NCGAT")
        self.assertEqual(self.mer_pcr._reverse_complement("RWYS"), "SRWY")
        self.assertEqual(self.mer_pcr._reverse_complement("BDHV"), "BDHV")  # These are self-complementary sets
    
    def test_case_preservation(self):
        """Test case preservation in reverse complement."""
        self.assertEqual(self.mer_pcr._reverse_complement("atcg"), "cgat")
        self.assertEqual(self.mer_pcr._reverse_complement("AtCg"), "cGaT")


@pytest.mark.unit
class TestParameterValidation(unittest.TestCase):
    """Test parameter validation."""
    
    def test_valid_parameters(self):
        """Test creation with valid parameters."""
        # Should not raise any exceptions
        MerPCR(wordsize=11, margin=50, mismatches=2, three_prime_match=1)
        MerPCR(wordsize=3, margin=0, mismatches=0, three_prime_match=0)
        MerPCR(wordsize=16, margin=10000, mismatches=10)
    
    def test_invalid_wordsize(self):
        """Test invalid wordsize parameters."""
        with self.assertRaises(ValueError):
            MerPCR(wordsize=2)  # Too small
        with self.assertRaises(ValueError):
            MerPCR(wordsize=17)  # Too large
    
    def test_invalid_mismatches(self):
        """Test invalid mismatch parameters."""
        with self.assertRaises(ValueError):
            MerPCR(mismatches=-1)  # Negative
        with self.assertRaises(ValueError):
            MerPCR(mismatches=11)  # Too large
    
    def test_invalid_margin(self):
        """Test invalid margin parameters."""
        with self.assertRaises(ValueError):
            MerPCR(margin=-1)  # Negative
        with self.assertRaises(ValueError):
            MerPCR(margin=20000)  # Too large
    
    def test_invalid_pcr_size(self):
        """Test invalid PCR size parameters."""
        with self.assertRaises(ValueError):
            MerPCR(default_pcr_size=0)  # Too small
        with self.assertRaises(ValueError):
            MerPCR(default_pcr_size=15000)  # Too large


if __name__ == "__main__":
    unittest.main()