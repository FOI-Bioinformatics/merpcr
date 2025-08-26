#!/usr/bin/env python3
"""
Comprehensive tests for merPCR using actual test data from me-PCR.
"""

import unittest
import tempfile
import os
import sys
from pathlib import Path
import pytest

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'src'))

from merpcr import MerPCR, FASTARecord, STSRecord


@pytest.mark.integration
class TestMerPCRComprehensive(unittest.TestCase):
    """Comprehensive tests using real data from me-PCR."""

    def setUp(self):
        """Set up test data."""
        self.test_dir = Path(__file__).parent / "data"
        self.sts_file = self.test_dir / "test.sts"
        self.fasta_file = self.test_dir / "test.fa"
        
        # Ensure test files exist
        if not self.sts_file.exists():
            self.skipTest("test.sts file not found")
        if not self.fasta_file.exists():
            self.skipTest("test.fa file not found")

    def test_load_real_sts_file(self):
        """Test loading the real STS file from me-PCR."""
        mer_pcr = MerPCR()
        success = mer_pcr.load_sts_file(str(self.sts_file))
        
        self.assertTrue(success)
        # Should have 6 records (3 STSs * 2 directions each)
        self.assertEqual(len(mer_pcr.sts_records), 6)
        
        # Check specific STSs
        sts_ids = [sts.id for sts in mer_pcr.sts_records]
        self.assertIn("AFM256vb9", sts_ids)
        self.assertIn("AFM248yg9", sts_ids)
        self.assertIn("AFM186xa1", sts_ids)
        
        # Check that we have both directions
        directions = [sts.direct for sts in mer_pcr.sts_records]
        self.assertIn('+', directions)
        self.assertIn('-', directions)

    def test_load_real_fasta_file(self):
        """Test loading the real FASTA file from me-PCR."""
        mer_pcr = MerPCR()
        fasta_records = mer_pcr.load_fasta_file(str(self.fasta_file))
        
        self.assertEqual(len(fasta_records), 1)
        self.assertEqual(fasta_records[0].label, "L78833")
        # Check sequence length
        self.assertGreater(len(fasta_records[0].sequence), 100000)

    def test_search_real_data(self):
        """Test searching with real data - should match me-PCR output."""
        mer_pcr = MerPCR()
        
        # Load data
        success = mer_pcr.load_sts_file(str(self.sts_file))
        self.assertTrue(success)
        
        fasta_records = mer_pcr.load_fasta_file(str(self.fasta_file))
        self.assertEqual(len(fasta_records), 1)
        
        # Create temp output file
        fd, output_file = tempfile.mkstemp(suffix='.txt')
        os.close(fd)
        
        try:
            # Run search
            hit_count = mer_pcr.search(fasta_records, output_file)
            
            # Should find 1 hit (same as original me-PCR)
            self.assertEqual(hit_count, 1)
            
            # Check output content
            with open(output_file, 'r') as f:
                output = f.read().strip()
            
            # Expected: L78833	75823..76023	AFM248yg9	(-)
            self.assertIn("L78833", output)
            self.assertIn("AFM248yg9", output)
            self.assertIn("75823..76023", output)
            self.assertIn("(-)", output)
            
        finally:
            os.unlink(output_file)

    def test_parameters_match_original(self):
        """Test that default parameters match original me-PCR."""
        mer_pcr = MerPCR()
        
        # Check default values match C++ version
        self.assertEqual(mer_pcr.wordsize, 11)  # ePCR_WDSIZE_DEFAULT
        self.assertEqual(mer_pcr.margin, 50)    # ePCR_MARGIN_DEFAULT
        self.assertEqual(mer_pcr.mismatches, 0) # ePCR_MMATCH_DEFAULT
        self.assertEqual(mer_pcr.three_prime_match, 1)  # ePCR_THREE_PRIME_MATCH_DEFAULT
        self.assertEqual(mer_pcr.default_pcr_size, 240) # ePCR_DEFAULT_PCR_SIZE_DEFAULT

    def test_parameter_validation(self):
        """Test parameter validation."""
        # Test valid parameters
        mer_pcr = MerPCR(wordsize=12, margin=100, mismatches=2)
        self.assertEqual(mer_pcr.wordsize, 12)
        
        # Test invalid wordsize
        with self.assertRaises(ValueError):
            MerPCR(wordsize=2)  # Too small
        
        with self.assertRaises(ValueError):
            MerPCR(wordsize=17)  # Too large
            
        # Test invalid mismatches
        with self.assertRaises(ValueError):
            MerPCR(mismatches=11)  # Too many
            
        # Test invalid margin
        with self.assertRaises(ValueError):
            MerPCR(margin=20000)  # Too large

    def test_iupac_support(self):
        """Test IUPAC ambiguity support."""
        mer_pcr = MerPCR(iupac_mode=1)
        
        # Test IUPAC matching
        # N should match any base
        self.assertTrue(mer_pcr._compare_seqs("ACGT", "NCGT", '+'))
        self.assertTrue(mer_pcr._compare_seqs("ACGT", "ACNT", '+'))
        
        # R (A or G) should match A and G
        self.assertTrue(mer_pcr._compare_seqs("ACGT", "RCGT", '+'))
        
        # Without IUPAC mode, should not match
        mer_pcr_no_iupac = MerPCR(iupac_mode=0)
        self.assertFalse(mer_pcr_no_iupac._compare_seqs("ACGT", "NCGT", '+'))

    def test_reverse_complement_accuracy(self):
        """Test reverse complement function accuracy."""
        mer_pcr = MerPCR()
        
        test_cases = [
            ("ATGC", "GCAT"),
            ("AAAA", "TTTT"),
            ("CGCG", "CGCG"),
            ("ATCGATCG", "CGATCGAT"),
            ("AGTCAGTC", "GACTGACT"),
        ]
        
        for seq, expected in test_cases:
            result = mer_pcr._reverse_complement(seq)
            self.assertEqual(result, expected, f"RC of {seq} should be {expected}, got {result}")

    def test_hash_function(self):
        """Test hash function accuracy."""
        mer_pcr = MerPCR(wordsize=11)
        
        # Test valid sequence
        offset, hash_val = mer_pcr._hash_value("GCTAAAAATACACGGATGG")
        self.assertGreaterEqual(offset, 0)
        self.assertGreater(hash_val, 0)
        
        # Test sequence with ambiguities
        offset, hash_val = mer_pcr._hash_value("GCTNNNNNNNNGATGG")
        self.assertEqual(offset, -1)  # No valid hash possible
        
        # Test short sequence
        offset, hash_val = mer_pcr._hash_value("ACGT")  # Shorter than wordsize
        self.assertEqual(offset, -1)

    def test_threading_behavior(self):
        """Test that threading is correctly disabled for small sequences."""
        mer_pcr = MerPCR(threads=4)
        
        # Load test data
        mer_pcr.load_sts_file(str(self.sts_file))
        fasta_records = mer_pcr.load_fasta_file(str(self.fasta_file))
        
        # Should still work with threading parameter
        hit_count = mer_pcr.search(fasta_records)
        self.assertEqual(hit_count, 1)

    def test_output_format(self):
        """Test output format matches original me-PCR."""
        mer_pcr = MerPCR()
        mer_pcr.load_sts_file(str(self.sts_file))
        fasta_records = mer_pcr.load_fasta_file(str(self.fasta_file))
        
        # Create temp output file
        fd, output_file = tempfile.mkstemp(suffix='.txt')
        os.close(fd)
        
        try:
            mer_pcr.search(fasta_records, output_file)
            
            with open(output_file, 'r') as f:
                lines = f.readlines()
            
            self.assertEqual(len(lines), 1)  # Should have exactly one hit
            
            line = lines[0].strip()
            parts = line.split('\t')
            
            # Format: sequence_label    pos1..pos2    sts_id    (orientation)
            self.assertEqual(len(parts), 4)
            self.assertEqual(parts[0], "L78833")  # sequence label
            self.assertTrue(".." in parts[1])     # position range
            self.assertEqual(parts[2], "AFM248yg9")  # STS ID
            self.assertEqual(parts[3], "(-)")     # orientation
            
        finally:
            os.unlink(output_file)

    def test_margin_effect(self):
        """Test that margin parameter affects results."""
        # Test with small margin
        mer_pcr_small = MerPCR(margin=10)
        mer_pcr_small.load_sts_file(str(self.sts_file))
        fasta_records = mer_pcr_small.load_fasta_file(str(self.fasta_file))
        hits_small = mer_pcr_small.search(fasta_records)
        
        # Test with large margin
        mer_pcr_large = MerPCR(margin=100)
        mer_pcr_large.load_sts_file(str(self.sts_file))
        fasta_records = mer_pcr_large.load_fasta_file(str(self.fasta_file))
        hits_large = mer_pcr_large.search(fasta_records)
        
        # Large margin should find at least as many hits as small margin
        self.assertGreaterEqual(hits_large, hits_small)


if __name__ == "__main__":
    unittest.main()