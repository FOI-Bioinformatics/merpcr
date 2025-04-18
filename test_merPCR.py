#!/usr/bin/env python3
"""
Unit tests for merPCR
"""

import unittest
import tempfile
import os
from merPCR import MerPCR, FASTARecord, STSRecord


class TestMerPCR(unittest.TestCase):
    """Tests for merPCR functionality."""

    def setUp(self):
        """Set up test data."""
        # Create a simple STS file
        self.sts_content = (
            "TEST001\tGCTAAAAATACACGGATGG\tTGCAAGACTGCGTCTC\t193\tTest STS\n"
            "TEST002\tCTTCTGCATGGATCAGAAAA\tACTCACCCAGATGATTGCTT\t103\tTest STS 2\n"
        )

        self.sts_fd, self.sts_file = tempfile.mkstemp(suffix='.sts')
        with os.fdopen(self.sts_fd, 'w') as f:
            f.write(self.sts_content)

        # Create a simple FASTA file with a match for the first STS
        self.fasta_content = (
            ">test_sequence\n"
            "ACGTACGTACGTACGTACGTGCTAAAAATACACGGATGGACGTACGTACGTACGTACGT"
            "ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACG"
            "TACGTACGTACGTACGTGAGACGCAGTCTTGCAACGTACGTACGTACGTACGTACGTAC\n"
        )

        self.fasta_fd, self.fasta_file = tempfile.mkstemp(suffix='.fa')
        with os.fdopen(self.fasta_fd, 'w') as f:
            f.write(self.fasta_content)

        # Initialize merPCR with standard parameters
        self.mer_pcr = MerPCR(
            wordsize=11,
            margin=50,
            mismatches=0,
            three_prime_match=1,
            iupac_mode=0,
            default_pcr_size=240,
            threads=1
        )

    def tearDown(self):
        """Clean up test files."""
        os.unlink(self.sts_file)
        os.unlink(self.fasta_file)

    def test_load_sts_file(self):
        """Test loading STS file."""
        success = self.mer_pcr.load_sts_file(self.sts_file)
        self.assertTrue(success)
        self.assertEqual(len(self.mer_pcr.sts_records), 4)  # 2 STSs * 2 directions

    def test_load_fasta_file(self):
        """Test loading FASTA file."""
        fasta_records = self.mer_pcr.load_fasta_file(self.fasta_file)
        self.assertEqual(len(fasta_records), 1)
        self.assertEqual(fasta_records[0].label, "test_sequence")

    def test_search(self):
        """Test searching for STSs in a FASTA file."""
        self.mer_pcr.load_sts_file(self.sts_file)
        fasta_records = self.mer_pcr.load_fasta_file(self.fasta_file)

        # Create a temporary output file
        fd, output_file = tempfile.mkstemp(suffix='.txt')
        os.close(fd)

        hit_count = self.mer_pcr.search(fasta_records, output_file)
        self.assertEqual(hit_count, 1)  # Should find one match

        # Check output file content
        with open(output_file, 'r') as f:
            output = f.read()

        self.assertIn("test_sequence", output)
        self.assertIn("TEST001", output)

        os.unlink(output_file)

    def test_reverse_complement(self):
        """Test reverse complement function."""
        self.assertEqual(self.mer_pcr._reverse_complement("ACGT"), "ACGT")
        self.assertEqual(self.mer_pcr._reverse_complement("AATT"), "AATT")
        self.assertEqual(self.mer_pcr._reverse_complement("GATC"), "GATC")
        self.assertEqual(self.mer_pcr._reverse_complement("GCTA"), "TAGC")
        self.assertEqual(self.mer_pcr._reverse_complement("ATGCN"), "NGCAT")

    def test_hash_value(self):
        """Test hash value computation."""
        offset, hash_val = self.mer_pcr._hash_value("GCTAAAAATACACGGATGG")
        self.assertGreaterEqual(offset, 0)

        # Test with ambiguous sequence
        offset, hash_val = self.mer_pcr._hash_value("NNNNNNNNNNN")
        self.assertEqual(offset, -1)

    def test_compare_seqs(self):
        """Test sequence comparison with mismatches."""
        # Exact match
        self.assertTrue(self.mer_pcr._compare_seqs("ACGT", "ACGT", '+'))

        # One mismatch - should fail with default parameters
        self.assertFalse(self.mer_pcr._compare_seqs("ACGT", "ACGA", '+'))

        # One mismatch - should pass with mismatches=1
        self.mer_pcr.mismatches = 1
        self.assertTrue(self.mer_pcr._compare_seqs("ACGT", "ACGA", '+'))

        # One mismatch at 3' end - should fail even with mismatches=1
        self.assertFalse(self.mer_pcr._compare_seqs("ACGT", "TCGT", '+'))

        # Reset to default
        self.mer_pcr.mismatches = 0


if __name__ == "__main__":
    unittest.main()