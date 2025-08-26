#!/usr/bin/env python3
"""
Tests for IO modules.
"""

import unittest
import tempfile
import os
import sys
import pytest

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'src'))

from merpcr.io.fasta import FASTALoader
from merpcr.core.models import FASTARecord


@pytest.mark.unit
class TestFASTALoader(unittest.TestCase):
    """Tests for FASTA file loading."""
    
    def setUp(self):
        """Set up test files."""
        self.temp_files = []
    
    def tearDown(self):
        """Clean up test files."""
        for temp_file in self.temp_files:
            if os.path.exists(temp_file):
                os.unlink(temp_file)
    
    def create_temp_fasta(self, content):
        """Create a temporary FASTA file."""
        fd, temp_file = tempfile.mkstemp(suffix='.fa')
        with os.fdopen(fd, 'w') as f:
            f.write(content)
        self.temp_files.append(temp_file)
        return temp_file
    
    def test_single_sequence(self):
        """Test loading a single sequence."""
        content = ">seq1\nATCGATCG\n"
        temp_file = self.create_temp_fasta(content)
        
        records = FASTALoader.load_file(temp_file)
        
        self.assertEqual(len(records), 1)
        self.assertEqual(records[0].label, "seq1")
        self.assertEqual(records[0].sequence, "ATCGATCG")
    
    def test_multiple_sequences(self):
        """Test loading multiple sequences."""
        content = """>seq1
ATCGATCG
>seq2 Description
GGCCTTAA
>seq3
AAAACCCC
TTTTGGGG
"""
        temp_file = self.create_temp_fasta(content)
        
        records = FASTALoader.load_file(temp_file)
        
        self.assertEqual(len(records), 3)
        self.assertEqual(records[0].label, "seq1")
        self.assertEqual(records[0].sequence, "ATCGATCG")
        self.assertEqual(records[1].label, "seq2")
        self.assertEqual(records[1].sequence, "GGCCTTAA")
        self.assertEqual(records[2].label, "seq3")
        self.assertEqual(records[2].sequence, "AAAACCCCTTTTGGGG")
    
    def test_multiline_sequence(self):
        """Test sequences split across multiple lines."""
        content = """>multiline
AAAATTTT
CCCCGGGG
TTTTAAAA
"""
        temp_file = self.create_temp_fasta(content)
        
        records = FASTALoader.load_file(temp_file)
        
        self.assertEqual(len(records), 1)
        self.assertEqual(records[0].sequence, "AAAATTTTCCCCGGGGTTTTAAAA")
    
    def test_sequence_filtering(self):
        """Test that non-nucleotide characters are filtered."""
        content = """>filtered
ATCG123NNNN456ATCG
WXYZ789GCTA
"""
        temp_file = self.create_temp_fasta(content)
        
        records = FASTALoader.load_file(temp_file)
        
        self.assertEqual(len(records), 1)
        # Should keep valid nucleotides including ambiguity codes
        self.assertEqual(records[0].sequence, "ATCGNNNATCGWXYGCTA")
    
    def test_empty_file(self):
        """Test loading empty file."""
        temp_file = self.create_temp_fasta("")
        
        records = FASTALoader.load_file(temp_file)
        
        self.assertEqual(len(records), 0)
    
    def test_no_sequences(self):
        """Test file with headers but no sequences."""
        content = ">seq1\n>seq2\n"
        temp_file = self.create_temp_fasta(content)
        
        records = FASTALoader.load_file(temp_file)
        
        self.assertEqual(len(records), 1)
        self.assertEqual(records[0].label, "seq1")
        self.assertEqual(records[0].sequence, "")
    
    def test_blank_lines(self):
        """Test handling of blank lines."""
        content = """>seq1

ATCG

GCTA

>seq2


AAAA

"""
        temp_file = self.create_temp_fasta(content)
        
        records = FASTALoader.load_file(temp_file)
        
        self.assertEqual(len(records), 2)
        self.assertEqual(records[0].sequence, "ATCGGCTA")
        self.assertEqual(records[1].sequence, "AAAA")
    
    def test_nonexistent_file(self):
        """Test loading nonexistent file."""
        with self.assertLogs(level='ERROR') as cm:
            records = FASTALoader.load_file("/nonexistent/file.fa")
        
        self.assertEqual(len(records), 0)


@pytest.mark.integration
class TestSTSLoader(unittest.TestCase):
    """Tests for STS file loading functionality."""
    
    def setUp(self):
        """Set up test environment."""
        # We'll test STS loading through the main MerPCR class
        # since STSLoader is tightly integrated
        sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'src'))
        from merpcr import MerPCR
        self.mer_pcr = MerPCR()
        self.temp_files = []
    
    def tearDown(self):
        """Clean up test files."""
        for temp_file in self.temp_files:
            if os.path.exists(temp_file):
                os.unlink(temp_file)
    
    def create_temp_sts(self, content):
        """Create a temporary STS file."""
        fd, temp_file = tempfile.mkstemp(suffix='.sts')
        with os.fdopen(fd, 'w') as f:
            f.write(content)
        self.temp_files.append(temp_file)
        return temp_file
    
    def test_basic_sts_loading(self):
        """Test basic STS file loading."""
        content = "TEST001\tATCGATCGATCG\tCGATCGATCGAT\t200\tTest STS\n"
        temp_file = self.create_temp_sts(content)
        
        success = self.mer_pcr.load_sts_file(temp_file)
        
        self.assertTrue(success)
        # Should have 2 records (forward and reverse)
        self.assertEqual(len(self.mer_pcr.sts_records), 2)
    
    def test_sts_range_parsing(self):
        """Test STS with size ranges."""
        content = "TEST001\tATCGATCGATCG\tCGATCGATCGAT\t150-250\tTest STS\n"
        temp_file = self.create_temp_sts(content)
        
        success = self.mer_pcr.load_sts_file(temp_file)
        
        self.assertTrue(success)
        # Should use average of range: (150+250)/2 = 200
        for record in self.mer_pcr.sts_records:
            if record.direct == '+':
                self.assertEqual(record.pcr_size, 200)
    
    def test_invalid_sts_format(self):
        """Test handling of invalid STS format."""
        content = "INVALID\tONLY_TWO_FIELDS\n"
        temp_file = self.create_temp_sts(content)
        
        with self.assertLogs(level='ERROR') as cm:
            success = self.mer_pcr.load_sts_file(temp_file)
        
        self.assertFalse(success)
    
    def test_short_primers(self):
        """Test handling of primers too short for word size."""
        content = "TEST001\tAT\tGC\t100\tToo short\n"
        temp_file = self.create_temp_sts(content)
        
        with self.assertLogs(level='WARNING') as cm:
            success = self.mer_pcr.load_sts_file(temp_file)
        
        self.assertTrue(success)
        # Should have 0 records due to short primers
        self.assertEqual(len(self.mer_pcr.sts_records), 0)
    
    def test_comments_and_blank_lines(self):
        """Test handling of comments and blank lines."""
        content = """# This is a comment
TEST001\tATCGATCGATCG\tCGATCGATCGAT\t200\tTest STS

# Another comment
TEST002\tGGCCTTAAGGCC\tGGCCTTAAGGCC\t180\tAnother STS
"""
        temp_file = self.create_temp_sts(content)
        
        success = self.mer_pcr.load_sts_file(temp_file)
        
        self.assertTrue(success)
        # Should have 4 records (2 STSs × 2 directions)
        self.assertEqual(len(self.mer_pcr.sts_records), 4)


if __name__ == "__main__":
    unittest.main()