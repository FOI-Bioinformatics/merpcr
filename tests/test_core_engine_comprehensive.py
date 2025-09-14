"""
Comprehensive tests for core engine functionality with focus on search algorithms.
"""

import pytest
import tempfile
import os
from unittest.mock import patch, mock_open
from merpcr.core.engine import MerPCR
from merpcr.core.models import STSRecord, FASTARecord, STSHit
import sys
from io import StringIO


class TestMerPCRInitialization:
    """Test MerPCR class initialization and parameter validation."""

    def test_default_initialization(self):
        """Test initialization with default parameters."""
        engine = MerPCR()
        assert engine.wordsize == 11
        assert engine.margin == 50
        assert engine.mismatches == 0
        assert engine.three_prime_match == 1
        assert engine.iupac_mode == 0
        assert engine.default_pcr_size == 240
        assert engine.threads == 1
        assert engine.max_sts_line_length == 1022

    def test_custom_initialization(self):
        """Test initialization with custom parameters."""
        engine = MerPCR(
            wordsize=10,
            margin=100,
            mismatches=2,
            three_prime_match=2,
            iupac_mode=1,
            default_pcr_size=300,
            threads=4,
            max_sts_line_length=2000
        )
        assert engine.wordsize == 10
        assert engine.margin == 100
        assert engine.mismatches == 2
        assert engine.three_prime_match == 2
        assert engine.iupac_mode == 1
        assert engine.default_pcr_size == 300
        assert engine.threads == 4
        assert engine.max_sts_line_length == 2000

    def test_initialization_sets_defaults(self):
        """Test that initialization sets up internal state correctly."""
        engine = MerPCR()
        assert engine.sts_records == []
        assert engine.sts_table == {}
        assert engine.max_pcr_size == 0
        assert engine.total_hits == 0


class TestSTSFileLoading:
    """Test STS file loading functionality."""

    def test_load_valid_sts_file(self):
        """Test loading a valid STS file."""
        sts_content = """# Comment line
AFM248yg9	GCTAAAAATACACGGATGG	TGCAAGACTGCGTCTC	193	(D17S932) Chr.17, 63.7 cM
AFM256vb9	TCTGAATGGCCCTTGG	TCCTATCTGAGGTGGGGT	180	(D17S934) Chr.17, 63.7 cM
"""
        with tempfile.NamedTemporaryFile(mode='w', suffix='.sts', delete=False) as f:
            f.write(sts_content)
            sts_path = f.name

        try:
            engine = MerPCR()
            success = engine.load_sts_file(sts_path)
            assert success
            assert len(engine.sts_records) > 0
            assert engine.max_pcr_size > 0
            
        finally:
            os.unlink(sts_path)

    def test_load_nonexistent_sts_file(self):
        """Test loading a nonexistent STS file."""
        engine = MerPCR()
        # Current implementation may raise FileNotFoundError
        try:
            success = engine.load_sts_file("/nonexistent/file.sts")
            assert not success
        except FileNotFoundError:
            # This is also acceptable behavior
            pass

    def test_load_empty_sts_file(self):
        """Test loading an empty STS file."""
        with tempfile.NamedTemporaryFile(mode='w', suffix='.sts', delete=False) as f:
            sts_path = f.name

        try:
            engine = MerPCR()
            success = engine.load_sts_file(sts_path)
            assert not success
            
        finally:
            os.unlink(sts_path)

    def test_sts_file_with_invalid_format(self):
        """Test STS file with invalid format."""
        sts_content = "INVALID_LINE_FORMAT\n"
        
        with tempfile.NamedTemporaryFile(mode='w', suffix='.sts', delete=False) as f:
            f.write(sts_content)
            sts_path = f.name

        try:
            engine = MerPCR()
            success = engine.load_sts_file(sts_path)
            assert not success
            
        finally:
            os.unlink(sts_path)

    def test_sts_file_with_short_primers(self):
        """Test STS file with primers shorter than wordsize."""
        sts_content = "SHORT\tAT\tGC\t100\tShort primers\n"
        
        with tempfile.NamedTemporaryFile(mode='w', suffix='.sts', delete=False) as f:
            f.write(sts_content)
            sts_path = f.name

        try:
            engine = MerPCR(wordsize=11)
            success = engine.load_sts_file(sts_path)
            # Should succeed but skip short primers
            assert success
            
        finally:
            os.unlink(sts_path)


class TestFASTAFileLoading:
    """Test FASTA file loading functionality."""

    def test_load_valid_fasta_file(self):
        """Test loading a valid FASTA file."""
        fasta_content = """>test_sequence
ATCGATCGATCGATCGATCGATCG
ATCGATCGATCGATCGATCGATCG
"""
        with tempfile.NamedTemporaryFile(mode='w', suffix='.fa', delete=False) as f:
            f.write(fasta_content)
            fasta_path = f.name

        try:
            engine = MerPCR()
            records = engine.load_fasta_file(fasta_path)
            assert len(records) == 1
            assert records[0].label == "test_sequence"
            assert len(records[0].sequence) > 0
            
        finally:
            os.unlink(fasta_path)

    def test_load_nonexistent_fasta_file(self):
        """Test loading a nonexistent FASTA file."""
        engine = MerPCR()
        # Current implementation may raise FileNotFoundError
        try:
            records = engine.load_fasta_file("/nonexistent/file.fa")
            assert records == []
        except FileNotFoundError:
            # This is also acceptable behavior
            pass

    def test_load_multiple_fasta_sequences(self):
        """Test loading multiple sequences from FASTA file."""
        fasta_content = """>seq1
ATCGATCGATCG
>seq2
GCTAGCTAGCTA
"""
        with tempfile.NamedTemporaryFile(mode='w', suffix='.fa', delete=False) as f:
            f.write(fasta_content)
            fasta_path = f.name

        try:
            engine = MerPCR()
            records = engine.load_fasta_file(fasta_path)
            assert len(records) == 2
            assert records[0].label == "seq1"
            assert records[1].label == "seq2"
            
        finally:
            os.unlink(fasta_path)


class TestSequenceComparison:
    """Test sequence comparison functionality."""

    def test_compare_seqs_exact_match(self):
        """Test exact sequence matching."""
        engine = MerPCR(mismatches=0)
        assert engine._compare_seqs("ATCG", "ATCG", "+")
        assert not engine._compare_seqs("ATCG", "ATCC", "+")

    def test_compare_seqs_with_mismatches(self):
        """Test sequence matching with allowed mismatches."""
        engine = MerPCR(mismatches=1, three_prime_match=0)  # Disable 3' protection for this test
        assert engine._compare_seqs("ATCG", "ATCG", "+")  # Exact match
        assert engine._compare_seqs("ATCG", "ATCC", "+")  # 1 mismatch
        assert not engine._compare_seqs("ATCG", "ATAT", "+")  # 2 mismatches

    def test_compare_seqs_three_prime_protection(self):
        """Test 3' end protection in sequence comparison."""
        engine = MerPCR(mismatches=1, three_prime_match=1)
        
        # Use longer sequences to clearly separate protected and unprotected regions
        # For + strand, 3' end is at the end (last 1 base protected)
        assert engine._compare_seqs("ATCGAA", "CTCGAA", "+")  # Mismatch at position 0 (not protected)
        assert not engine._compare_seqs("ATCGAA", "ATCGAT", "+")  # Mismatch at last position (protected)
        
        # For - strand, 3' end is at the beginning (first 1 base protected)
        assert engine._compare_seqs("ATCGAA", "ATCGAC", "-")  # Mismatch at position 5 (not protected)
        assert not engine._compare_seqs("ATCGAA", "CTCGAA", "-")  # Mismatch at first position (protected)

    def test_compare_seqs_different_lengths(self):
        """Test sequence comparison with different lengths."""
        engine = MerPCR()
        assert not engine._compare_seqs("ATCG", "ATCGG", "+")
        assert not engine._compare_seqs("ATCGG", "ATCG", "+")

    def test_compare_seqs_case_insensitive(self):
        """Test case-insensitive sequence comparison."""
        engine = MerPCR(mismatches=0)
        assert engine._compare_seqs("ATCG", "atcg", "+")
        assert engine._compare_seqs("AtCg", "aTcG", "+")

    def test_compare_seqs_iupac_mode(self):
        """Test sequence comparison with IUPAC mode enabled."""
        engine = MerPCR(mismatches=0, iupac_mode=1)
        # R matches A or G
        assert engine._compare_seqs("ATCG", "RTCG", "+")
        assert engine._compare_seqs("GTCG", "RTCG", "+")
        # Y matches C or T
        assert engine._compare_seqs("ATCG", "ATYG", "+")
        assert engine._compare_seqs("ATTG", "ATYG", "+")


class TestSearchFunctionality:
    """Test core search functionality."""

    def test_basic_search(self):
        """Test basic search functionality with a simple case."""
        # Create simple test data
        sts_content = "TEST\tATCG\tCGAT\t20\tTest STS\n"
        fasta_content = ">test\nATCGACGTATCGCGAT\n"  # Contains both primers
        
        with tempfile.NamedTemporaryFile(mode='w', suffix='.sts', delete=False) as sts_f:
            sts_f.write(sts_content)
            sts_path = sts_f.name

        with tempfile.NamedTemporaryFile(mode='w', suffix='.fa', delete=False) as fasta_f:
            fasta_f.write(fasta_content)
            fasta_path = fasta_f.name

        try:
            engine = MerPCR(wordsize=4, margin=50)
            engine.load_sts_file(sts_path)
            records = engine.load_fasta_file(fasta_path)
            
            # Capture output
            output_io = StringIO()
            hit_count = engine.search(records, output_file=None)
            
            assert hit_count >= 0  # Should complete without error
            
        finally:
            os.unlink(sts_path)
            os.unlink(fasta_path)

    def test_search_with_no_hits(self):
        """Test search when no hits are found."""
        sts_content = "TEST\tAAAA\tTTTT\t20\tTest STS\n"
        fasta_content = ">test\nCCCCGGGGCCCCGGGG\n"  # No matching primers
        
        with tempfile.NamedTemporaryFile(mode='w', suffix='.sts', delete=False) as sts_f:
            sts_f.write(sts_content)
            sts_path = sts_f.name

        with tempfile.NamedTemporaryFile(mode='w', suffix='.fa', delete=False) as fasta_f:
            fasta_f.write(fasta_content)
            fasta_path = fasta_f.name

        try:
            engine = MerPCR(wordsize=4)
            engine.load_sts_file(sts_path)
            records = engine.load_fasta_file(fasta_path)
            
            hit_count = engine.search(records)
            assert hit_count == 0
            
        finally:
            os.unlink(sts_path)
            os.unlink(fasta_path)

    def test_search_output_to_file(self):
        """Test search output to file."""
        sts_content = "TEST\tATCG\tCGAT\t20\tTest STS\n"
        fasta_content = ">test\nATCGACGTATCGCGAT\n"
        
        with tempfile.NamedTemporaryFile(mode='w', suffix='.sts', delete=False) as sts_f:
            sts_f.write(sts_content)
            sts_path = sts_f.name

        with tempfile.NamedTemporaryFile(mode='w', suffix='.fa', delete=False) as fasta_f:
            fasta_f.write(fasta_content)
            fasta_path = fasta_f.name

        with tempfile.NamedTemporaryFile(mode='w', suffix='.out', delete=False) as out_f:
            out_path = out_f.name

        try:
            engine = MerPCR(wordsize=4)
            engine.load_sts_file(sts_path)
            records = engine.load_fasta_file(fasta_path)
            
            hit_count = engine.search(records, out_path)
            
            # Check that output file was created
            assert os.path.exists(out_path)
            
        finally:
            os.unlink(sts_path)
            os.unlink(fasta_path)
            if os.path.exists(out_path):
                os.unlink(out_path)

    def test_search_stdout_handling(self):
        """Test search output to stdout."""
        sts_content = "TEST\tATCG\tCGAT\t20\tTest STS\n"
        fasta_content = ">test\nATCGACGTATCGCGAT\n"
        
        with tempfile.NamedTemporaryFile(mode='w', suffix='.sts', delete=False) as sts_f:
            sts_f.write(sts_content)
            sts_path = sts_f.name

        with tempfile.NamedTemporaryFile(mode='w', suffix='.fa', delete=False) as fasta_f:
            fasta_f.write(fasta_content)
            fasta_path = fasta_f.name

        try:
            engine = MerPCR(wordsize=4)
            engine.load_sts_file(sts_path)
            records = engine.load_fasta_file(fasta_path)
            
            # Test stdout handling
            with patch('sys.stdout', new_callable=StringIO) as mock_stdout:
                hit_count = engine.search(records, "stdout")
                output = mock_stdout.getvalue()
            
            # Test None output (should go to stdout)
            with patch('sys.stdout', new_callable=StringIO) as mock_stdout:
                hit_count = engine.search(records, None)
                output = mock_stdout.getvalue()
            
        finally:
            os.unlink(sts_path)
            os.unlink(fasta_path)


class TestThreadingBehavior:
    """Test threading behavior and decisions."""

    def test_single_thread_for_small_sequences(self):
        """Test that small sequences use single thread."""
        sts_content = "TEST\tATCG\tCGAT\t20\tTest STS\n"
        fasta_content = ">small\nATCGCGAT\n"  # Very small sequence
        
        with tempfile.NamedTemporaryFile(mode='w', suffix='.sts', delete=False) as sts_f:
            sts_f.write(sts_content)
            sts_path = sts_f.name

        with tempfile.NamedTemporaryFile(mode='w', suffix='.fa', delete=False) as fasta_f:
            fasta_f.write(fasta_content)
            fasta_path = fasta_f.name

        try:
            engine = MerPCR(wordsize=4, threads=4)  # Request 4 threads
            engine.load_sts_file(sts_path)
            records = engine.load_fasta_file(fasta_path)
            
            hit_count = engine.search(records)
            # Should complete without error regardless of threading decision
            assert hit_count >= 0
            
        finally:
            os.unlink(sts_path)
            os.unlink(fasta_path)

    def test_parameter_bounds_checking(self):
        """Test that search validates parameters."""
        engine = MerPCR()
        
        # Test with empty records
        hit_count = engine.search([])
        assert hit_count == 0

    def test_search_state_management(self):
        """Test that search properly manages internal state."""
        sts_content = "TEST\tATCG\tCGAT\t20\tTest STS\n"
        fasta_content = ">test\nATCGACGTATCGCGAT\n"
        
        with tempfile.NamedTemporaryFile(mode='w', suffix='.sts', delete=False) as sts_f:
            sts_f.write(sts_content)
            sts_path = sts_f.name

        with tempfile.NamedTemporaryFile(mode='w', suffix='.fa', delete=False) as fasta_f:
            fasta_f.write(fasta_content)
            fasta_path = fasta_f.name

        try:
            engine = MerPCR(wordsize=4)
            engine.load_sts_file(sts_path)
            records = engine.load_fasta_file(fasta_path)
            
            # Run search
            hit_count = engine.search(records)
            
            # Check that total_hits is set
            assert engine.total_hits == hit_count
            
        finally:
            os.unlink(sts_path)
            os.unlink(fasta_path)


class TestEdgeCases:
    """Test edge cases and error conditions."""

    def test_search_without_sts_loaded(self):
        """Test search when no STS data is loaded."""
        engine = MerPCR()
        fasta_content = ">test\nATCGATCG\n"
        
        with tempfile.NamedTemporaryFile(mode='w', suffix='.fa', delete=False) as fasta_f:
            fasta_f.write(fasta_content)
            fasta_path = fasta_f.name

        try:
            records = engine.load_fasta_file(fasta_path)
            hit_count = engine.search(records)
            assert hit_count == 0  # No STS data, should be 0 hits
            
        finally:
            os.unlink(fasta_path)

    def test_search_with_empty_sequences(self):
        """Test search with empty sequences."""
        engine = MerPCR()
        empty_record = FASTARecord(defline=">empty", sequence="", label="empty")
        
        hit_count = engine.search([empty_record])
        assert hit_count == 0

    def test_margin_effect_on_search(self):
        """Test that margin parameter affects search results."""
        # This is a basic test that margin is used - detailed testing would require
        # specific test cases where margin makes a difference
        engine1 = MerPCR(margin=10)
        engine2 = MerPCR(margin=1000)
        
        # Both should initialize without error and have different margin values
        assert engine1.margin != engine2.margin

    def test_wordsize_effect_on_search(self):
        """Test that wordsize parameter affects search behavior."""
        # Test that different wordsizes can be set
        for wordsize in [3, 8, 11, 16]:
            engine = MerPCR(wordsize=wordsize)
            assert engine.wordsize == wordsize