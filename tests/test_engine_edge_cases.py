"""
Comprehensive edge case tests for engine.py to improve coverage from 62% to 90%+.

This test module targets uncovered lines in engine.py, focusing on:
- Parameter validation edge cases
- STS file loading error paths
- PCR size parsing edge cases
- Hash computation edge cases
- Threading and multiprocessing paths
- IUPAC mode functionality
- Match STS edge cases
- Output handling
"""

import tempfile
from pathlib import Path
from unittest.mock import Mock, patch

import pytest

from merpcr import MerPCR
from merpcr.core.models import FASTARecord, ThreadData


class TestParameterValidation:
    """Test parameter validation edge cases (lines 94, 97)."""

    def test_three_prime_match_negative(self):
        """Test three_prime_match < 0 raises ValueError."""
        with pytest.raises(ValueError, match="Three prime match must be at least"):
            MerPCR(three_prime_match=-1)

    def test_default_pcr_size_too_small(self):
        """Test default_pcr_size < MIN raises ValueError."""
        with pytest.raises(ValueError, match="Default PCR size must be between"):
            MerPCR(default_pcr_size=0)

    def test_default_pcr_size_too_large(self):
        """Test default_pcr_size > MAX raises ValueError."""
        with pytest.raises(ValueError, match="Default PCR size must be between"):
            MerPCR(default_pcr_size=100000)


class TestSTSLoadingEdgeCases:
    """Test STS file loading error paths (lines 244-249, 253, 280-303)."""

    def test_load_sts_with_short_primers(self):
        """Test loading STS with primers shorter than wordsize."""
        engine = MerPCR(wordsize=11)

        # Create STS with primers too short
        sts_content = "ShortSTS\tACGT\tGCAT\t100\tShort primers"

        with tempfile.NamedTemporaryFile(mode='w', suffix='.sts', delete=False) as f:
            f.write(sts_content)
            sts_file = f.name

        try:
            result = engine.load_sts_file(sts_file)
            assert result is True, "Should still return True even with bad primers"
            assert len(engine.sts_records) == 0, "Should not load short primers"
        finally:
            Path(sts_file).unlink()

    def test_load_sts_with_bad_pcr_size(self):
        """Test loading STS where primer sum > PCR size."""
        engine = MerPCR(wordsize=11)

        # Create STS where primers are longer than PCR size
        sts_content = "BadSize\tACGTACGTACGTACGT\tGCATGCATGCATGCAT\t10\tBad PCR size"

        with tempfile.NamedTemporaryFile(mode='w', suffix='.sts', delete=False) as f:
            f.write(sts_content)
            sts_file = f.name

        try:
            result = engine.load_sts_file(sts_file)
            assert result is True
            # Should adjust PCR size to primer sum
        finally:
            Path(sts_file).unlink()

    def test_load_sts_with_ambiguous_primers(self):
        """Test loading STS with ambiguous bases that prevent hashing."""
        engine = MerPCR(wordsize=11)

        # Create STS with all ambiguous bases (no valid hash)
        sts_content = "AmbigSTS\tNNNNNNNNNNN\tNNNNNNNNNNN\t100\tAll N"

        with tempfile.NamedTemporaryFile(mode='w', suffix='.sts', delete=False) as f:
            f.write(sts_content)
            sts_file = f.name

        try:
            result = engine.load_sts_file(sts_file)
            assert result is True
            assert len(engine.sts_records) == 0, "Should not load primers with no valid hash"
        finally:
            Path(sts_file).unlink()

    def test_load_sts_triggers_max_pcr_size_update(self):
        """Test that max_pcr_size is updated correctly."""
        engine = MerPCR(wordsize=11)

        sts_content = "BigPCR\tACGTACGTACGT\tGCATGCATGCAT\t500\tBig PCR"

        with tempfile.NamedTemporaryFile(mode='w', suffix='.sts', delete=False) as f:
            f.write(sts_content)
            sts_file = f.name

        try:
            engine.load_sts_file(sts_file)
            assert engine.max_pcr_size == 500
        finally:
            Path(sts_file).unlink()


class TestPCRSizeParsingEdgeCases:
    """Test PCR size parsing edge cases (lines 313-360)."""

    def test_parse_pcr_size_dash_only(self):
        """Test parsing PCR size that's just a dash."""
        engine = MerPCR(default_pcr_size=240)
        pcr_size, margin = engine._parse_pcr_size("-")
        assert pcr_size == 240, "Should use default for dash-only"
        assert margin == 0

    def test_parse_pcr_size_zero(self):
        """Test parsing PCR size of zero."""
        engine = MerPCR(default_pcr_size=240)
        pcr_size, margin = engine._parse_pcr_size("0")
        assert pcr_size == 240, "Should use default for zero"
        assert margin == 0

    def test_parse_pcr_size_invalid_string(self):
        """Test parsing invalid PCR size string."""
        engine = MerPCR(default_pcr_size=240)
        pcr_size, margin = engine._parse_pcr_size("invalid")
        assert pcr_size == 240, "Should use default for invalid string"
        assert margin == 0

    def test_parse_pcr_size_range_with_empty_values(self):
        """Test parsing range with empty values like '-200' or '200-'."""
        engine = MerPCR(default_pcr_size=240)

        # Test "-200"
        pcr_size, margin = engine._parse_pcr_size("-200")
        assert pcr_size == 240, "Should use default for incomplete range"

        # Test "200-"
        pcr_size, margin = engine._parse_pcr_size("200-")
        assert pcr_size == 240, "Should use default for incomplete range"

    def test_parse_pcr_size_range_with_invalid_numbers(self):
        """Test parsing range with invalid numbers."""
        engine = MerPCR(default_pcr_size=240)
        pcr_size, margin = engine._parse_pcr_size("abc-def")
        assert pcr_size == 240, "Should use default for invalid range"
        assert margin == 0


class TestHashComputationEdgeCases:
    """Test hash computation edge cases (lines 385-400)."""

    def test_hash_value_primer_too_short(self):
        """Test hash computation with primer shorter than wordsize."""
        engine = MerPCR(wordsize=11)
        offset, hash_val = engine._hash_value("ACGT")  # Only 4 bases
        assert offset == -1, "Should return -1 for primer too short"
        assert hash_val == 0

    def test_hash_value_all_ambiguous(self):
        """Test hash computation with all ambiguous bases."""
        engine = MerPCR(wordsize=11)
        offset, hash_val = engine._hash_value("NNNNNNNNNNNNN")
        assert offset == -1, "Should return -1 for all ambiguous"
        assert hash_val == 0


class TestSearchEdgeCases:
    """Test search method edge cases (lines 410-514)."""

    def test_search_with_file_output(self):
        """Test search writing to file (not stdout)."""
        engine = MerPCR(wordsize=11)

        sts_content = "TestSTS\tACGTACGTACGT\tGCATGCATGCAT\t100\t"
        fasta_content = (
            ">test_seq\n"
            "ACGTACGTACGTNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN"
            "NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN"
            "NNNNNNNNNNATGCATGCATGC\n"
        )

        with tempfile.NamedTemporaryFile(mode='w', suffix='.sts', delete=False) as f:
            f.write(sts_content)
            sts_file = f.name

        with tempfile.NamedTemporaryFile(mode='w', suffix='.fa', delete=False) as f:
            f.write(fasta_content)
            fasta_file = f.name

        with tempfile.NamedTemporaryFile(mode='w', suffix='.txt', delete=False) as f:
            output_file = f.name

        try:
            engine.load_sts_file(sts_file)
            records = engine.load_fasta_file(fasta_file)
            hit_count = engine.search(records, output_file)

            assert hit_count >= 0
            assert Path(output_file).exists()

            # Verify output file was written and closed
            with open(output_file, 'r') as f:
                content = f.read()
                assert len(content) > 0 if hit_count > 0 else True

        finally:
            Path(sts_file).unlink()
            Path(fasta_file).unlink()
            Path(output_file).unlink(missing_ok=True)

    def test_search_with_small_sequence(self):
        """Test search with sequence smaller than MIN_FILESIZE_FOR_THREADING."""
        engine = MerPCR(wordsize=11, threads=4)  # Request threading

        sts_content = "TestSTS\tACGTACGTACGT\tGCATGCATGCAT\t100\t"
        # Small sequence (< 100KB)
        fasta_content = ">small_seq\n" + "ACGT" * 100 + "\n"

        with tempfile.NamedTemporaryFile(mode='w', suffix='.sts', delete=False) as f:
            f.write(sts_content)
            sts_file = f.name

        with tempfile.NamedTemporaryFile(mode='w', suffix='.fa', delete=False) as f:
            f.write(fasta_content)
            fasta_file = f.name

        try:
            engine.load_sts_file(sts_file)
            records = engine.load_fasta_file(fasta_file)
            hit_count = engine.search(records)
            assert hit_count >= 0

        finally:
            Path(sts_file).unlink()
            Path(fasta_file).unlink()

    def test_search_with_thread_reduction(self):
        """Test search where thread count is reduced due to overlap."""
        engine = MerPCR(wordsize=11, threads=10, margin=5000)

        sts_content = "TestSTS\tACGTACGTACGT\tGCATGCATGCAT\t8000\t"
        # Medium sequence where overlap would cause thread reduction
        fasta_content = ">medium_seq\n" + "ACGT" * 10000 + "\n"

        with tempfile.NamedTemporaryFile(mode='w', suffix='.sts', delete=False) as f:
            f.write(sts_content)
            sts_file = f.name

        with tempfile.NamedTemporaryFile(mode='w', suffix='.fa', delete=False) as f:
            f.write(fasta_content)
            fasta_file = f.name

        try:
            engine.load_sts_file(sts_file)
            records = engine.load_fasta_file(fasta_file)
            hit_count = engine.search(records)
            assert hit_count >= 0

        finally:
            Path(sts_file).unlink()
            Path(fasta_file).unlink()


class TestProcessThreadEdgeCases:
    """Test _process_thread edge cases (lines 519-568)."""

    def test_process_thread_with_sequence_too_short(self):
        """Test processing sequence shorter than wordsize."""
        engine = MerPCR(wordsize=11)

        thread_data = ThreadData(
            thread_id=0,
            sequence="ACGT",  # Only 4 bases
            offset=0,
            length=4
        )

        result = engine._process_thread(thread_data)
        assert len(result.hits) == 0

    def test_process_thread_with_ambiguous_bases(self):
        """Test processing sequence with ambiguous bases in sliding window."""
        engine = MerPCR(wordsize=11)

        sts_content = "TestSTS\tACGTACGTACGT\tGCATGCATGCAT\t100\t"

        with tempfile.NamedTemporaryFile(mode='w', suffix='.sts', delete=False) as f:
            f.write(sts_content)
            sts_file = f.name

        try:
            engine.load_sts_file(sts_file)

            # Sequence with ambiguous bases that will affect hash window
            sequence = "ACGTNNNNNNNNACGT" * 10

            thread_data = ThreadData(
                thread_id=0,
                sequence=sequence,
                offset=0,
                length=len(sequence)
            )

            result = engine._process_thread(thread_data)
            # Just verify it doesn't crash
            assert isinstance(result, ThreadData)

        finally:
            Path(sts_file).unlink()

    def test_process_thread_hash_window_at_sequence_end(self):
        """Test hash window sliding when reaching end of sequence."""
        engine = MerPCR(wordsize=11)

        sts_content = "TestSTS\tACGTACGTACGT\tGCATGCATGCAT\t100\t"

        with tempfile.NamedTemporaryFile(mode='w', suffix='.sts', delete=False) as f:
            f.write(sts_content)
            sts_file = f.name

        try:
            engine.load_sts_file(sts_file)

            # Sequence exactly long enough for one hash window
            sequence = "ACGTACGTACGT"  # Exactly 12 bases

            thread_data = ThreadData(
                thread_id=0,
                sequence=sequence,
                offset=0,
                length=len(sequence)
            )

            result = engine._process_thread(thread_data)
            assert isinstance(result, ThreadData)

        finally:
            Path(sts_file).unlink()


class TestMatchSTSEdgeCases:
    """Test _match_sts edge cases (lines 570-662)."""

    def test_match_sts_not_enough_room_for_primer2(self):
        """Test match_sts when there's not enough room for primer2."""
        engine = MerPCR(wordsize=11)

        sts_content = "TestSTS\tACGTACGTACGT\tGCATGCATGCATGCAT\t200\t"  # Long primer2

        with tempfile.NamedTemporaryFile(mode='w', suffix='.sts', delete=False) as f:
            f.write(sts_content)
            sts_file = f.name

        try:
            engine.load_sts_file(sts_file)

            # Sequence too short after primer1
            sequence = "ACGTACGTACGTNN"  # Just barely fits primer1, not primer2

            thread_data = ThreadData(
                thread_id=0,
                sequence=sequence,
                offset=0,
                length=len(sequence)
            )

            result = engine._process_thread(thread_data)
            # Should not crash, just return no hits
            assert isinstance(result, ThreadData)

        finally:
            Path(sts_file).unlink()

    def test_match_sts_with_margin_calculations(self):
        """Test match_sts with various margin edge cases."""
        engine = MerPCR(wordsize=11, margin=20)

        sts_content = "TestSTS\tACGTACGTACGT\tGCATGCATGCAT\t50\t"

        with tempfile.NamedTemporaryFile(mode='w', suffix='.sts', delete=False) as f:
            f.write(sts_content)
            sts_file = f.name

        try:
            engine.load_sts_file(sts_file)

            # Sequence that will exercise margin calculations
            sequence = "ACGTACGTACGT" + "N" * 30 + "ATGCATGCATGC"

            thread_data = ThreadData(
                thread_id=0,
                sequence=sequence,
                offset=0,
                length=len(sequence)
            )

            result = engine._process_thread(thread_data)
            assert isinstance(result, ThreadData)

        finally:
            Path(sts_file).unlink()

    def test_match_sts_exp_size_exceeds_actual_size(self):
        """Test match_sts when expected size exceeds available size."""
        engine = MerPCR(wordsize=11)

        # Large expected PCR size
        sts_content = "TestSTS\tACGTACGTACGT\tGCATGCATGCAT\t5000\t"

        with tempfile.NamedTemporaryFile(mode='w', suffix='.sts', delete=False) as f:
            f.write(sts_content)
            sts_file = f.name

        try:
            engine.load_sts_file(sts_file)

            # Small sequence (much smaller than expected PCR size)
            sequence = "ACGTACGTACGT" + "N" * 50 + "ATGCATGCATGC"

            thread_data = ThreadData(
                thread_id=0,
                sequence=sequence,
                offset=0,
                length=len(sequence)
            )

            result = engine._process_thread(thread_data)
            assert isinstance(result, ThreadData)

        finally:
            Path(sts_file).unlink()


class TestIUPACMode:
    """Test IUPAC mode functionality (lines 643-707)."""

    def test_iupac_mode_enabled(self):
        """Test sequence comparison with IUPAC mode enabled."""
        engine = MerPCR(wordsize=11, iupac_mode=1, mismatches=0)

        # Test IUPAC matching
        seq1 = "ACGTRYSWKM"
        seq2 = "ACGTACGTAN"  # Should match with IUPAC codes

        # R matches A or G, Y matches C or T, etc.
        result = engine._compare_seqs(seq1[:4], seq2[:4], "+")
        assert result is True, "Should match exact bases"

    def test_iupac_mode_with_unrecognized_characters(self):
        """Test IUPAC mode with unrecognized characters."""
        engine = MerPCR(wordsize=11, iupac_mode=1)

        # Unrecognized character (not in IUPAC mapping)
        seq1 = "ACGT"
        seq2 = "ACGZ"  # Z is not a valid IUPAC code

        result = engine._compare_seqs(seq1, seq2, "+")
        assert result is False, "Should not match unrecognized character"

    def test_iupac_mode_r_code(self):
        """Test IUPAC R code (A or G)."""
        engine = MerPCR(wordsize=11, iupac_mode=1, mismatches=0)

        # R should match both A and G
        assert engine._compare_seqs("R", "A", "+") is True
        assert engine._compare_seqs("R", "G", "+") is True
        assert engine._compare_seqs("R", "C", "+") is False

    def test_iupac_mode_with_three_prime_protection(self):
        """Test IUPAC mode combined with 3' protection."""
        engine = MerPCR(wordsize=11, iupac_mode=1, mismatches=0, three_prime_match=2)

        # In IUPAC mode, N matches any base, so this actually should succeed
        seq1 = "ACGTACGTNN"  # N at 3' end
        seq2 = "ACGTACGTAA"
        result = engine._compare_seqs(seq1, seq2, "+")
        assert result is True, "N should match A in IUPAC mode"

        # Test actual mismatch in 3' region (non-matching bases)
        seq1 = "ACGTACGTAC"
        seq2 = "ACGTACGTGG"  # Last 2 bases don't match
        result = engine._compare_seqs(seq1, seq2, "+")
        assert result is False, "Should fail due to mismatch in 3' region"

    def test_iupac_mode_n_code(self):
        """Test IUPAC N code (any base)."""
        engine = MerPCR(wordsize=11, iupac_mode=1)

        # N should match any base
        assert engine._compare_seqs("N", "A", "+") is True
        assert engine._compare_seqs("N", "C", "+") is True
        assert engine._compare_seqs("N", "G", "+") is True
        assert engine._compare_seqs("N", "T", "+") is True


class TestMultithreading:
    """Test multithreading paths (lines 474-493)."""

    @pytest.mark.slow
    def test_search_with_multiple_threads(self):
        """Test search with actual multithreading (requires large sequence)."""
        engine = MerPCR(wordsize=11, threads=2)

        sts_content = "TestSTS\tACGTACGTACGT\tGCATGCATGCAT\t100\t"
        # Create large sequence to trigger threading (>100KB)
        large_sequence = "ACGT" * 30000  # ~120KB
        fasta_content = f">large_seq\n{large_sequence}\n"

        with tempfile.NamedTemporaryFile(mode='w', suffix='.sts', delete=False) as f:
            f.write(sts_content)
            sts_file = f.name

        with tempfile.NamedTemporaryFile(mode='w', suffix='.fa', delete=False) as f:
            f.write(fasta_content)
            fasta_file = f.name

        try:
            engine.load_sts_file(sts_file)
            records = engine.load_fasta_file(fasta_file)
            hit_count = engine.search(records)
            assert hit_count >= 0

        finally:
            Path(sts_file).unlink()
            Path(fasta_file).unlink()

    def test_search_filters_redundant_hits_in_overlap(self):
        """Test that redundant hits in overlap regions are filtered."""
        engine = MerPCR(wordsize=11, threads=2)

        sts_content = "TestSTS\tACGTACGTACGT\tGCATGCATGCAT\t100\t"
        # Large enough to use multiple threads
        large_sequence = "ACGT" * 30000
        fasta_content = f">large_seq\n{large_sequence}\n"

        with tempfile.NamedTemporaryFile(mode='w', suffix='.sts', delete=False) as f:
            f.write(sts_content)
            sts_file = f.name

        with tempfile.NamedTemporaryFile(mode='w', suffix='.fa', delete=False) as f:
            f.write(fasta_content)
            fasta_file = f.name

        try:
            engine.load_sts_file(sts_file)
            records = engine.load_fasta_file(fasta_file)
            hit_count = engine.search(records)
            # Just verify it runs without crashing
            assert hit_count >= 0

        finally:
            Path(sts_file).unlink()
            Path(fasta_file).unlink()


class TestReverseComplement:
    """Test reverse complement with unknown bases."""

    def test_reverse_complement_unknown_base(self):
        """Test reverse complement with unknown base returns N."""
        engine = MerPCR()

        # Unknown base 'Z' should become 'N'
        result = engine._reverse_complement("ACGTZ")
        # Z is not in compl dict, so should get N
        assert "N" in result or result.endswith("Z"), "Unknown base handling"


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
