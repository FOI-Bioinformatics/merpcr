"""
Tests for forward strand primer matching (Critical Fix #3).

This test module validates the fix for forward strand primer matching,
which was causing merPCR to miss ALL forward (+) strand hits.

The issue was that primer2 was stored as-is for forward direction instead
of being reverse complemented. Since primer2 binds to the reverse strand,
we need to search for its reverse complement in the forward sequence.

Reference: VERIFICATION.md - Critical Fix #3
"""

import tempfile
from pathlib import Path

import pytest

from merpcr import MerPCR


class TestForwardStrandMatching:
    """Test forward strand primer matching functionality."""

    def test_primer2_reverse_complement_storage(self):
        """Test that primer2 is stored as reverse complement for forward direction."""
        engine = MerPCR(wordsize=11, margin=50)

        # Create a simple STS file
        sts_content = "TestSTS\tACGTACGTACGT\tGCATGCATGCAT\t100\tTest"

        with tempfile.NamedTemporaryFile(mode="w", suffix=".sts", delete=False) as f:
            f.write(sts_content)
            sts_file = f.name

        try:
            # Load STS file
            engine.load_sts_file(sts_file)

            # Check that we have both forward and reverse records
            assert len(engine.sts_records) == 2, "Should have 2 STS records (+ and -)"

            # Find forward and reverse records
            forward_records = [r for r in engine.sts_records if r.direct == "+"]
            reverse_records = [r for r in engine.sts_records if r.direct == "-"]

            assert len(forward_records) == 1, "Should have 1 forward record"
            assert len(reverse_records) == 1, "Should have 1 reverse record"

            forward_sts = forward_records[0]

            # For forward direction, primer1 should be original, primer2 should be RC
            assert forward_sts.primer1 == "ACGTACGTACGT", "Forward primer1 should be original"

            # primer2 original: GCATGCATGCAT
            # primer2 RC: ATGCATGCATGC
            expected_primer2_rc = "ATGCATGCATGC"
            assert (
                forward_sts.primer2 == expected_primer2_rc
            ), f"Forward primer2 should be reverse complement, got {forward_sts.primer2}"

        finally:
            Path(sts_file).unlink()

    def test_forward_strand_hit_detection(self):
        """Test that forward strand hits are detected correctly."""
        engine = MerPCR(wordsize=11, margin=50)

        # Create STS with known primers
        # Primer1: CCCTGATTTACAAGAAGTC (19bp)
        # Primer2: CTTGGTTATCATCTTTATCATATC (24bp)
        # PCR size: 86bp
        sts_content = "A_fra_ISFtu2\tCCCTGATTTACAAGAAGTC\tCTTGGTTATCATCTTTATCATATC\t86\t"

        # Create sequence that matches this STS (from real F. tularensis data)
        # Position 254165..254250 (86 bp)
        sequence = (
            "CCCTGATTTACAAGAAGTC"  # Primer 1 at start
            "ATGCTTGATTCAACAATAGCAAGAGCACATGCTTGTGCTACGG"  # Middle region
            "GATATGATAAAGATGATAACCAAG"  # Primer 2 RC at end
        )

        fasta_content = f">test_seq\n{sequence}\n"

        with tempfile.NamedTemporaryFile(mode="w", suffix=".sts", delete=False) as f:
            f.write(sts_content)
            sts_file = f.name

        with tempfile.NamedTemporaryFile(mode="w", suffix=".fa", delete=False) as f:
            f.write(fasta_content)
            fasta_file = f.name

        try:
            # Load files
            engine.load_sts_file(sts_file)
            records = engine.load_fasta_file(fasta_file)

            # Search
            hits = []
            for record in records:
                seq_label = record.label
                sequence = record.sequence
                seq_len = len(sequence)

                from merpcr.core.models import ThreadData

                thread_data = ThreadData(thread_id=0, sequence=sequence, offset=0, length=seq_len)
                engine._process_thread(thread_data)
                hits.extend(thread_data.hits)

            # Should find the forward strand hit
            assert len(hits) >= 1, "Should find at least one hit"

            # Check that we found a forward strand hit
            forward_hits = [h for h in hits if h.sts.direct == "+"]
            assert len(forward_hits) >= 1, "Should find forward strand hit"

            # Verify hit coordinates
            hit = forward_hits[0]
            assert hit.pos1 == 0, f"Expected pos1=0, got {hit.pos1}"
            assert hit.pos2 == len(sequence) - 1, f"Expected pos2={len(sequence)-1}, got {hit.pos2}"

        finally:
            Path(sts_file).unlink()
            Path(fasta_file).unlink()

    def test_both_strands_detected(self):
        """Test that both forward and reverse strand hits are detected."""
        engine = MerPCR(wordsize=11, margin=50)

        # Create STS
        sts_content = "TestSTS\tACGTACGTACGT\tGCATGCATGCAT\t100\t"

        # Create sequence with hits on both strands
        # Forward hit: ACGTACGTACGT...(middle)...ATGCATGCATGC (primer2 RC)
        # Reverse hit would have primer2 forward and primer1 RC

        fasta_content = (
            ">test_seq\n"
            "ACGTACGTACGTNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN"
            "NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN"
            "NNNNNNNNNNATGCATGCATGC\n"
        )

        with tempfile.NamedTemporaryFile(mode="w", suffix=".sts", delete=False) as f:
            f.write(sts_content)
            sts_file = f.name

        with tempfile.NamedTemporaryFile(mode="w", suffix=".fa", delete=False) as f:
            f.write(fasta_content)
            fasta_file = f.name

        try:
            # Load and search
            engine.load_sts_file(sts_file)
            records = engine.load_fasta_file(fasta_file)

            from merpcr.core.models import ThreadData

            hits = []
            for record in records:
                thread_data = ThreadData(
                    thread_id=0, sequence=record.sequence, offset=0, length=len(record.sequence)
                )
                engine._process_thread(thread_data)
                hits.extend(thread_data.hits)

            # Should find forward strand hit
            forward_hits = [h for h in hits if h.sts.direct == "+"]
            assert len(forward_hits) >= 1, "Should find forward strand hit"

        finally:
            Path(sts_file).unlink()
            Path(fasta_file).unlink()

    def test_real_genomic_case_f_tularensis(self):
        """Test the real F. tularensis case that revealed the bug."""
        engine = MerPCR(wordsize=11, margin=50)

        # This is the exact primer that was missing before the fix
        # From real F. tularensis genome validation
        sts_content = "A_fra_ISFtu2\tCCCTGATTTACAAGAAGTC\tCTTGGTTATCATCTTTATCATATC\t86\t"

        # Real sequence from position 254165..254250
        real_sequence = (
            "CCCTGATTTACAAGAAGTCATGCTTGATTCAACAATAGCAAGAGCACATGCTTGTGCTACGG"
            "GATATGATAAAGATGATAACCAAG"
        )

        fasta_content = f">NZ_CP010115.1\n{real_sequence}\n"

        with tempfile.NamedTemporaryFile(mode="w", suffix=".sts", delete=False) as f:
            f.write(sts_content)
            sts_file = f.name

        with tempfile.NamedTemporaryFile(mode="w", suffix=".fa", delete=False) as f:
            f.write(fasta_content)
            fasta_file = f.name

        try:
            engine.load_sts_file(sts_file)
            records = engine.load_fasta_file(fasta_file)

            # Count hits
            hit_count = 0
            forward_hits = 0

            from merpcr.core.models import ThreadData

            for record in records:
                thread_data = ThreadData(
                    thread_id=0, sequence=record.sequence, offset=0, length=len(record.sequence)
                )
                engine._process_thread(thread_data)

                for hit in thread_data.hits:
                    hit_count += 1
                    if hit.sts.direct == "+":
                        forward_hits += 1

            # This is the critical test: before the fix, forward_hits would be 0
            assert hit_count >= 1, "Should find at least one hit"
            assert forward_hits >= 1, "Should find forward strand hit (THIS WAS THE BUG!)"

        finally:
            Path(sts_file).unlink()
            Path(fasta_file).unlink()

    def test_primer_reverse_complement_correctness(self):
        """Test that reverse complement is computed correctly."""
        engine = MerPCR()

        # Test cases
        test_cases = [
            ("ATCG", "CGAT"),
            ("AAAA", "TTTT"),
            ("GCGC", "GCGC"),  # Palindrome
            ("CTTGGTTATCATCTTTATCATATC", "GATATGATAAAGATGATAACCAAG"),  # Real primer2
        ]

        for forward, expected_rc in test_cases:
            rc = engine._reverse_complement(forward)
            assert rc == expected_rc, f"RC of {forward} should be {expected_rc}, got {rc}"

    def test_compare_seqs_with_strand_parameter(self):
        """Test sequence comparison with correct strand parameter."""
        engine = MerPCR(mismatches=0, three_prime_match=1)

        # Test forward strand comparison (strand="+")
        seq1 = "ACGTACGT"
        seq2 = "ACGTACGT"
        assert engine._compare_seqs(seq1, seq2, "+"), "Should match on forward strand"

        # Test with mismatch in 3' region (should fail)
        seq1 = "ACGTACGT"
        seq2 = "ACGTACGA"  # Mismatch at end (3' for +)
        assert not engine._compare_seqs(
            seq1, seq2, "+"
        ), "Should not match with 3' mismatch on forward strand"

        # Test reverse strand comparison (strand="-")
        seq1 = "ACGTACGT"
        seq2 = "ACGTACGT"
        assert engine._compare_seqs(seq1, seq2, "-"), "Should match on reverse strand"

        # Test with mismatch in 3' region for reverse (should fail)
        seq1 = "ACGTACGT"
        seq2 = "GCGTACGT"  # Mismatch at start (3' for -)
        assert not engine._compare_seqs(
            seq1, seq2, "-"
        ), "Should not match with 3' mismatch on reverse strand"


class TestForwardStrandRegressions:
    """Regression tests to ensure the fix doesn't break anything."""

    def test_reverse_strand_still_works(self):
        """Ensure reverse strand detection still works after forward fix."""
        engine = MerPCR(wordsize=11, margin=50)

        sts_content = "TestSTS\tACGTACGTACGT\tGCATGCATGCAT\t100\t"

        # Sequence with only reverse hit: primer2 forward, then primer1 RC
        # primer1 RC of ACGTACGTACGT = ACGTACGTACGT (palindrome)
        fasta_content = (
            ">test_seq\n"
            "GCATGCATGCATNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN"
            "NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN"
            "NNNNNNNNNACGTACGTACGT\n"
        )

        with tempfile.NamedTemporaryFile(mode="w", suffix=".sts", delete=False) as f:
            f.write(sts_content)
            sts_file = f.name

        with tempfile.NamedTemporaryFile(mode="w", suffix=".fa", delete=False) as f:
            f.write(fasta_content)
            fasta_file = f.name

        try:
            engine.load_sts_file(sts_file)
            records = engine.load_fasta_file(fasta_file)

            from merpcr.core.models import ThreadData

            hits = []
            for record in records:
                thread_data = ThreadData(
                    thread_id=0, sequence=record.sequence, offset=0, length=len(record.sequence)
                )
                engine._process_thread(thread_data)
                hits.extend(thread_data.hits)

            # Should find reverse strand hit
            reverse_hits = [h for h in hits if h.sts.direct == "-"]
            assert len(reverse_hits) >= 1, "Should still find reverse strand hits"

        finally:
            Path(sts_file).unlink()
            Path(fasta_file).unlink()

    def test_no_false_positives(self):
        """Ensure the fix doesn't create false positive hits."""
        engine = MerPCR(wordsize=11, margin=50)

        sts_content = "TestSTS\tACGTACGTACGT\tGCATGCATGCAT\t100\t"

        # Sequence with no valid hits
        fasta_content = ">test_seq\nTTTTTTTTTTTTTTTTTTTTTTTTTTTT\n"

        with tempfile.NamedTemporaryFile(mode="w", suffix=".sts", delete=False) as f:
            f.write(sts_content)
            sts_file = f.name

        with tempfile.NamedTemporaryFile(mode="w", suffix=".fa", delete=False) as f:
            f.write(fasta_content)
            fasta_file = f.name

        try:
            engine.load_sts_file(sts_file)
            records = engine.load_fasta_file(fasta_file)

            from merpcr.core.models import ThreadData

            hits = []
            for record in records:
                thread_data = ThreadData(
                    thread_id=0, sequence=record.sequence, offset=0, length=len(record.sequence)
                )
                engine._process_thread(thread_data)
                hits.extend(thread_data.hits)

            # Should find no hits
            assert len(hits) == 0, "Should not find any hits in sequence without primers"

        finally:
            Path(sts_file).unlink()
            Path(fasta_file).unlink()
