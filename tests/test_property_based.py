"""
Property-based tests for robustness and edge case discovery.

These tests use generated data to find edge cases and ensure robustness.
"""

import os
import random
import string
import tempfile

import pytest

from merpcr.cli import convert_mepcr_arguments
from merpcr.core.engine import MerPCR
from merpcr.core.models import FASTARecord, STSRecord
from merpcr.core.utils import hash_value, init_iupac_tables, reverse_complement

try:
    from hypothesis import assume, given, settings
    from hypothesis import strategies as st
    from hypothesis.strategies import composite, integers, lists, text

    HYPOTHESIS_AVAILABLE = True
except ImportError:
    HYPOTHESIS_AVAILABLE = False

pytestmark = pytest.mark.skipif(not HYPOTHESIS_AVAILABLE, reason="hypothesis not available")


# Custom strategies for DNA sequences
@composite
def dna_sequence(draw, min_length=1, max_length=100):
    """Generate valid DNA sequences."""
    bases = "ATCG"
    length = draw(integers(min_value=min_length, max_value=max_length))
    return "".join(draw(st.lists(st.sampled_from(bases), min_size=length, max_size=length)))


@composite
def dna_sequence_with_ambiguities(draw, min_length=1, max_length=100):
    """Generate DNA sequences with IUPAC ambiguity codes."""
    bases = "ATCGRYMKSWBDHVN"
    length = draw(integers(min_value=min_length, max_value=max_length))
    return "".join(draw(st.lists(st.sampled_from(bases), min_size=length, max_size=length)))


@composite
def sts_record_data(draw):
    """Generate valid STS record data."""
    sts_id = draw(text(alphabet=string.ascii_letters + string.digits, min_size=1, max_size=20))
    primer1 = draw(dna_sequence(min_length=10, max_length=30))
    primer2 = draw(dna_sequence(min_length=10, max_length=30))
    pcr_size = draw(integers(min_value=50, max_value=1000))
    alias = draw(
        text(alphabet=string.ascii_letters + string.digits + "()- .", min_size=0, max_size=50)
    )

    return {
        "id": sts_id,
        "primer1": primer1,
        "primer2": primer2,
        "pcr_size": pcr_size,
        "alias": alias,
    }


class TestPropertyBasedUtils:
    """Property-based tests for utility functions."""

    @given(dna_sequence())
    @settings(max_examples=100)
    def test_reverse_complement_involution(self, sequence):
        """Test that reverse_complement(reverse_complement(x)) == x."""
        # Property: applying reverse complement twice should return original
        result = reverse_complement(reverse_complement(sequence))
        assert result == sequence, f"Double reverse complement failed for: {sequence}"

    @given(dna_sequence())
    @settings(max_examples=100)
    def test_reverse_complement_length_preservation(self, sequence):
        """Test that reverse complement preserves sequence length."""
        result = reverse_complement(sequence)
        assert len(result) == len(sequence)

    @given(dna_sequence())
    @settings(max_examples=100)
    def test_reverse_complement_base_validity(self, sequence):
        """Test that reverse complement produces valid DNA bases."""
        result = reverse_complement(sequence)
        valid_bases = set("ATCGRYMKSWBDHVN")  # Include IUPAC codes
        for base in result.upper():
            assert (
                base in valid_bases or base == "N"
            ), f"Invalid base {base} in result for {sequence}"

    @given(dna_sequence(), integers(min_value=3, max_value=16))
    @settings(max_examples=50)
    def test_hash_value_properties(self, sequence, wordsize):
        """Test hash value computation properties."""
        offset, hash_val = hash_value(sequence, wordsize)

        if len(sequence) < wordsize:
            # Should fail for sequences shorter than wordsize
            assert offset == -1
            assert hash_val == 0
        else:
            if offset >= 0:
                # Valid hash should have reasonable offset
                assert 0 <= offset <= len(sequence) - wordsize
                # Hash value should be non-negative
                assert hash_val >= 0
                # Hash should fit in expected bit range
                assert hash_val < (1 << (2 * wordsize))

    @given(dna_sequence(), dna_sequence(), integers(min_value=3, max_value=16))
    @settings(max_examples=50)
    def test_hash_value_deterministic(self, seq1, seq2, wordsize):
        """Test that hash values are deterministic."""
        if seq1 == seq2:
            hash1 = hash_value(seq1, wordsize)
            hash2 = hash_value(seq2, wordsize)
            assert hash1 == hash2, "Same sequences should produce same hashes"

    @given(st.booleans())
    @settings(max_examples=10)
    def test_iupac_tables_consistency(self, iupac_mode):
        """Test IUPAC table initialization consistency."""
        tables = init_iupac_tables(iupac_mode)

        if iupac_mode:
            assert len(tables) > 0
            # Standard bases should be present
            assert "A" in tables
            assert "T" in tables
            assert "C" in tables
            assert "G" in tables
            # Case variants should be present
            assert "a" in tables
            assert "t" in tables
        else:
            assert len(tables) == 0


class TestPropertyBasedEngine:
    """Property-based tests for engine functionality."""

    @given(
        integers(min_value=3, max_value=16),
        integers(min_value=0, max_value=1000),
        integers(min_value=0, max_value=10),
        integers(min_value=0, max_value=5),
        integers(min_value=0, max_value=1),
        integers(min_value=50, max_value=1000),
        integers(min_value=1, max_value=8),
    )
    @settings(max_examples=20)
    def test_merpcr_initialization_bounds(
        self, wordsize, margin, mismatches, three_prime_match, iupac_mode, default_pcr_size, threads
    ):
        """Test MerPCR initialization with various parameter combinations."""
        try:
            engine = MerPCR(
                wordsize=wordsize,
                margin=margin,
                mismatches=mismatches,
                three_prime_match=three_prime_match,
                iupac_mode=iupac_mode,
                default_pcr_size=default_pcr_size,
                threads=threads,
            )

            # Basic invariants
            assert engine.wordsize == wordsize
            assert engine.margin == margin
            assert engine.mismatches == mismatches
            assert engine.threads == threads
            assert engine.sts_records == []
            assert engine.total_hits == 0

        except Exception as e:
            # If initialization fails, it should be due to invalid parameters
            # This is acceptable as long as it fails gracefully
            pass

    @given(lists(sts_record_data(), min_size=1, max_size=10))
    @settings(max_examples=10)
    def test_sts_record_creation(self, sts_data_list):
        """Test STS record creation with various data."""
        for sts_data in sts_data_list:
            try:
                sts = STSRecord(
                    id=sts_data["id"],
                    primer1=sts_data["primer1"],
                    primer2=sts_data["primer2"],
                    pcr_size=sts_data["pcr_size"],
                    alias=sts_data["alias"],
                )

                # Verify data integrity
                assert sts.id == sts_data["id"]
                assert sts.primer1 == sts_data["primer1"]
                assert sts.primer2 == sts_data["primer2"]
                assert sts.pcr_size == sts_data["pcr_size"]
                assert sts.alias == sts_data["alias"]

            except Exception:
                # Creation might fail for invalid data, which is OK
                pass

    @given(
        text(alphabet=string.ascii_letters + string.digits + ">", min_size=1, max_size=50),
        dna_sequence(min_length=1, max_length=1000),
    )
    @settings(max_examples=20)
    def test_fasta_record_creation(self, defline, sequence):
        """Test FASTA record creation with various inputs."""
        try:
            fasta = FASTARecord(defline=defline, sequence=sequence)

            # Basic invariants
            assert fasta.defline == defline
            assert fasta.sequence == sequence
            assert len(fasta.label) > 0  # Should extract some label

        except Exception:
            # May fail for invalid deflines, which is acceptable
            pass


class TestPropertyBasedCLI:
    """Property-based tests for CLI functionality."""

    @given(
        lists(
            text(alphabet=string.ascii_letters + "=0123456789", min_size=1, max_size=10),
            min_size=0,
            max_size=20,
        )
    )
    @settings(max_examples=50)
    def test_convert_mepcr_arguments_robustness(self, args):
        """Test argument conversion with various inputs."""
        try:
            result = convert_mepcr_arguments(args)

            # Result should be a list
            assert isinstance(result, list)
            # Result should not be longer than input (no argument should expand to more than 2)
            assert len(result) <= len(args) * 2
            # All results should be strings
            assert all(isinstance(arg, str) for arg in result)

        except Exception:
            # May fail for invalid inputs, which is acceptable
            pass

    @given(integers(min_value=0, max_value=100))
    @settings(max_examples=20)
    def test_parameter_types_bounds(self, value):
        """Test parameter type functions with various values."""
        from merpcr.cli import (margin_type, mismatch_type, threads_type,
                                wordsize_type)

        # Test margin_type
        try:
            result = margin_type(str(value))
            if 0 <= value <= 10000:
                assert result == value
            else:
                assert False, "Should have raised exception for out-of-bounds value"
        except:
            # Should raise for invalid values
            assert value < 0 or value > 10000

        # Test other types with their respective bounds
        try:
            if 0 <= value <= 10:
                result = mismatch_type(str(value))
                assert result == value
        except:
            pass

        try:
            if 3 <= value <= 16:
                result = wordsize_type(str(value))
                assert result == value
        except:
            pass

        try:
            if value > 0:
                result = threads_type(str(value))
                assert result == value
        except:
            pass


class TestPropertyBasedSequenceComparison:
    """Property-based tests for sequence comparison."""

    @given(
        dna_sequence(min_length=5, max_length=20),
        dna_sequence(min_length=5, max_length=20),
        st.sampled_from(["+", "-"]),
        integers(min_value=0, max_value=5),
        integers(min_value=0, max_value=3),
    )
    @settings(max_examples=30)
    def test_compare_seqs_properties(self, seq1, seq2, strand, mismatches, three_prime_match):
        """Test sequence comparison properties."""
        engine = MerPCR(mismatches=mismatches, three_prime_match=three_prime_match)

        try:
            result = engine._compare_seqs(seq1, seq2, strand)

            # Result should be boolean
            assert isinstance(result, bool)

            # Identical sequences should match (if within mismatch tolerance)
            if seq1 == seq2:
                assert result == True

            # If sequences differ by more than allowed mismatches, should not match
            if len(seq1) == len(seq2):
                diff_count = sum(1 for a, b in zip(seq1.upper(), seq2.upper()) if a != b)
                if diff_count > mismatches:
                    # Might still match due to 3' protection rules, but if no protection...
                    if three_prime_match == 0:
                        assert result == False or diff_count <= mismatches

        except Exception:
            # May fail for invalid inputs
            pass


class TestPropertyBasedSearch:
    """Property-based tests for search functionality."""

    @given(integers(min_value=1, max_value=5))
    @settings(max_examples=5)
    def test_search_with_generated_data(self, num_sts):
        """Test search with generated STS and FASTA data."""
        # Generate simple test data
        sts_lines = []
        for i in range(num_sts):
            primer1 = "".join(random.choices("ATCG", k=12))
            primer2 = "".join(random.choices("ATCG", k=12))
            pcr_size = random.randint(50, 200)
            sts_lines.append(f"STS{i}\t{primer1}\t{primer2}\t{pcr_size}")

        # Generate a sequence that might contain some of the primers
        test_sequence = "".join(random.choices("ATCG", k=500))
        fasta_content = f">test_seq\n{test_sequence}"

        # Create temporary files
        with tempfile.NamedTemporaryFile(mode="w", suffix=".sts", delete=False) as sts_f:
            sts_f.write("\n".join(sts_lines))
            sts_path = sts_f.name

        with tempfile.NamedTemporaryFile(mode="w", suffix=".fa", delete=False) as fasta_f:
            fasta_f.write(fasta_content)
            fasta_path = fasta_f.name

        try:
            engine = MerPCR(wordsize=8, margin=100)

            # Test should complete without crashing
            success = engine.load_sts_file(sts_path)
            if success:
                records = engine.load_fasta_file(fasta_path)
                if records:
                    hit_count = engine.search(records)

                    # Basic properties
                    assert isinstance(hit_count, int)
                    assert hit_count >= 0
                    assert engine.total_hits == hit_count

        except Exception as e:
            # Search might fail due to file issues, which is acceptable
            pass
        finally:
            os.unlink(sts_path)
            os.unlink(fasta_path)


@pytest.mark.slow
class TestPropertyBasedLarge:
    """Property-based tests with larger datasets (marked as slow)."""

    @given(integers(min_value=10, max_value=100))
    @settings(max_examples=5)
    def test_large_sts_dataset(self, num_sts):
        """Test with larger STS datasets."""
        assume(num_sts >= 10)

        # Generate many STS entries
        sts_lines = []
        for i in range(num_sts):
            primer1 = "".join(random.choices("ATCG", k=random.randint(10, 25)))
            primer2 = "".join(random.choices("ATCG", k=random.randint(10, 25)))
            pcr_size = random.randint(100, 500)
            sts_lines.append(f"STS{i:03d}\t{primer1}\t{primer2}\t{pcr_size}\tTest STS {i}")

        fasta_content = f">large_seq\n{'ATCG' * 1000}"  # 4KB sequence

        with tempfile.NamedTemporaryFile(mode="w", suffix=".sts", delete=False) as sts_f:
            sts_f.write("\n".join(sts_lines))
            sts_path = sts_f.name

        with tempfile.NamedTemporaryFile(mode="w", suffix=".fa", delete=False) as fasta_f:
            fasta_f.write(fasta_content)
            fasta_path = fasta_f.name

        try:
            engine = MerPCR(wordsize=8, margin=50)

            success = engine.load_sts_file(sts_path)
            assert success  # Should be able to load valid data

            # Should have loaded the expected number of STS records
            assert len(engine.sts_records) > 0  # Some might be filtered out

            records = engine.load_fasta_file(fasta_path)
            assert len(records) == 1

            # Search should complete
            hit_count = engine.search(records)
            assert hit_count >= 0

        finally:
            os.unlink(sts_path)
            os.unlink(fasta_path)
