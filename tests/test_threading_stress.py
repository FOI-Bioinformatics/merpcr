"""
Threading and concurrency stress tests.

These tests verify that the multi-threading implementation is robust
and handles various threading scenarios correctly.
"""

import logging
import os
import random
import tempfile
import threading
import time
from concurrent.futures import ThreadPoolExecutor, as_completed

import pytest

from merpcr.core.engine import MerPCR
from merpcr.core.models import FASTARecord

# Disable logging for stress tests to reduce noise
logging.getLogger("merpcr").setLevel(logging.ERROR)


class TestThreadingBehavior:
    """Test threading behavior and thread safety."""

    def test_single_vs_multi_thread_consistency(self):
        """Test that single-threaded and multi-threaded results are consistent."""
        # Create test data
        sts_content = """TEST1\tATCGATCGATCG\tGCTAGCTAGCTA\t100
TEST2\tTACGTACGTACG\tCGTACGTACGTA\t120
TEST3\tGGGGAAAATTTT\tCCCCTTTTAAAA\t150"""

        # Large sequence to force threading
        large_sequence = (
            "ATCGATCGATCG"
            + "N" * 200
            + "GCTAGCTAGCTA"
            + "N" * 1000
            + "TACGTACGTACG"
            + "N" * 200
            + "CGTACGTACGTA"
            + "N" * 1000
            + "GGGGAAAATTTT"
            + "N" * 200
            + "CCCCTTTTAAAA"
        ) * 100  # ~350KB

        fasta_content = f">large_test\n{large_sequence}"

        with tempfile.NamedTemporaryFile(mode="w", suffix=".sts", delete=False) as sts_f:
            sts_f.write(sts_content)
            sts_path = sts_f.name

        with tempfile.NamedTemporaryFile(mode="w", suffix=".fa", delete=False) as fasta_f:
            fasta_f.write(fasta_content)
            fasta_path = fasta_f.name

        try:
            # Single-threaded run
            engine1 = MerPCR(wordsize=8, threads=1)
            engine1.load_sts_file(sts_path)
            records1 = engine1.load_fasta_file(fasta_path)
            hits1 = engine1.search(records1)

            # Multi-threaded run
            engine2 = MerPCR(wordsize=8, threads=4)
            engine2.load_sts_file(sts_path)
            records2 = engine2.load_fasta_file(fasta_path)
            hits2 = engine2.search(records2)

            # Results should be identical
            assert hits1 == hits2, f"Single-threaded: {hits1}, Multi-threaded: {hits2}"

        finally:
            os.unlink(sts_path)
            os.unlink(fasta_path)

    def test_thread_count_scaling(self):
        """Test behavior with different thread counts."""
        sts_content = "TEST\tATCGATCGATCG\tGCTAGCTAGCTA\t100"
        # Medium-sized sequence
        sequence = ("ATCGATCGATCG" + "N" * 1000 + "GCTAGCTAGCTA") * 50  # ~50KB
        fasta_content = f">test\n{sequence}"

        with tempfile.NamedTemporaryFile(mode="w", suffix=".sts", delete=False) as sts_f:
            sts_f.write(sts_content)
            sts_path = sts_f.name

        with tempfile.NamedTemporaryFile(mode="w", suffix=".fa", delete=False) as fasta_f:
            fasta_f.write(fasta_content)
            fasta_path = fasta_f.name

        try:
            results = {}
            thread_counts = [1, 2, 4, 8]

            for threads in thread_counts:
                engine = MerPCR(wordsize=8, threads=threads)
                engine.load_sts_file(sts_path)
                records = engine.load_fasta_file(fasta_path)

                start_time = time.time()
                hits = engine.search(records)
                end_time = time.time()

                results[threads] = {"hits": hits, "time": end_time - start_time}

            # All thread counts should produce same results
            hit_counts = [r["hits"] for r in results.values()]
            assert all(h == hit_counts[0] for h in hit_counts), f"Inconsistent results: {results}"

            # More threads should generally not be slower (though this can vary)
            # At minimum, they should all complete successfully
            assert all(r["time"] > 0 for r in results.values())

        finally:
            os.unlink(sts_path)
            os.unlink(fasta_path)

    def test_concurrent_merpcr_instances(self):
        """Test running multiple MerPCR instances concurrently."""
        num_instances = 4

        def run_merpcr_instance(instance_id):
            """Run a single MerPCR instance."""
            sts_content = f"TEST{instance_id}\tATCGATCGATCG\tGCTAGCTAGCTA\t100"
            sequence = ("ATCGATCGATCG" + "N" * 500 + "GCTAGCTAGCTA") * 20
            fasta_content = f">test{instance_id}\n{sequence}"

            with tempfile.NamedTemporaryFile(mode="w", suffix=".sts", delete=False) as sts_f:
                sts_f.write(sts_content)
                sts_path = sts_f.name

            with tempfile.NamedTemporaryFile(mode="w", suffix=".fa", delete=False) as fasta_f:
                fasta_f.write(fasta_content)
                fasta_path = fasta_f.name

            try:
                engine = MerPCR(wordsize=8, threads=2)
                engine.load_sts_file(sts_path)
                records = engine.load_fasta_file(fasta_path)
                hits = engine.search(records)

                return {"instance": instance_id, "hits": hits, "success": True}

            except Exception as e:
                return {"instance": instance_id, "error": str(e), "success": False}
            finally:
                os.unlink(sts_path)
                os.unlink(fasta_path)

        # Run instances concurrently
        with ThreadPoolExecutor(max_workers=num_instances) as executor:
            futures = [executor.submit(run_merpcr_instance, i) for i in range(num_instances)]
            results = [future.result() for future in as_completed(futures)]

        # All instances should complete successfully
        assert len(results) == num_instances
        successful = [r for r in results if r["success"]]
        assert len(successful) == num_instances, f"Some instances failed: {results}"

    def test_thread_safety_shared_data(self):
        """Test thread safety when multiple threads access shared data structures."""
        # This test checks if there are any race conditions in shared data access

        sts_content = """TEST1\tATCGATCGATCG\tGCTAGCTAGCTA\t100
TEST2\tTACGTACGTACG\tCGTACGTACGTA\t120
TEST3\tGGGGAAAATTTT\tCCCCTTTTAAAA\t150"""

        sequence = (
            "ATCGATCGATCG"
            + "GCTAGCTAGCTA"
            + "TACGTACGTACG"
            + "CGTACGTACGTA"
            + "GGGGAAAATTTT"
            + "CCCCTTTTAAAA"
        ) * 1000  # Large sequence

        fasta_content = f">shared_test\n{sequence}"

        with tempfile.NamedTemporaryFile(mode="w", suffix=".sts", delete=False) as sts_f:
            sts_f.write(sts_content)
            sts_path = sts_f.name

        with tempfile.NamedTemporaryFile(mode="w", suffix=".fa", delete=False) as fasta_f:
            fasta_f.write(fasta_content)
            fasta_path = fasta_f.name

        try:
            # Run the same search multiple times concurrently
            def run_search():
                engine = MerPCR(wordsize=8, threads=4)
                engine.load_sts_file(sts_path)
                records = engine.load_fasta_file(fasta_path)
                return engine.search(records)

            num_concurrent = 6
            with ThreadPoolExecutor(max_workers=num_concurrent) as executor:
                futures = [executor.submit(run_search) for _ in range(num_concurrent)]
                results = [future.result() for future in as_completed(futures)]

            # All runs should produce the same result
            assert len(set(results)) == 1, f"Inconsistent results from concurrent runs: {results}"

        finally:
            os.unlink(sts_path)
            os.unlink(fasta_path)


@pytest.mark.slow
class TestStressTesting:
    """Stress tests that push the system to its limits."""

    def test_memory_pressure_large_files(self):
        """Test behavior under memory pressure with large files."""
        # Generate a large STS dataset
        large_sts_lines = []
        for i in range(1000):  # 1000 STS entries
            primer1 = "".join(random.choices("ATCG", k=20))
            primer2 = "".join(random.choices("ATCG", k=20))
            pcr_size = random.randint(100, 500)
            large_sts_lines.append(f"STS{i:04d}\t{primer1}\t{primer2}\t{pcr_size}")

        # Large FASTA sequence (~1MB)
        large_sequence = "".join(random.choices("ATCG", k=250000))
        fasta_content = f">large_memory_test\n{large_sequence}"

        with tempfile.NamedTemporaryFile(mode="w", suffix=".sts", delete=False) as sts_f:
            sts_f.write("\n".join(large_sts_lines))
            sts_path = sts_f.name

        with tempfile.NamedTemporaryFile(mode="w", suffix=".fa", delete=False) as fasta_f:
            fasta_f.write(fasta_content)
            fasta_path = fasta_f.name

        try:
            # Test with maximum threading
            engine = MerPCR(wordsize=10, threads=8, margin=200)

            # Should be able to load large datasets
            success = engine.load_sts_file(sts_path)
            assert success, "Failed to load large STS file"

            records = engine.load_fasta_file(fasta_path)
            assert len(records) == 1, "Failed to load large FASTA file"

            # Search should complete without running out of memory
            start_time = time.time()
            hit_count = engine.search(records)
            end_time = time.time()

            # Basic success criteria
            assert isinstance(hit_count, int)
            assert hit_count >= 0
            assert end_time - start_time < 300  # Should complete within 5 minutes

        finally:
            os.unlink(sts_path)
            os.unlink(fasta_path)

    def test_rapid_consecutive_searches(self):
        """Test rapid consecutive searches for stability."""
        sts_content = "TEST\tATCGATCGATCG\tGCTAGCTAGCTA\t100"
        sequence = ("ATCGATCGATCG" + "N" * 1000 + "GCTAGCTAGCTA") * 100
        fasta_content = f">rapid_test\n{sequence}"

        with tempfile.NamedTemporaryFile(mode="w", suffix=".sts", delete=False) as sts_f:
            sts_f.write(sts_content)
            sts_path = sts_f.name

        with tempfile.NamedTemporaryFile(mode="w", suffix=".fa", delete=False) as fasta_f:
            fasta_f.write(fasta_content)
            fasta_path = fasta_f.name

        try:
            results = []
            num_searches = 20

            for i in range(num_searches):
                engine = MerPCR(wordsize=8, threads=4)
                engine.load_sts_file(sts_path)
                records = engine.load_fasta_file(fasta_path)
                hits = engine.search(records)
                results.append(hits)

            # All searches should produce consistent results
            assert len(set(results)) == 1, f"Inconsistent results in rapid searches: {set(results)}"

        finally:
            os.unlink(sts_path)
            os.unlink(fasta_path)

    def test_thread_pool_exhaustion(self):
        """Test behavior when thread pool is exhausted."""
        sts_content = "TEST\tATCGATCGATCG\tGCTAGCTAGCTA\t100"
        # Small sequence that forces single-threading
        sequence = "ATCGATCGATCG" + "N" * 100 + "GCTAGCTAGCTA"  # Small sequence
        fasta_content = f">small_test\n{sequence}"

        with tempfile.NamedTemporaryFile(mode="w", suffix=".sts", delete=False) as sts_f:
            sts_f.write(sts_content)
            sts_path = sts_f.name

        with tempfile.NamedTemporaryFile(mode="w", suffix=".fa", delete=False) as fasta_f:
            fasta_f.write(fasta_content)
            fasta_path = fasta_f.name

        try:
            # Request many threads for a small file (should be limited to 1)
            engine = MerPCR(wordsize=8, threads=100)  # Excessive thread count
            engine.load_sts_file(sts_path)
            records = engine.load_fasta_file(fasta_path)

            # Should handle thread pool management gracefully
            hits = engine.search(records)
            assert isinstance(hits, int)
            assert hits >= 0

        finally:
            os.unlink(sts_path)
            os.unlink(fasta_path)


class TestThreadingEdgeCases:
    """Test threading edge cases and error conditions."""

    def test_threading_with_zero_hits(self):
        """Test threading behavior when no hits are found."""
        sts_content = "TEST\tAAAAAAAAAAAA\tTTTTTTTTTTTT\t100"  # Won't match
        sequence = "CCCCGGGGCCCCGGGG" * 1000  # No A's or T's
        fasta_content = f">no_hits\n{sequence}"

        with tempfile.NamedTemporaryFile(mode="w", suffix=".sts", delete=False) as sts_f:
            sts_f.write(sts_content)
            sts_path = sts_f.name

        with tempfile.NamedTemporaryFile(mode="w", suffix=".fa", delete=False) as fasta_f:
            fasta_f.write(fasta_content)
            fasta_path = fasta_f.name

        try:
            engine = MerPCR(wordsize=8, threads=4)
            engine.load_sts_file(sts_path)
            records = engine.load_fasta_file(fasta_path)
            hits = engine.search(records)

            assert hits == 0
            assert engine.total_hits == 0

        finally:
            os.unlink(sts_path)
            os.unlink(fasta_path)

    def test_threading_with_many_hits(self):
        """Test threading behavior with many expected hits."""
        # STS that will match many times
        sts_content = "REPEAT\tATCG\tCGAT\t20"
        # Sequence with many matches
        sequence = "ATCGCGAT" * 1000  # Should produce many hits
        fasta_content = f">many_hits\n{sequence}"

        with tempfile.NamedTemporaryFile(mode="w", suffix=".sts", delete=False) as sts_f:
            sts_f.write(sts_content)
            sts_path = sts_f.name

        with tempfile.NamedTemporaryFile(mode="w", suffix=".fa", delete=False) as fasta_f:
            fasta_f.write(fasta_content)
            fasta_path = fasta_f.name

        try:
            # Single-threaded
            engine1 = MerPCR(wordsize=4, threads=1, margin=50)
            engine1.load_sts_file(sts_path)
            records1 = engine1.load_fasta_file(fasta_path)
            hits1 = engine1.search(records1)

            # Multi-threaded
            engine2 = MerPCR(wordsize=4, threads=4, margin=50)
            engine2.load_sts_file(sts_path)
            records2 = engine2.load_fasta_file(fasta_path)
            hits2 = engine2.search(records2)

            # Should find the same number of hits
            assert hits1 == hits2, f"Hit count mismatch: single={hits1}, multi={hits2}"
            assert hits1 > 0  # Should find some hits

        finally:
            os.unlink(sts_path)
            os.unlink(fasta_path)

    def test_interruption_resilience(self):
        """Test resilience to interruption (basic version)."""
        sts_content = "TEST\tATCGATCGATCG\tGCTAGCTAGCTA\t100"
        sequence = ("ATCGATCGATCG" + "N" * 2000 + "GCTAGCTAGCTA") * 200
        fasta_content = f">interrupt_test\n{sequence}"

        with tempfile.NamedTemporaryFile(mode="w", suffix=".sts", delete=False) as sts_f:
            sts_f.write(sts_content)
            sts_path = sts_f.name

        with tempfile.NamedTemporaryFile(mode="w", suffix=".fa", delete=False) as fasta_f:
            fasta_f.write(fasta_content)
            fasta_path = fasta_f.name

        try:
            engine = MerPCR(wordsize=8, threads=4)
            engine.load_sts_file(sts_path)
            records = engine.load_fasta_file(fasta_path)

            # This test mainly ensures the search can complete
            # In a real interruption scenario, proper cleanup would be tested
            hits = engine.search(records)
            assert isinstance(hits, int)
            assert hits >= 0

        finally:
            os.unlink(sts_path)
            os.unlink(fasta_path)
