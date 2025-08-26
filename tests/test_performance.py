#!/usr/bin/env python3
"""
Performance and benchmarking tests for merPCR.
"""

import unittest
import time
import tempfile
import os
import sys
from pathlib import Path
import pytest

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'src'))

from merpcr import MerPCR
from merpcr.core.models import FASTARecord


@pytest.mark.performance
class TestPerformance(unittest.TestCase):
    """Performance benchmarking tests."""
    
    def setUp(self):
        """Set up performance test environment."""
        self.mer_pcr = MerPCR(threads=1)  # Single thread for consistent timing
        
    def create_large_sequence(self, length=100000, pattern="ATCGATCG"):
        """Create a large test sequence."""
        repeats = length // len(pattern)
        remainder = length % len(pattern)
        return pattern * repeats + pattern[:remainder]
    
    def create_test_sts(self, num_sts=100):
        """Create multiple test STSs."""
        sts_content = []
        for i in range(num_sts):
            sts_id = f"STS_{i:03d}"
            primer1 = "ATCGATCGATCG"
            primer2 = "CGATCGATCGAT" 
            pcr_size = 200 + i  # Vary the PCR size
            alias = f"Test STS {i}"
            sts_content.append(f"{sts_id}\t{primer1}\t{primer2}\t{pcr_size}\t{alias}")
        
        return "\n".join(sts_content) + "\n"
    
    @unittest.skipIf(os.getenv('SKIP_PERFORMANCE_TESTS'), "Performance tests skipped")
    def test_large_sequence_processing(self):
        """Test processing of large sequences."""
        # Create a smaller sequence for CI/testing 
        seq_size = 1000 if os.getenv('CI') or os.getenv('GITHUB_ACTIONS') else 100000
        large_seq = self.create_large_sequence(seq_size)
        fasta_record = FASTARecord(
            defline=">large_sequence",
            sequence=large_seq
        )
        
        # Create test STS data that will actually match the repeating pattern
        sts_content = "TEST001\tATCGATCG\tATCGATCG\t50\tTest STS\n"  # Simple matching STS
        with tempfile.NamedTemporaryFile(mode='w', suffix='.sts', delete=False) as f:
            f.write(sts_content)
            sts_file = f.name
        
        try:
            # Time STS loading
            start_time = time.time()
            success = self.mer_pcr.load_sts_file(sts_file)
            sts_load_time = time.time() - start_time
            
            self.assertTrue(success)
            self.assertLess(sts_load_time, 1.0, "STS loading should be under 1 second")
            
            # Time sequence processing
            start_time = time.time()
            hit_count = self.mer_pcr.search([fasta_record])
            search_time = time.time() - start_time
            
            expected_time = 2.0 if os.getenv('CI') or os.getenv('GITHUB_ACTIONS') else 5.0
            self.assertLess(search_time, expected_time, f"Search should complete under {expected_time} seconds for {seq_size} bp")
            print(f"Large sequence search took {search_time:.3f} seconds for {seq_size} bp")
            
        finally:
            os.unlink(sts_file)
    
    @unittest.skipIf(os.getenv('SKIP_PERFORMANCE_TESTS'), "Performance tests skipped")
    def test_many_sts_loading(self):
        """Test loading many STS records."""
        # Create fewer STS records for CI
        sts_count = 100 if os.getenv('CI') or os.getenv('GITHUB_ACTIONS') else 1000
        sts_content = self.create_test_sts(sts_count)
        
        with tempfile.NamedTemporaryFile(mode='w', suffix='.sts', delete=False) as f:
            f.write(sts_content)
            sts_file = f.name
        
        try:
            start_time = time.time()
            success = self.mer_pcr.load_sts_file(sts_file)
            load_time = time.time() - start_time
            
            self.assertTrue(success)
            expected_records = sts_count * 2  # Each STS creates 2 records (forward and reverse)
            self.assertEqual(len(self.mer_pcr.sts_records), expected_records)
            expected_load_time = 1.0 if os.getenv('CI') or os.getenv('GITHUB_ACTIONS') else 2.0
            self.assertLess(load_time, expected_load_time, f"Loading {sts_count} STSs should be under {expected_load_time} seconds")
            print(f"Loading {sts_count} STSs took {load_time:.3f} seconds")
            
        finally:
            os.unlink(sts_file)
    
    def test_threading_performance(self):
        """Test that threading improves performance on large sequences."""
        # Create a large sequence to ensure threading is used
        seq_size = 1000 if os.getenv('CI') or os.getenv('GITHUB_ACTIONS') else 500000  # 1KB vs 500KB
        large_seq = self.create_large_sequence(seq_size)
        fasta_record = FASTARecord(
            defline=">large_sequence",
            sequence=large_seq
        )
        
        # Create simple test STS data
        sts_content = "TEST001\tATCGATCG\tATCGATCG\t50\tTest STS\n"
        with tempfile.NamedTemporaryFile(mode='w', suffix='.sts', delete=False) as f:
            f.write(sts_content)
            sts_file = f.name
        
        try:
            # Test single-threaded
            mer_pcr_single = MerPCR(threads=1)
            mer_pcr_single.load_sts_file(sts_file)
            
            start_time = time.time()
            hits_single = mer_pcr_single.search([fasta_record])
            time_single = time.time() - start_time
            
            # Test multi-threaded
            mer_pcr_multi = MerPCR(threads=4)
            mer_pcr_multi.load_sts_file(sts_file)
            
            start_time = time.time()
            hits_multi = mer_pcr_multi.search([fasta_record])
            time_multi = time.time() - start_time
            
            # Results should be identical
            self.assertEqual(hits_single, hits_multi)
            
            print(f"Single-threaded: {time_single:.3f}s, Multi-threaded: {time_multi:.3f}s")
            
            # Multi-threaded should be faster (or at least not significantly slower)
            # Allow some variance due to overhead
            self.assertLessEqual(time_multi, time_single * 1.5, 
                               "Multi-threading should not be significantly slower")
            
        finally:
            os.unlink(sts_file)
    
    def test_memory_efficiency(self):
        """Test memory efficiency with large datasets."""
        import psutil
        import os
        
        # Get initial memory usage
        process = psutil.Process(os.getpid())
        initial_memory = process.memory_info().rss / 1024 / 1024  # MB
        
        # Create test data - smaller for CI
        seq_size = 1000 if os.getenv('CI') or os.getenv('GITHUB_ACTIONS') else 1000000  # 1KB vs 1MB
        large_seq = self.create_large_sequence(seq_size)
        fasta_record = FASTARecord(
            defline=">large_sequence",
            sequence=large_seq
        )
        
        # Use simple STS for CI to avoid too many matches
        sts_content = "TEST001\tATCGATCG\tATCGATCG\t50\tTest STS\n" if os.getenv('CI') or os.getenv('GITHUB_ACTIONS') else self.create_test_sts(1000)
        with tempfile.NamedTemporaryFile(mode='w', suffix='.sts', delete=False) as f:
            f.write(sts_content)
            sts_file = f.name
        
        try:
            mer_pcr = MerPCR()
            mer_pcr.load_sts_file(sts_file)
            
            # Redirect output to prevent massive console output in CI  
            import os
            with open(os.devnull, 'w') as devnull:
                mer_pcr.search([fasta_record], output_file=os.devnull)
            
            # Check memory usage after processing
            final_memory = process.memory_info().rss / 1024 / 1024  # MB
            memory_increase = final_memory - initial_memory
            
            sts_count = 1 if os.getenv('CI') or os.getenv('GITHUB_ACTIONS') else 1000
            print(f"Memory usage increased by {memory_increase:.1f} MB for {seq_size} bp sequence with {sts_count} STSs")
            
            # Memory increase should be reasonable 
            expected_memory_limit = 100 if os.getenv('CI') or os.getenv('GITHUB_ACTIONS') else 500
            self.assertLess(memory_increase, expected_memory_limit, 
                          f"Memory usage increase should be under {expected_memory_limit} MB")
            
        finally:
            os.unlink(sts_file)


@pytest.mark.performance
class TestScalability(unittest.TestCase):
    """Test scalability with increasing data sizes."""
    
    def test_sequence_length_scaling(self):
        """Test how performance scales with sequence length."""
        sizes = [10000, 50000, 100000]  # 10KB, 50KB, 100KB
        times = []
        
        # Use non-matching STS for scalability test to avoid excessive output
        sts_content = "TEST001\tGGGGGGGGGGGG\tCCCCCCCCCCCC\t200\tTest STS\n"
        with tempfile.NamedTemporaryFile(mode='w', suffix='.sts', delete=False) as f:
            f.write(sts_content)
            sts_file = f.name
        
        try:
            for size in sizes:
                mer_pcr = MerPCR(threads=1)
                mer_pcr.load_sts_file(sts_file)
                
                # Create sequence of specified size
                sequence = "ATCG" * (size // 4)
                fasta_record = FASTARecord(
                    defline=f">test_seq_{size}",
                    sequence=sequence
                )
                
                start_time = time.time()
                # Redirect output to avoid console spam
                import os
                with open(os.devnull, 'w') as devnull:
                    mer_pcr.search([fasta_record], output_file=os.devnull)
                elapsed = time.time() - start_time
                times.append(elapsed)
                
                print(f"Size: {size:6d} bp, Time: {elapsed:.3f}s")
            
            # Check that time scaling is reasonable (should be roughly linear)
            # Allow for some variance due to system factors
            for i in range(1, len(times)):
                ratio = sizes[i] / sizes[i-1]
                time_ratio = times[i] / times[i-1]
                
                # Time ratio should not be much worse than size ratio
                self.assertLess(time_ratio, ratio * 2, 
                              f"Time scaling should be reasonable: {time_ratio:.2f}x vs {ratio:.2f}x size increase")
        
        finally:
            os.unlink(sts_file)


if __name__ == "__main__":
    # Allow skipping performance tests with environment variable
    if os.getenv('SKIP_PERFORMANCE_TESTS'):
        print("Performance tests skipped due to SKIP_PERFORMANCE_TESTS environment variable")
    
    unittest.main()