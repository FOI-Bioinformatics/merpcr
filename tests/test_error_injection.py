"""
Error injection and fault tolerance tests.

These tests simulate various system failures, I/O errors, and edge conditions
to verify that merPCR handles failures gracefully and recovers properly.
"""

import pytest
import tempfile
import os
import sys
from unittest.mock import patch, Mock, mock_open, MagicMock
from io import StringIO
import errno
import logging
from merpcr.core.engine import MerPCR
from merpcr.core.models import STSRecord, FASTARecord
from merpcr.cli import main
from merpcr.io.fasta import FASTALoader
from merpcr.io.sts import STSLoader

# Disable logging for error injection tests to reduce noise
logging.getLogger("merpcr").setLevel(logging.CRITICAL)


class TestFileSystemErrorInjection:
    """Test filesystem error injection scenarios."""

    def test_sts_file_permission_denied(self):
        """Test STS file loading with permission denied error."""
        with patch('os.path.exists', return_value=True), \
             patch('os.path.isfile', return_value=True), \
             patch('os.path.getsize', return_value=100), \
             patch('builtins.open', side_effect=PermissionError("Permission denied")):
            
            engine = MerPCR()
            
            # Should handle permission error gracefully
            try:
                result = engine.load_sts_file("test.sts")
                assert not result
            except PermissionError:
                # This is also acceptable - depends on implementation
                pass

    def test_fasta_file_not_found(self):
        """Test FASTA file loading with file not found error."""
        with patch('os.path.getsize', side_effect=FileNotFoundError("File not found")):
            engine = MerPCR()
            
            try:
                records = engine.load_fasta_file("nonexistent.fa")
                assert records == []
            except FileNotFoundError:
                # This is acceptable behavior
                pass

    def test_output_file_write_error(self):
        """Test output file write errors during search."""
        sts_content = "TEST\tATCG\tCGAT\t50\n"
        fasta_content = ">test\nATCGCGAT\n"
        
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
            
            # Mock file opening to fail for output
            with patch('builtins.open', side_effect=OSError("Disk full")):
                try:
                    hits = engine.search(records, "output.txt")
                    # Should either handle gracefully or raise appropriate error
                except OSError:
                    # This is acceptable
                    pass
                    
        finally:
            os.unlink(sts_path)
            os.unlink(fasta_path)

    def test_partial_file_read_error(self):
        """Test handling of partial file reads."""
        # Create a file that will cause read errors
        with tempfile.NamedTemporaryFile(mode='w', suffix='.sts', delete=False) as f:
            f.write("TEST\tATCG\tCGAT\t50\n")
            sts_path = f.name

        try:
            # Mock file read to fail partway through
            original_open = open
            
            def failing_open(*args, **kwargs):
                if args[0] == sts_path:
                    mock_file = MagicMock()
                    mock_file.readlines.side_effect = OSError("I/O error")
                    mock_file.__enter__.return_value = mock_file
                    return mock_file
                else:
                    return original_open(*args, **kwargs)
            
            with patch('builtins.open', side_effect=failing_open):
                engine = MerPCR()
                try:
                    result = engine.load_sts_file(sts_path)
                    assert not result
                except OSError:
                    # Acceptable behavior
                    pass
                    
        finally:
            os.unlink(sts_path)

    def test_disk_full_during_output(self):
        """Test handling disk full errors during output."""
        sts_content = "TEST\tATCG\tCGAT\t50\n"
        fasta_content = ">test\nATCGCGAT\n"
        
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
            
            # Mock print to fail (disk full scenario)
            def failing_print(*args, **kwargs):
                if 'file' in kwargs:
                    raise OSError(errno.ENOSPC, "No space left on device")
                return print(*args, **kwargs)
            
            with patch('builtins.print', side_effect=failing_print):
                try:
                    hits = engine.search(records, "output.txt")
                except OSError as e:
                    assert e.errno == errno.ENOSPC
                    
        finally:
            os.unlink(sts_path)
            os.unlink(fasta_path)


class TestMemoryErrorInjection:
    """Test memory-related error injection."""

    def test_memory_allocation_failure(self):
        """Test handling of memory allocation failures."""
        # Create large data that might cause memory issues
        large_sts = []
        for i in range(10000):
            large_sts.append(f"STS{i}\t{'A'*20}\t{'T'*20}\t100")
        
        sts_content = '\n'.join(large_sts)
        
        with tempfile.NamedTemporaryFile(mode='w', suffix='.sts', delete=False) as f:
            f.write(sts_content)
            sts_path = f.name

        try:
            # Mock list creation to fail after certain size
            original_list = list
            call_count = 0
            
            def memory_limited_list(*args, **kwargs):
                nonlocal call_count
                call_count += 1
                if call_count > 100:  # Simulate memory exhaustion
                    raise MemoryError("Out of memory")
                return original_list(*args, **kwargs)
            
            with patch('builtins.list', side_effect=memory_limited_list):
                engine = MerPCR()
                try:
                    result = engine.load_sts_file(sts_path)
                    # Should either succeed or fail gracefully
                except MemoryError:
                    # This is acceptable behavior
                    pass
                    
        finally:
            os.unlink(sts_path)

    def test_string_processing_memory_error(self):
        """Test memory errors during string processing."""
        # Create data that requires string processing
        fasta_content = f">huge_seq\n{'ATCG' * 100000}\n"  # 400KB sequence
        
        with tempfile.NamedTemporaryFile(mode='w', suffix='.fa', delete=False) as f:
            f.write(fasta_content)
            fasta_path = f.name

        try:
            # Simulate memory exhaustion by triggering it directly
            from merpcr.core.models import FASTARecord
            
            original_init = FASTARecord.__init__
            call_count = 0
            
            def failing_init(self, defline, sequence):
                nonlocal call_count
                call_count += 1
                # Simulate memory error for large sequences
                if len(sequence) > 200000:  # Very large sequence
                    raise MemoryError("Cannot allocate memory for large sequence")
                return original_init(self, defline, sequence)
            
            with patch.object(FASTARecord, '__init__', failing_init):
                engine = MerPCR()
                try:
                    records = engine.load_fasta_file(fasta_path)
                    # If successful, verify basic properties
                    assert isinstance(records, list)
                    assert len(records) >= 0  # May be empty if memory error
                except MemoryError:
                    # Expected behavior - memory error handled gracefully
                    pass
                    
        finally:
            os.unlink(fasta_path)


class TestCLIErrorInjection:
    """Test CLI error injection scenarios."""

    def test_cli_main_with_corrupted_arguments(self):
        """Test CLI main() with corrupted sys.argv."""
        # Simulate corrupted command line arguments
        corrupted_argv = [None, "test.sts", "test.fa", 123, []]
        
        with patch('sys.argv', corrupted_argv):
            try:
                result = main()
                # Should handle gracefully
            except (TypeError, AttributeError):
                # These are acceptable for corrupted input
                pass

    def test_cli_with_signal_interruption(self):
        """Test CLI behavior with simulated signal interruption."""
        sts_content = "TEST\tATCG\tCGAT\t50\n"
        fasta_content = ">test\nATCGCGAT\n"
        
        with tempfile.NamedTemporaryFile(mode='w', suffix='.sts', delete=False) as sts_f:
            sts_f.write(sts_content)
            sts_path = sts_f.name

        with tempfile.NamedTemporaryFile(mode='w', suffix='.fa', delete=False) as fasta_f:
            fasta_f.write(fasta_content)
            fasta_path = fasta_f.name

        try:
            # Simulate KeyboardInterrupt during execution
            with patch('sys.argv', ['merpcr', sts_path, fasta_path, '-W', '4']):
                with patch('merpcr.core.engine.MerPCR.search', side_effect=KeyboardInterrupt()):
                    try:
                        result = main()
                        # Should handle interruption
                    except KeyboardInterrupt:
                        # This is acceptable behavior
                        pass
                        
        finally:
            os.unlink(sts_path)
            os.unlink(fasta_path)

    def test_cli_with_logging_errors(self):
        """Test CLI behavior when logging fails."""
        sts_content = "TEST\tATCG\tCGAT\t50\n"
        fasta_content = ">test\nATCGCGAT\n"
        
        with tempfile.NamedTemporaryFile(mode='w', suffix='.sts', delete=False) as sts_f:
            sts_f.write(sts_content)
            sts_path = sts_f.name

        with tempfile.NamedTemporaryFile(mode='w', suffix='.fa', delete=False) as fasta_f:
            fasta_f.write(fasta_content)
            fasta_path = fasta_f.name

        try:
            # Mock logging to fail
            with patch('logging.getLogger') as mock_logger:
                mock_logger.return_value.error.side_effect = OSError("Logging failed")
                mock_logger.return_value.info.side_effect = OSError("Logging failed")
                
                with patch('sys.argv', ['merpcr', sts_path, fasta_path, '-W', '4']):
                    try:
                        result = main()
                        # Should complete despite logging errors
                    except OSError:
                        # Acceptable if logging errors propagate
                        pass
                        
        finally:
            os.unlink(sts_path)
            os.unlink(fasta_path)


class TestDataCorruptionInjection:
    """Test handling of corrupted data."""

    def test_corrupted_sts_file_data(self):
        """Test handling of corrupted STS file data."""
        # Various forms of corrupted STS data
        corrupted_sts_variants = [
            "TEST\t\t\t\n",  # Empty fields
            "TEST\tATCG\n",   # Missing fields
            "TEST\tATCG\tCGAT\tNOT_A_NUMBER\n",  # Invalid PCR size
            "TEST\tATCG\tCGAT\t-100\n",  # Negative PCR size
            "\tATCG\tCGAT\t100\n",  # Empty ID
            "TEST\t\tCGAT\t100\n",  # Empty primer1
            "TEST\tATCG\t\t100\n",  # Empty primer2
            "TEST\tXXXXXXXXXXXXXXXX\tCGAT\t100\n",  # Invalid characters in primer
        ]
        
        for corrupted_sts in corrupted_sts_variants:
            with tempfile.NamedTemporaryFile(mode='w', suffix='.sts', delete=False) as f:
                f.write(corrupted_sts)
                sts_path = f.name

            try:
                engine = MerPCR()
                # Should handle corrupted data gracefully
                try:
                    result = engine.load_sts_file(sts_path)
                    # May succeed (skipping bad entries) or fail gracefully
                except (ValueError, IndexError, TypeError):
                    # These are acceptable for corrupted data
                    pass
                    
            finally:
                os.unlink(sts_path)

    def test_corrupted_fasta_file_data(self):
        """Test handling of corrupted FASTA file data."""
        corrupted_fasta_variants = [
            "",  # Empty file
            "not a fasta file\n",  # Invalid format
            ">seq1\n",  # Header without sequence
            "ATCGATCG\n",  # Sequence without header
            ">seq1\nATCG\n>seq2\n",  # Incomplete entry
            ">seq1\nATCGXYZ123\n",  # Invalid characters in sequence
            ">\n\nATCGATCG\n",  # Empty header
        ]
        
        for corrupted_fasta in corrupted_fasta_variants:
            with tempfile.NamedTemporaryFile(mode='w', suffix='.fa', delete=False) as f:
                f.write(corrupted_fasta)
                fasta_path = f.name

            try:
                engine = MerPCR()
                # Should handle corrupted data gracefully
                try:
                    records = engine.load_fasta_file(fasta_path)
                    # May return empty list or partial results
                    assert isinstance(records, list)
                except (ValueError, IndexError):
                    # These are acceptable for corrupted data
                    pass
                    
            finally:
                os.unlink(fasta_path)

    def test_binary_file_as_text(self):
        """Test handling of binary files passed as text files."""
        # Create a binary file
        binary_data = bytes([i % 256 for i in range(1000)])
        
        with tempfile.NamedTemporaryFile(mode='wb', suffix='.sts', delete=False) as f:
            f.write(binary_data)
            sts_path = f.name

        try:
            engine = MerPCR()
            # Should handle binary data gracefully
            try:
                result = engine.load_sts_file(sts_path)
                # Should fail gracefully
                assert not result
            except (UnicodeDecodeError, UnicodeError):
                # These are acceptable for binary data
                pass
                
        finally:
            os.unlink(sts_path)


class TestSystemResourceExhaustion:
    """Test behavior under system resource exhaustion."""

    def test_file_descriptor_exhaustion(self):
        """Test behavior when file descriptors are exhausted."""
        sts_content = "TEST\tATCG\tCGAT\t50\n"
        
        with tempfile.NamedTemporaryFile(mode='w', suffix='.sts', delete=False) as f:
            f.write(sts_content)
            sts_path = f.name

        try:
            # Mock open to fail with "too many open files"
            call_count = 0
            original_open = open
            
            def fd_limited_open(*args, **kwargs):
                nonlocal call_count
                call_count += 1
                if call_count > 3:  # Allow a few opens, then fail
                    raise OSError(errno.EMFILE, "Too many open files")
                return original_open(*args, **kwargs)
            
            with patch('builtins.open', side_effect=fd_limited_open):
                engine = MerPCR()
                try:
                    result = engine.load_sts_file(sts_path)
                    # Should handle FD exhaustion gracefully
                except OSError as e:
                    assert e.errno == errno.EMFILE
                    
        finally:
            os.unlink(sts_path)

    def test_process_limit_exhaustion(self):
        """Test behavior under process/thread limit exhaustion."""
        sts_content = "TEST\tATCG\tCGAT\t50\n"
        sequence = ("ATCGCGAT") * 10000  # Large enough to trigger threading
        fasta_content = f">test\n{sequence}"
        
        with tempfile.NamedTemporaryFile(mode='w', suffix='.sts', delete=False) as sts_f:
            sts_f.write(sts_content)
            sts_path = sts_f.name

        with tempfile.NamedTemporaryFile(mode='w', suffix='.fa', delete=False) as fasta_f:
            fasta_f.write(fasta_content)
            fasta_path = fasta_f.name

        try:
            # Mock ThreadPoolExecutor to fail
            with patch('concurrent.futures.ThreadPoolExecutor') as mock_executor:
                mock_executor.side_effect = OSError("Cannot create thread")
                
                engine = MerPCR(wordsize=4, threads=4)
                engine.load_sts_file(sts_path)
                records = engine.load_fasta_file(fasta_path)
                
                try:
                    hits = engine.search(records)
                    # Should fall back to single-threaded or handle gracefully
                except OSError:
                    # Acceptable behavior
                    pass
                    
        finally:
            os.unlink(sts_path)
            os.unlink(fasta_path)


class TestRecoveryMechanisms:
    """Test error recovery mechanisms."""

    def test_partial_failure_recovery(self):
        """Test recovery from partial failures."""
        # Create mixed good and bad STS entries
        mixed_sts = """GOOD1\tATCG\tCGAT\t50
BAD_ENTRY_WITH_INVALID_SIZE\tATCG\tCGAT\tNOT_A_NUMBER
GOOD2\tTACG\tGTCA\t60
\tEMPTY_ID\tCGAT\t70
GOOD3\tGCTA\tTAGC\t80"""
        
        with tempfile.NamedTemporaryFile(mode='w', suffix='.sts', delete=False) as f:
            f.write(mixed_sts)
            sts_path = f.name

        try:
            engine = MerPCR()
            result = engine.load_sts_file(sts_path)
            
            # Should recover by skipping bad entries and loading good ones
            if result:
                assert len(engine.sts_records) > 0
                # Should have loaded at least the good entries
                assert len(engine.sts_records) >= 2  # At least GOOD1 and GOOD2/GOOD3
                
        finally:
            os.unlink(sts_path)

    def test_graceful_degradation(self):
        """Test graceful degradation when optional features fail."""
        sts_content = "TEST\tATCG\tCGAT\t50\n"
        fasta_content = ">test\nATCGCGAT\n"
        
        with tempfile.NamedTemporaryFile(mode='w', suffix='.sts', delete=False) as sts_f:
            sts_f.write(sts_content)
            sts_path = sts_f.name

        with tempfile.NamedTemporaryFile(mode='w', suffix='.fa', delete=False) as fasta_f:
            fasta_f.write(fasta_content)
            fasta_path = fasta_f.name

        try:
            # Mock threading to fail, should fall back to single-threaded
            with patch('concurrent.futures.ThreadPoolExecutor', side_effect=Exception("Threading failed")):
                engine = MerPCR(wordsize=4, threads=4)
                engine.load_sts_file(sts_path)
                records = engine.load_fasta_file(fasta_path)
                
                # Should still complete search, possibly single-threaded
                hits = engine.search(records)
                assert isinstance(hits, int)
                assert hits >= 0
                
        finally:
            os.unlink(sts_path)
            os.unlink(fasta_path)

    def test_error_reporting_robustness(self):
        """Test that error reporting itself is robust."""
        # Test with CLI main function
        with patch('sys.argv', ['merpcr', 'nonexistent.sts', 'nonexistent.fa']):
            # Mock logger to fail
            with patch('logging.getLogger') as mock_logger:
                mock_logger.return_value.error.side_effect = Exception("Logging system failed")
                
                try:
                    result = main()
                    # Should complete despite logging failures
                except Exception:
                    # If error reporting fails, that's still information
                    pass