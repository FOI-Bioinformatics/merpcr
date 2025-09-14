"""
Tests for module entry points and CLI integration.
"""

import logging
import os
import subprocess
import sys
import tempfile
from unittest.mock import Mock, patch

import pytest

from merpcr.cli import main, setup_logging


class TestModuleEntryPoint:
    """Test module entry points."""

    def test_python_module_execution(self):
        """Test that python -m merpcr works correctly."""
        # Create minimal test files
        sts_content = "TEST\tAAAA\tTTTT\t50\n"
        fasta_content = ">test\nAAAATTTT\n"

        with tempfile.NamedTemporaryFile(mode="w", suffix=".sts", delete=False) as sts_f:
            sts_f.write(sts_content)
            sts_path = sts_f.name

        with tempfile.NamedTemporaryFile(mode="w", suffix=".fa", delete=False) as fasta_f:
            fasta_f.write(fasta_content)
            fasta_path = fasta_f.name

        try:
            # Test module execution
            result = subprocess.run(
                [sys.executable, "-m", "merpcr", sts_path, fasta_path, "-W", "4"],
                capture_output=True,
                text=True,
                timeout=10,
            )

            # Should complete without error (exit code 0 expected for successful run)
            assert result.returncode == 0

        finally:
            os.unlink(sts_path)
            os.unlink(fasta_path)

    def test_module_help(self):
        """Test module help output."""
        result = subprocess.run(
            [sys.executable, "-m", "merpcr", "--help"], capture_output=True, text=True, timeout=10
        )

        assert result.returncode == 0
        assert "merPCR - Modern Electronic Rapid PCR" in result.stdout

    def test_module_version(self):
        """Test module version output."""
        result = subprocess.run(
            [sys.executable, "-m", "merpcr", "--version"],
            capture_output=True,
            text=True,
            timeout=10,
        )

        assert result.returncode == 0
        assert "merPCR version" in result.stdout


class TestCLIMainFunction:
    """Test CLI main() function comprehensively."""

    def test_main_successful_execution(self):
        """Test successful main() execution."""
        sts_content = "TEST\tATCG\tCGAT\t50\n"
        fasta_content = ">test\nATCGCGAT\n"

        with tempfile.NamedTemporaryFile(mode="w", suffix=".sts", delete=False) as sts_f:
            sts_f.write(sts_content)
            sts_path = sts_f.name

        with tempfile.NamedTemporaryFile(mode="w", suffix=".fa", delete=False) as fasta_f:
            fasta_f.write(fasta_content)
            fasta_path = fasta_f.name

        try:
            with patch("sys.argv", ["merpcr", sts_path, fasta_path, "-W", "4"]):
                result = main()
                assert result == 0  # Success

        finally:
            os.unlink(sts_path)
            os.unlink(fasta_path)

    def test_main_sts_file_load_failure(self):
        """Test main() when STS file loading fails."""
        fasta_content = ">test\nATCGCGAT\n"

        with tempfile.NamedTemporaryFile(mode="w", suffix=".fa", delete=False) as fasta_f:
            fasta_f.write(fasta_content)
            fasta_path = fasta_f.name

        try:
            # Use nonexistent STS file
            with patch("sys.argv", ["merpcr", "/nonexistent.sts", fasta_path, "-W", "4"]):
                result = main()
                assert result == 1  # Failure

        finally:
            os.unlink(fasta_path)

    def test_main_fasta_file_load_failure(self):
        """Test main() when FASTA file loading fails."""
        sts_content = "TEST\tATCG\tCGAT\t50\n"

        with tempfile.NamedTemporaryFile(mode="w", suffix=".sts", delete=False) as sts_f:
            sts_f.write(sts_content)
            sts_path = sts_f.name

        try:
            # Use nonexistent FASTA file
            with patch("sys.argv", ["merpcr", sts_path, "/nonexistent.fa", "-W", "4"]):
                result = main()
                assert result == 1  # Failure

        finally:
            os.unlink(sts_path)

    def test_main_empty_fasta_file(self):
        """Test main() with empty FASTA file."""
        sts_content = "TEST\tATCG\tCGAT\t50\n"

        with tempfile.NamedTemporaryFile(mode="w", suffix=".sts", delete=False) as sts_f:
            sts_f.write(sts_content)
            sts_path = sts_f.name

        with tempfile.NamedTemporaryFile(mode="w", suffix=".fa", delete=False) as fasta_f:
            # Empty FASTA file
            fasta_path = fasta_f.name

        try:
            with patch("sys.argv", ["merpcr", sts_path, fasta_path, "-W", "4"]):
                result = main()
                assert result == 1  # Should fail due to empty FASTA

        finally:
            os.unlink(sts_path)
            os.unlink(fasta_path)

    def test_main_exception_handling(self):
        """Test main() exception handling."""
        with patch("sys.argv", ["merpcr", "test.sts", "test.fa"]):
            with patch("merpcr.cli.MerPCR") as mock_merpcr:
                # Make MerPCR initialization raise an exception
                mock_merpcr.side_effect = RuntimeError("Test error")

                result = main()
                assert result == 1  # Should handle exception and return 1

    def test_main_exception_with_debug(self):
        """Test main() exception handling with debug mode."""
        with patch("sys.argv", ["merpcr", "test.sts", "test.fa", "--debug"]):
            with patch("merpcr.cli.MerPCR") as mock_merpcr:
                mock_merpcr.side_effect = RuntimeError("Test error")

                with patch("traceback.print_exc") as mock_traceback:
                    result = main()
                    assert result == 1
                    mock_traceback.assert_called_once()  # Debug mode should print traceback

    def test_main_with_mepcr_arguments(self):
        """Test main() with me-PCR style arguments."""
        sts_content = "TEST\tATCG\tCGAT\t50\n"
        fasta_content = ">test\nATCGCGAT\n"

        with tempfile.NamedTemporaryFile(mode="w", suffix=".sts", delete=False) as sts_f:
            sts_f.write(sts_content)
            sts_path = sts_f.name

        with tempfile.NamedTemporaryFile(mode="w", suffix=".fa", delete=False) as fasta_f:
            fasta_f.write(fasta_content)
            fasta_path = fasta_f.name

        try:
            # Use me-PCR style arguments
            with patch("sys.argv", ["merpcr", sts_path, fasta_path, "W=4", "M=100"]):
                result = main()
                assert result == 0  # Should succeed with argument conversion

        finally:
            os.unlink(sts_path)
            os.unlink(fasta_path)

    def test_main_with_output_file(self):
        """Test main() with output file specified."""
        sts_content = "TEST\tATCG\tCGAT\t50\n"
        fasta_content = ">test\nATCGCGAT\n"

        with tempfile.NamedTemporaryFile(mode="w", suffix=".sts", delete=False) as sts_f:
            sts_f.write(sts_content)
            sts_path = sts_f.name

        with tempfile.NamedTemporaryFile(mode="w", suffix=".fa", delete=False) as fasta_f:
            fasta_f.write(fasta_content)
            fasta_path = fasta_f.name

        with tempfile.NamedTemporaryFile(mode="w", suffix=".out", delete=False) as out_f:
            out_path = out_f.name

        try:
            with patch("sys.argv", ["merpcr", sts_path, fasta_path, "-W", "4", "-O", out_path]):
                result = main()
                assert result == 0

                # Check that output file was created
                assert os.path.exists(out_path)

        finally:
            os.unlink(sts_path)
            os.unlink(fasta_path)
            if os.path.exists(out_path):
                os.unlink(out_path)


class TestLoggingSetup:
    """Test logging setup comprehensively."""

    def test_setup_logging_debug_mode(self):
        """Test logging setup with debug mode."""
        # Clear any existing handlers
        logger = logging.getLogger("merpcr")
        for handler in logger.handlers[:]:
            logger.removeHandler(handler)

        setup_logging(quiet=0, debug=True)

        logger = logging.getLogger("merpcr")
        assert logger.level == logging.DEBUG

    def test_setup_logging_verbose_mode(self):
        """Test logging setup with verbose mode (quiet=0)."""
        logger = logging.getLogger("merpcr")
        for handler in logger.handlers[:]:
            logger.removeHandler(handler)

        setup_logging(quiet=0, debug=False)

        logger = logging.getLogger("merpcr")
        assert logger.level == logging.INFO

    def test_setup_logging_quiet_mode(self):
        """Test logging setup with quiet mode (quiet=1)."""
        logger = logging.getLogger("merpcr")
        for handler in logger.handlers[:]:
            logger.removeHandler(handler)

        setup_logging(quiet=1, debug=False)

        logger = logging.getLogger("merpcr")
        assert logger.level == logging.WARNING

    def test_setup_logging_debug_overrides_quiet(self):
        """Test that debug mode overrides quiet setting."""
        logger = logging.getLogger("merpcr")
        for handler in logger.handlers[:]:
            logger.removeHandler(handler)

        setup_logging(quiet=1, debug=True)

        logger = logging.getLogger("merpcr")
        assert logger.level == logging.DEBUG  # Debug should override quiet


class TestCLIErrorPaths:
    """Test CLI error paths and edge cases."""

    def test_invalid_arguments(self):
        """Test CLI with invalid arguments."""
        with patch("sys.argv", ["merpcr", "--invalid-arg"]):
            with pytest.raises(SystemExit):
                main()

    def test_missing_required_arguments(self):
        """Test CLI with missing required arguments."""
        with patch("sys.argv", ["merpcr"]):
            with pytest.raises(SystemExit):
                main()

    def test_invalid_parameter_values(self):
        """Test CLI with invalid parameter values."""
        with patch("sys.argv", ["merpcr", "test.sts", "test.fa", "-W", "100"]):  # Invalid wordsize
            with pytest.raises(SystemExit):
                main()

    def test_help_argument(self):
        """Test CLI help argument."""
        with patch("sys.argv", ["merpcr", "--help"]):
            with pytest.raises(SystemExit) as exc_info:
                main()
            assert exc_info.value.code == 0  # Help should exit with 0

    def test_version_argument(self):
        """Test CLI version argument."""
        with patch("sys.argv", ["merpcr", "--version"]):
            with pytest.raises(SystemExit) as exc_info:
                main()
            assert exc_info.value.code == 0  # Version should exit with 0


class TestCLIIntegrationComplex:
    """Complex integration tests for CLI functionality."""

    def test_complex_sts_file(self):
        """Test with complex STS file containing various edge cases."""
        complex_sts = """# This is a comment
AFM256vb9	TCTGAATGGCCCTTGG	TCCTATCTGAGGTGGGGT	180	(D17S934)  Chr.17, 63.7 cM
AFM248yg9	GCTAAAAATACACGGATGG	TGCAAGACTGCGTCTC	193	(D17S932)  Chr.17, 63.7 cM
SHORT\tAT\tGC\t20\tShort primers (should be skipped)

# Another comment
RANGE_TEST\tATCGATCG\tGCTAGCTA\t100-150\tRange test
"""

        fasta_content = ">test_seq\n" + "ATCG" * 1000 + "\n"  # Long sequence

        with tempfile.NamedTemporaryFile(mode="w", suffix=".sts", delete=False) as sts_f:
            sts_f.write(complex_sts)
            sts_path = sts_f.name

        with tempfile.NamedTemporaryFile(mode="w", suffix=".fa", delete=False) as fasta_f:
            fasta_f.write(fasta_content)
            fasta_path = fasta_f.name

        try:
            with patch("sys.argv", ["merpcr", sts_path, fasta_path, "-W", "8", "-M", "200"]):
                result = main()
                # Should succeed even with complex data
                assert result == 0

        finally:
            os.unlink(sts_path)
            os.unlink(fasta_path)

    def test_all_parameters_integration(self):
        """Test integration with all possible parameters."""
        sts_content = "TEST\tATCGATCGATCG\tGCTAGCTAGCTA\t100\n"
        fasta_content = ">test\nATCGATCGATCGGCTAGCTAGCTA\n"

        with tempfile.NamedTemporaryFile(mode="w", suffix=".sts", delete=False) as sts_f:
            sts_f.write(sts_content)
            sts_path = sts_f.name

        with tempfile.NamedTemporaryFile(mode="w", suffix=".fa", delete=False) as fasta_f:
            fasta_f.write(fasta_content)
            fasta_path = fasta_f.name

        with tempfile.NamedTemporaryFile(mode="w", suffix=".out", delete=False) as out_f:
            out_path = out_f.name

        try:
            # Test with all parameters
            argv = [
                "merpcr",
                sts_path,
                fasta_path,
                "-M",
                "50",
                "-N",
                "1",
                "-W",
                "8",
                "-T",
                "2",
                "-X",
                "1",
                "-Z",
                "240",
                "-I",
                "1",
                "-S",
                "1024",
                "-O",
                out_path,
                "-Q",
                "0",
                "--debug",
            ]

            with patch("sys.argv", argv):
                result = main()
                assert result == 0

        finally:
            os.unlink(sts_path)
            os.unlink(fasta_path)
            if os.path.exists(out_path):
                os.unlink(out_path)
