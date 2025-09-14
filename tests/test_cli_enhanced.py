"""
Enhanced CLI tests focusing on argument conversion and validation.
"""

import pytest
import sys
from unittest.mock import patch
from merpcr.cli import (
    convert_mepcr_arguments,
    margin_type,
    mismatch_type,
    wordsize_type,
    threads_type,
    pcr_size_type,
    sts_line_length_type,
    create_parser,
    main,
)
import argparse
import tempfile
import os


class TestArgumentConversion:
    """Test me-PCR argument format conversion."""

    def test_basic_mepcr_format(self):
        """Test basic M=50 format conversion."""
        args = ["M=50", "N=1", "W=11"]
        result = convert_mepcr_arguments(args)
        expected = ["-M", "50", "-N", "1", "-W", "11"]
        assert result == expected

    def test_mixed_formats(self):
        """Test mixing me-PCR and regular formats."""
        args = ["M=75", "file1.txt", "-v", "N=2"]
        result = convert_mepcr_arguments(args)
        expected = ["-M", "75", "file1.txt", "-v", "-N", "2"]
        assert result == expected

    def test_all_parameters(self):
        """Test conversion of all supported parameters."""
        args = ["M=100", "N=2", "W=10", "X=2", "T=4", "Q=0", "Z=300", "I=1", "S=2000", "O=out.txt"]
        result = convert_mepcr_arguments(args)
        expected = ["-M", "100", "-N", "2", "-W", "10", "-X", "2", "-T", "4", "-Q", "0", "-Z", "300", "-I", "1", "-S", "2000", "-O", "out.txt"]
        assert result == expected

    def test_help_conversion(self):
        """Test -help to --help conversion."""
        args = ["-help"]
        result = convert_mepcr_arguments(args)
        assert result == ["--help"]

    def test_priority_parameter_ignored(self):
        """Test that P parameter is ignored (Mac-specific)."""
        args = ["M=50", "P=15", "N=1"]
        result = convert_mepcr_arguments(args)
        expected = ["-M", "50", "-N", "1"]
        assert result == expected

    def test_regular_arguments_untouched(self):
        """Test that regular arguments pass through unchanged."""
        args = ["file1.txt", "file2.txt", "-v", "--debug"]
        result = convert_mepcr_arguments(args)
        assert result == args

    def test_empty_args(self):
        """Test empty argument list."""
        assert convert_mepcr_arguments([]) == []

    def test_invalid_mepcr_format(self):
        """Test that invalid formats pass through unchanged."""
        args = ["M50", "=50", "M=", "XX=50"]
        result = convert_mepcr_arguments(args)
        assert result == args


class TestParameterValidation:
    """Test parameter validation functions."""

    def test_margin_validation(self):
        """Test margin parameter validation."""
        assert margin_type("50") == 50
        assert margin_type("0") == 0
        assert margin_type("10000") == 10000
        
        with pytest.raises(argparse.ArgumentTypeError):
            margin_type("-1")
        with pytest.raises(argparse.ArgumentTypeError):
            margin_type("10001")

    def test_mismatch_validation(self):
        """Test mismatch parameter validation."""
        assert mismatch_type("0") == 0
        assert mismatch_type("10") == 10
        
        with pytest.raises(argparse.ArgumentTypeError):
            mismatch_type("-1")
        with pytest.raises(argparse.ArgumentTypeError):
            mismatch_type("11")

    def test_wordsize_validation(self):
        """Test wordsize parameter validation."""
        assert wordsize_type("3") == 3
        assert wordsize_type("16") == 16
        assert wordsize_type("11") == 11
        
        with pytest.raises(argparse.ArgumentTypeError):
            wordsize_type("2")
        with pytest.raises(argparse.ArgumentTypeError):
            wordsize_type("17")

    def test_threads_validation(self):
        """Test threads parameter validation."""
        assert threads_type("1") == 1
        assert threads_type("100") == 100
        
        with pytest.raises(argparse.ArgumentTypeError):
            threads_type("0")
        with pytest.raises(argparse.ArgumentTypeError):
            threads_type("-1")

    def test_pcr_size_validation(self):
        """Test PCR size parameter validation."""
        assert pcr_size_type("1") == 1
        assert pcr_size_type("10000") == 10000
        assert pcr_size_type("240") == 240
        
        with pytest.raises(argparse.ArgumentTypeError):
            pcr_size_type("0")
        with pytest.raises(argparse.ArgumentTypeError):
            pcr_size_type("10001")

    def test_sts_line_length_validation(self):
        """Test STS line length parameter validation."""
        assert sts_line_length_type("1") == 1
        assert sts_line_length_type("1022") == 1022
        assert sts_line_length_type("5000") == 5000
        
        with pytest.raises(argparse.ArgumentTypeError):
            sts_line_length_type("0")


class TestCLIIntegration:
    """Test CLI integration with argument conversion."""

    def test_mepcr_format_integration(self):
        """Test that me-PCR format works end-to-end."""
        # Create temporary test files
        with tempfile.NamedTemporaryFile(mode='w', suffix='.sts', delete=False) as sts_file:
            sts_file.write("TEST\tACGT\tTCGA\t100\ttest marker\n")
            sts_path = sts_file.name

        with tempfile.NamedTemporaryFile(mode='w', suffix='.fa', delete=False) as fasta_file:
            fasta_file.write(">test\nACGTACGTTCGATCGA\n")
            fasta_path = fasta_file.name

        try:
            # Test parsing with me-PCR format
            with patch('sys.argv', ['merpcr', sts_path, fasta_path, 'M=100', 'N=1']):
                parser = create_parser()
                # Simulate argument conversion
                converted_args = convert_mepcr_arguments([sts_path, fasta_path, 'M=100', 'N=1'])
                args = parser.parse_args(converted_args)
                
                assert args.sts_file == sts_path
                assert args.fasta_file == fasta_path
                assert args.margin == 100
                assert args.mismatches == 1

        finally:
            os.unlink(sts_path)
            os.unlink(fasta_path)

    def test_help_format_integration(self):
        """Test that -help works correctly."""
        with patch('sys.argv', ['merpcr', '-help']):
            converted_args = convert_mepcr_arguments(['-help'])
            assert converted_args == ['--help']

    def test_stdout_handling(self):
        """Test stdout output handling."""
        with tempfile.NamedTemporaryFile(mode='w', suffix='.sts', delete=False) as sts_file:
            sts_file.write("TEST\tACGT\tTCGA\t100\n")
            sts_path = sts_file.name

        with tempfile.NamedTemporaryFile(mode='w', suffix='.fa', delete=False) as fasta_file:
            fasta_file.write(">test\nACGTTCGA\n")
            fasta_path = fasta_file.name

        try:
            parser = create_parser()
            # Test O=stdout conversion
            converted_args = convert_mepcr_arguments([sts_path, fasta_path, 'O=stdout'])
            args = parser.parse_args(converted_args)
            assert args.output == "stdout"

        finally:
            os.unlink(sts_path)
            os.unlink(fasta_path)


class TestErrorHandling:
    """Test CLI error handling scenarios."""

    def test_invalid_parameter_values(self):
        """Test that invalid parameter values are caught."""
        parser = create_parser()
        
        # Test invalid margin
        with pytest.raises(SystemExit):
            parser.parse_args(['test.sts', 'test.fa', '-M', '15000'])
        
        # Test invalid wordsize
        with pytest.raises(SystemExit):
            parser.parse_args(['test.sts', 'test.fa', '-W', '20'])

    def test_missing_required_args(self):
        """Test that missing required arguments are caught."""
        parser = create_parser()
        
        with pytest.raises(SystemExit):
            parser.parse_args([])
        
        with pytest.raises(SystemExit):
            parser.parse_args(['test.sts'])  # Missing fasta file

    def test_boundary_values(self):
        """Test boundary values for parameters."""
        parser = create_parser()
        
        # Test boundary values that should work
        args = parser.parse_args(['test.sts', 'test.fa', '-M', '0', '-W', '3', '-N', '10'])
        assert args.margin == 0
        assert args.wordsize == 3
        assert args.mismatches == 10
        
        args = parser.parse_args(['test.sts', 'test.fa', '-M', '10000', '-W', '16'])
        assert args.margin == 10000
        assert args.wordsize == 16