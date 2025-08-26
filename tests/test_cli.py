#!/usr/bin/env python3
"""
Tests for command-line interface.
"""

import os
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

import pytest

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "src"))


@pytest.mark.cli
class TestCLI(unittest.TestCase):
    """Test command-line interface functionality."""

    def setUp(self):
        """Set up test environment."""
        self.script_path = Path(__file__).parent.parent / "scripts" / "merpcr"
        self.data_dir = Path(__file__).parent / "data"
        self.sts_file = self.data_dir / "test.sts"
        self.fasta_file = self.data_dir / "test.fa"

        # Check if test files exist
        if not self.sts_file.exists():
            self.skipTest("Test data files not available")

    def test_basic_cli_execution(self):
        """Test basic CLI execution."""
        result = subprocess.run(
            [str(self.script_path), str(self.sts_file), str(self.fasta_file)],
            capture_output=True,
            text=True,
        )

        self.assertEqual(result.returncode, 0)
        self.assertIn("L78833", result.stdout)
        self.assertIn("AFM248yg9", result.stdout)

    def test_verbose_output(self):
        """Test verbose output mode."""
        result = subprocess.run(
            [str(self.script_path), "-Q", "0", str(self.sts_file), str(self.fasta_file)],
            capture_output=True,
            text=True,
        )

        self.assertEqual(result.returncode, 0)
        self.assertIn("Reading STS file", result.stderr)
        self.assertIn("Processing sequence", result.stderr)
        self.assertIn("hits found", result.stderr)

    def test_output_file(self):
        """Test output to file."""
        with tempfile.NamedTemporaryFile(mode="w", delete=False, suffix=".txt") as f:
            output_file = f.name

        try:
            result = subprocess.run(
                [
                    str(self.script_path),
                    "-O",
                    output_file,
                    str(self.sts_file),
                    str(self.fasta_file),
                ],
                capture_output=True,
                text=True,
            )

            self.assertEqual(result.returncode, 0)

            # Check output file was created and contains expected content
            with open(output_file, "r") as f:
                content = f.read()

            self.assertIn("L78833", content)
            self.assertIn("AFM248yg9", content)

        finally:
            if os.path.exists(output_file):
                os.unlink(output_file)

    def test_parameter_parsing(self):
        """Test parameter parsing."""
        result = subprocess.run(
            [
                str(self.script_path),
                "-M",
                "100",  # margin
                "-N",
                "1",  # mismatches
                "-W",
                "10",  # wordsize
                "-T",
                "2",  # threads
                "-X",
                "2",  # three-prime-match
                "-I",
                "1",  # IUPAC mode
                "-Q",
                "1",  # quiet
                str(self.sts_file),
                str(self.fasta_file),
            ],
            capture_output=True,
            text=True,
        )

        self.assertEqual(result.returncode, 0)

    def test_help_message(self):
        """Test help message."""
        result = subprocess.run([str(self.script_path), "--help"], capture_output=True, text=True)

        self.assertEqual(result.returncode, 0)
        self.assertIn("merPCR", result.stdout)
        self.assertIn("usage:", result.stdout.lower())
        self.assertIn("margin", result.stdout.lower())

    def test_version_message(self):
        """Test version message."""
        result = subprocess.run(
            [str(self.script_path), "--version"], capture_output=True, text=True
        )

        self.assertEqual(result.returncode, 0)
        self.assertIn("1.0.0", result.stdout)

    def test_missing_files(self):
        """Test error handling for missing files."""
        result = subprocess.run(
            [str(self.script_path), "/nonexistent/sts.txt", "/nonexistent/fasta.fa"],
            capture_output=True,
            text=True,
        )

        self.assertNotEqual(result.returncode, 0)

    def test_debug_mode(self):
        """Test debug mode."""
        result = subprocess.run(
            [str(self.script_path), "--debug", str(self.sts_file), str(self.fasta_file)],
            capture_output=True,
            text=True,
        )

        self.assertEqual(result.returncode, 0)
        # Debug mode should produce verbose output
        self.assertIn("Reading STS file", result.stderr)


@pytest.mark.cli
class TestModuleExecution(unittest.TestCase):
    """Test execution as Python module."""

    def setUp(self):
        """Set up test environment."""
        self.data_dir = Path(__file__).parent / "data"
        self.sts_file = self.data_dir / "test.sts"
        self.fasta_file = self.data_dir / "test.fa"

        if not self.sts_file.exists():
            self.skipTest("Test data files not available")

    def test_module_execution(self):
        """Test running as module with python -m merpcr."""
        env = os.environ.copy()
        env["PYTHONPATH"] = str(Path(__file__).parent.parent / "src")

        result = subprocess.run(
            [sys.executable, "-m", "merpcr", str(self.sts_file), str(self.fasta_file)],
            capture_output=True,
            text=True,
            env=env,
        )

        self.assertEqual(result.returncode, 0)
        self.assertIn("L78833", result.stdout)
        self.assertIn("AFM248yg9", result.stdout)


if __name__ == "__main__":
    unittest.main()
