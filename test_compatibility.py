#!/usr/bin/env python3
"""
Comprehensive test suite to verify merPCR compatibility with me-PCR.

This script runs both me-PCR and merPCR with various parameter combinations
and verifies that the outputs match exactly.
"""

import os
import subprocess
import sys
import tempfile
from pathlib import Path

# Paths to executables and test data
MEPCR_PATH = "/Users/andreassjodin/Code/me-PCR/src/me-PCR"
MERPCR_CMD = ["python", "-m", "merpcr"]
STS_FILE = "/Users/andreassjodin/Code/me-PCR/test/test.sts"
FASTA_FILE = "/Users/andreassjodin/Code/me-PCR/test/test.fa"

def run_mepcr(args):
    """Run me-PCR with given arguments."""
    cmd = [MEPCR_PATH, STS_FILE, FASTA_FILE] + args
    result = subprocess.run(cmd, capture_output=True, text=True, cwd="/Users/andreassjodin/Code/me-PCR")
    return result.stdout.strip()

def run_merpcr(args):
    """Run merPCR with given arguments."""
    cmd = MERPCR_CMD + [STS_FILE, FASTA_FILE] + args
    result = subprocess.run(cmd, capture_output=True, text=True)
    return result.stdout.strip()

def compare_outputs(test_name, mepcr_args, merpcr_args):
    """Compare outputs from me-PCR and merPCR."""
    print(f"Testing: {test_name}")
    
    try:
        mepcr_output = run_mepcr(mepcr_args)
        merpcr_output = run_merpcr(merpcr_args)
        
        # Extract just the result lines (ignore verbose output)
        mepcr_lines = [line for line in mepcr_output.split('\n') if '\t' in line and '..' in line]
        merpcr_lines = [line for line in merpcr_output.split('\n') if '\t' in line and '..' in line]
        
        if mepcr_lines == merpcr_lines:
            print(f"✓ {test_name}: PASS")
            return True
        else:
            print(f"✗ {test_name}: FAIL")
            print(f"  me-PCR:  {mepcr_lines}")
            print(f"  merPCR:  {merpcr_lines}")
            return False
    except Exception as e:
        print(f"✗ {test_name}: ERROR - {e}")
        return False

def test_argument_formats():
    """Test different argument formats."""
    tests = [
        ("Default parameters", [], []),
        ("Custom margin (me-PCR format)", ["M=100"], ["M=100"]),
        ("Custom margin (argparse format)", [], ["-M", "100"]),
        ("Multiple parameters (me-PCR format)", ["M=75", "N=1"], ["M=75", "N=1"]),
        ("Multiple parameters (argparse format)", [], ["-M", "75", "-N", "1"]),
        ("Word size test", ["W=10"], ["W=10"]),
        ("IUPAC mode", ["I=1"], ["I=1"]),
        ("Custom PCR size", ["Z=300"], ["Z=300"]),
        ("STS line length", ["S=2000"], ["S=2000"]),
    ]
    
    passed = 0
    total = len(tests)
    
    for test_name, mepcr_args, merpcr_args in tests:
        if compare_outputs(test_name, mepcr_args, merpcr_args):
            passed += 1
    
    return passed, total

def test_output_formats():
    """Test different output formats."""
    with tempfile.NamedTemporaryFile(mode='w', suffix='.txt', delete=False) as tmp:
        tmp_path = tmp.name
    
    try:
        # Test file output
        mepcr_out = run_mepcr([f"O={tmp_path}"])
        with open(tmp_path, 'r') as f:
            mepcr_file_content = f.read().strip()
        
        merpcr_out = run_merpcr(["-O", tmp_path])
        with open(tmp_path, 'r') as f:
            merpcr_file_content = f.read().strip()
        
        if mepcr_file_content == merpcr_file_content:
            print("✓ File output: PASS")
            file_passed = True
        else:
            print("✗ File output: FAIL")
            file_passed = False
        
        # Test stdout output
        mepcr_stdout = run_mepcr(["O=stdout"])
        merpcr_stdout = run_merpcr(["-O", "stdout"])
        
        if mepcr_stdout == merpcr_stdout:
            print("✓ Stdout output: PASS")
            stdout_passed = True
        else:
            print("✗ Stdout output: FAIL")
            stdout_passed = False
        
        return (1 if file_passed else 0) + (1 if stdout_passed else 0), 2
        
    finally:
        if os.path.exists(tmp_path):
            os.unlink(tmp_path)

def test_error_handling():
    """Test error handling for invalid parameters."""
    # These should fail gracefully
    error_tests = [
        ("Invalid wordsize", ["-W", "20"]),
        ("Invalid margin", ["-M", "15000"]),
        ("Invalid mismatches", ["-N", "15"]),
        ("Zero threads", ["-T", "0"]),
    ]
    
    passed = 0
    total = len(error_tests)
    
    for test_name, args in error_tests:
        result = subprocess.run(MERPCR_CMD + [STS_FILE, FASTA_FILE] + args, 
                              capture_output=True, text=True)
        if result.returncode != 0:
            print(f"✓ {test_name}: PASS (correctly failed)")
            passed += 1
        else:
            print(f"✗ {test_name}: FAIL (should have failed)")
    
    return passed, total

def main():
    """Run all compatibility tests."""
    print("merPCR Compatibility Test Suite")
    print("=" * 40)
    
    if not os.path.exists(MEPCR_PATH):
        print(f"Error: me-PCR not found at {MEPCR_PATH}")
        sys.exit(1)
    
    if not os.path.exists(STS_FILE) or not os.path.exists(FASTA_FILE):
        print("Error: Test data files not found")
        sys.exit(1)
    
    total_passed = 0
    total_tests = 0
    
    # Test argument formats
    print("\n1. Testing argument formats:")
    passed, total = test_argument_formats()
    total_passed += passed
    total_tests += total
    
    # Test output formats
    print("\n2. Testing output formats:")
    passed, total = test_output_formats()
    total_passed += passed
    total_tests += total
    
    # Test error handling
    print("\n3. Testing error handling:")
    passed, total = test_error_handling()
    total_passed += passed
    total_tests += total
    
    # Summary
    print("\n" + "=" * 40)
    print(f"Results: {total_passed}/{total_tests} tests passed")
    
    if total_passed == total_tests:
        print("🎉 All tests passed! merPCR is compatible with me-PCR.")
        sys.exit(0)
    else:
        print("❌ Some tests failed. Please review the output above.")
        sys.exit(1)

if __name__ == "__main__":
    main()