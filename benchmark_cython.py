#!/usr/bin/env python3
"""
Benchmark script to compare Cython vs Pure Python performance.
"""

import time
import subprocess
import sys
from pathlib import Path

# Paths
MERPCR_DIR = Path("/Users/andreassjodin/Code/merpcr")
DATA_DIR = Path("/Users/andreassjodin/Downloads/merpcr_data")
PRIMERS_FILE = DATA_DIR / "primers.sts"
GENOMES_DIR = DATA_DIR / "genomes"
CYTHON_SO = MERPCR_DIR / "src/merpcr/core/engine_fast.cpython-313-darwin.so"

# Test genomes (varying sizes)
TEST_GENOMES = [
    "Francisellaceae_Francisella_tularensis_GCF_000833475.1.fasta",  # 1.9 MB
    "Francisellaceae_Cysteiniphilum_litorale_GCF_003428165.1.fasta",  # 3.2 MB
    "Francisellaceae_Fangia_hongkongensis_GCF_000379445.1.fasta",  # 2.9 MB
]


def run_merpcr(genome_file, output_file):
    """Run merPCR and measure time."""
    start = time.time()
    result = subprocess.run(
        ["merpcr", str(PRIMERS_FILE), str(genome_file), "-O", str(output_file)],
        capture_output=True,
        text=True
    )
    elapsed = time.time() - start

    if result.returncode != 0:
        print(f"Error running merPCR: {result.stderr}")
        return None, 0

    # Count hits
    try:
        with open(output_file) as f:
            hits = len(f.readlines())
    except:
        hits = 0

    return elapsed, hits


def main():
    print("=" * 70)
    print("merPCR Performance Benchmark: Cython vs Pure Python")
    print("=" * 70)
    print()

    results = []

    for genome_name in TEST_GENOMES:
        genome_file = GENOMES_DIR / genome_name
        if not genome_file.exists():
            print(f"⚠️  Skipping {genome_name} (not found)")
            continue

        genome_size = genome_file.stat().st_size / (1024 * 1024)  # MB
        print(f"Testing: {genome_name} ({genome_size:.1f} MB)")
        print("-" * 70)

        # Test with Cython
        if CYTHON_SO.exists():
            output_file = Path("/tmp/cython_test.out")
            time_cython, hits_cython = run_merpcr(genome_file, output_file)
            if time_cython is not None:
                print(f"  Cython:      {time_cython:.4f}s  ({hits_cython} hits)")
        else:
            time_cython = None
            hits_cython = 0
            print("  Cython:      Not available")

        # Test with Pure Python (rename SO temporarily)
        backup_so = CYTHON_SO.with_suffix(".so.bak")
        if CYTHON_SO.exists():
            CYTHON_SO.rename(backup_so)

        output_file = Path("/tmp/python_test.out")
        time_python, hits_python = run_merpcr(genome_file, output_file)
        if time_python is not None:
            print(f"  Pure Python: {time_python:.4f}s  ({hits_python} hits)")

        # Restore Cython
        if backup_so.exists():
            backup_so.rename(CYTHON_SO)

        # Verify identical results
        if hits_cython == hits_python:
            print(f"  ✅ Results match: {hits_cython} hits")
        else:
            print(f"  ⚠️  Hit count mismatch: Cython={hits_cython}, Python={hits_python}")

        # Calculate speedup
        if time_cython and time_python:
            speedup = time_python / time_cython
            print(f"  Speedup:     {speedup:.2f}x")
            results.append({
                'genome': genome_name,
                'size_mb': genome_size,
                'time_cython': time_cython,
                'time_python': time_python,
                'speedup': speedup,
                'hits': hits_cython
            })

        print()

    # Summary
    if results:
        print("=" * 70)
        print("Summary")
        print("=" * 70)
        total_time_cython = sum(r['time_cython'] for r in results)
        total_time_python = sum(r['time_python'] for r in results)
        avg_speedup = sum(r['speedup'] for r in results) / len(results)

        print(f"Total genomes tested:  {len(results)}")
        print(f"Total time (Cython):   {total_time_cython:.4f}s")
        print(f"Total time (Python):   {total_time_python:.4f}s")
        print(f"Overall speedup:       {total_time_python / total_time_cython:.2f}x")
        print(f"Average speedup:       {avg_speedup:.2f}x")
        print()

        # Performance table
        print("Detailed Results:")
        print(f"{'Genome':<50} {'Size':>8} {'Cython':>10} {'Python':>10} {'Speedup':>8}")
        print("-" * 92)
        for r in results:
            genome_short = r['genome'][:47] + "..." if len(r['genome']) > 50 else r['genome']
            print(f"{genome_short:<50} {r['size_mb']:>6.1f}MB "
                  f"{r['time_cython']:>9.4f}s {r['time_python']:>9.4f}s "
                  f"{r['speedup']:>7.2f}x")


if __name__ == "__main__":
    main()
