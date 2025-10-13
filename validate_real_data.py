#!/usr/bin/env python3
"""
Comprehensive validation script for merPCR vs me-PCR on real genomic data.

This script:
1. Runs both me-PCR and merPCR on all genomes with multiple parameter configs
2. Compares outputs line-by-line
3. Generates detailed validation report
4. Does NOT commit any real data to the repository

Usage:
    python validate_real_data.py [--quick]  # --quick runs only default config
"""

import argparse
import json
import os
import subprocess
import sys
import tempfile
import time
from collections import defaultdict
from dataclasses import dataclass, field
from pathlib import Path
from typing import Dict, List, Optional, Tuple

# Configuration
MEPCR_BIN = "/Users/andreassjodin/Code/me-PCR/src/me-PCR"
MERPCR_CMD = ["python", "-m", "merpcr"]
DATA_DIR = "/Users/andreassjodin/Downloads/merpcr_data"
PRIMERS = f"{DATA_DIR}/primers.sts"
GENOMES_DIR = f"{DATA_DIR}/genomes"

# Test configurations (name, parameters)
TEST_CONFIGS = [
    ("default", []),
    ("high_margin", ["-M", "100"]),
    ("mismatches_1", ["-N", "1"]),
    ("mismatches_2", ["-N", "2"]),
    ("wordsize_9", ["-W", "9"]),
    ("wordsize_13", ["-W", "13"]),
    ("iupac_mode", ["-I", "1"]),
    ("stress_test", ["-M", "75", "-N", "1", "-W", "10"]),
]


@dataclass
class Hit:
    """Represents a single PCR hit."""

    seq_label: str
    pos1: int
    pos2: int
    sts_id: str
    alias: str
    direction: str

    @classmethod
    def from_line(cls, line: str) -> Optional["Hit"]:
        """Parse a hit from output line."""
        if not line.strip():
            return None

        try:
            parts = line.strip().split("\t")
            if len(parts) < 5:
                return None

            seq_label = parts[0]
            positions = parts[1].split("..")
            pos1 = int(positions[0])
            pos2 = int(positions[1])
            sts_id = parts[2]
            alias = parts[3] if len(parts) > 3 else ""
            direction = parts[4].strip("()") if len(parts) > 4 else ""

            return cls(seq_label, pos1, pos2, sts_id, alias, direction)
        except (ValueError, IndexError) as e:
            print(f"Warning: Failed to parse line: {line.strip()} - {e}")
            return None

    def to_line(self) -> str:
        """Convert hit back to output format."""
        return f"{self.seq_label}\t{self.pos1}..{self.pos2}\t{self.sts_id}\t{self.alias}\t({self.direction})"

    def __hash__(self):
        return hash((self.seq_label, self.pos1, self.pos2, self.sts_id, self.direction))

    def __eq__(self, other):
        if not isinstance(other, Hit):
            return False
        return (
            self.seq_label == other.seq_label
            and self.pos1 == other.pos1
            and self.pos2 == other.pos2
            and self.sts_id == other.sts_id
            and self.direction == other.direction
        )


@dataclass
class ComparisonResult:
    """Results of comparing two output files."""

    genome: str
    config: str
    identical: bool
    mepcr_hits: int
    merpcr_hits: int
    mepcr_time: float
    merpcr_time: float
    missing_in_merpcr: List[Hit] = field(default_factory=list)
    extra_in_merpcr: List[Hit] = field(default_factory=list)
    position_diffs: List[Tuple[Hit, Hit]] = field(default_factory=list)


class Validator:
    """Main validation orchestrator."""

    def __init__(self, output_dir: str, quick: bool = False):
        self.output_dir = Path(output_dir)
        self.quick = quick
        self.results = []

        # Create output directory structure
        self.output_dir.mkdir(parents=True, exist_ok=True)

        # Verify tools exist
        if not Path(MEPCR_BIN).exists():
            raise FileNotFoundError(f"me-PCR not found at: {MEPCR_BIN}")
        if not Path(PRIMERS).exists():
            raise FileNotFoundError(f"Primers file not found at: {PRIMERS}")
        if not Path(GENOMES_DIR).exists():
            raise FileNotFoundError(f"Genomes directory not found at: {GENOMES_DIR}")

    def get_genomes(self) -> List[Path]:
        """Get list of genome files to test."""
        genomes = sorted(Path(GENOMES_DIR).glob("*.fasta"))
        print(f"Found {len(genomes)} genome files")
        return genomes

    def run_mepcr(self, genome: Path, params: List[str], output_file: Path) -> float:
        """Run me-PCR and return execution time."""
        # me-PCR writes to file specified by O= parameter, not stdout
        cmd = [MEPCR_BIN, PRIMERS, str(genome)] + params + [f"O={output_file}", "Q=1"]

        start = time.time()
        try:
            result = subprocess.run(
                cmd, capture_output=True, text=True, timeout=300, check=False
            )
            elapsed = time.time() - start

            # me-PCR writes directly to the file, no need to capture stdout
            # Just verify the file exists and has content
            if not output_file.exists():
                print(f"Warning: me-PCR didn't create output file {output_file}")
                # Create empty file so comparison works
                output_file.touch()

            if result.returncode != 0:
                print(f"Warning: me-PCR returned {result.returncode}")
                if result.stderr:
                    print(f"stderr: {result.stderr}")

            return elapsed
        except subprocess.TimeoutExpired:
            print(f"ERROR: me-PCR timed out for {genome.name}")
            return -1
        except Exception as e:
            print(f"ERROR: me-PCR failed for {genome.name}: {e}")
            return -1

    def run_merpcr(self, genome: Path, params: List[str], output_file: Path) -> float:
        """Run merPCR and return execution time."""
        cmd = MERPCR_CMD + [PRIMERS, str(genome)] + params + ["-O", str(output_file), "-Q", "1"]

        start = time.time()
        try:
            result = subprocess.run(
                cmd, capture_output=True, text=True, timeout=300, check=False
            )
            elapsed = time.time() - start

            if result.returncode != 0:
                print(f"Warning: merPCR returned {result.returncode}")
                print(f"stderr: {result.stderr}")

            return elapsed
        except subprocess.TimeoutExpired:
            print(f"ERROR: merPCR timed out for {genome.name}")
            return -1
        except Exception as e:
            print(f"ERROR: merPCR failed for {genome.name}: {e}")
            return -1

    def load_hits(self, output_file: Path) -> List[Hit]:
        """Load hits from output file."""
        if not output_file.exists():
            return []

        hits = []
        with open(output_file) as f:
            for line in f:
                hit = Hit.from_line(line)
                if hit:
                    hits.append(hit)
        return hits

    def compare_hits(
        self, mepcr_hits: List[Hit], merpcr_hits: List[Hit]
    ) -> Tuple[List[Hit], List[Hit], List[Tuple[Hit, Hit]]]:
        """Compare two lists of hits."""
        mepcr_set = set(mepcr_hits)
        merpcr_set = set(merpcr_hits)

        missing = list(mepcr_set - merpcr_set)
        extra = list(merpcr_set - mepcr_set)

        # Look for position differences (same STS, different position)
        position_diffs = []
        for mep_hit in mepcr_hits:
            for mer_hit in merpcr_hits:
                if (
                    mep_hit.seq_label == mer_hit.seq_label
                    and mep_hit.sts_id == mer_hit.sts_id
                    and mep_hit.direction == mer_hit.direction
                    and (mep_hit.pos1 != mer_hit.pos1 or mep_hit.pos2 != mer_hit.pos2)
                ):
                    position_diffs.append((mep_hit, mer_hit))

        return missing, extra, position_diffs

    def validate_genome(
        self, genome: Path, config_name: str, params: List[str]
    ) -> ComparisonResult:
        """Validate one genome with one parameter configuration."""
        print(f"  Testing {genome.name} with config '{config_name}'...", end=" ", flush=True)

        # Create output directories
        config_dir = self.output_dir / config_name
        config_dir.mkdir(exist_ok=True)
        mepcr_dir = config_dir / "mepcr"
        merpcr_dir = config_dir / "merpcr"
        mepcr_dir.mkdir(exist_ok=True)
        merpcr_dir.mkdir(exist_ok=True)

        # Output files
        mepcr_output = mepcr_dir / f"{genome.stem}.txt"
        merpcr_output = merpcr_dir / f"{genome.stem}.txt"

        # Run both tools
        mepcr_time = self.run_mepcr(genome, params, mepcr_output)
        merpcr_time = self.run_merpcr(genome, params, merpcr_output)

        # Load hits
        mepcr_hits = self.load_hits(mepcr_output)
        merpcr_hits = self.load_hits(merpcr_output)

        # Compare
        missing, extra, position_diffs = self.compare_hits(mepcr_hits, merpcr_hits)

        identical = len(missing) == 0 and len(extra) == 0 and len(position_diffs) == 0

        result = ComparisonResult(
            genome=genome.name,
            config=config_name,
            identical=identical,
            mepcr_hits=len(mepcr_hits),
            merpcr_hits=len(merpcr_hits),
            mepcr_time=mepcr_time,
            merpcr_time=merpcr_time,
            missing_in_merpcr=missing,
            extra_in_merpcr=extra,
            position_diffs=position_diffs,
        )

        status = "✅ PASS" if identical else "❌ FAIL"
        print(
            f"{status} (me-PCR: {len(mepcr_hits)} hits, merPCR: {len(merpcr_hits)} hits)"
        )

        return result

    def run_validation(self):
        """Run full validation suite."""
        genomes = self.get_genomes()
        configs = [TEST_CONFIGS[0]] if self.quick else TEST_CONFIGS

        print(f"\n{'='*70}")
        print(f"Real Data Validation: merPCR vs me-PCR")
        print(f"{'='*70}")
        print(f"Genomes: {len(genomes)}")
        print(f"Configs: {len(configs)}")
        print(f"Total tests: {len(genomes) * len(configs)}")
        print(f"{'='*70}\n")

        for config_name, params in configs:
            print(f"\nConfiguration: {config_name} {' '.join(params)}")
            print(f"{'-'*70}")

            for genome in genomes:
                result = self.validate_genome(genome, config_name, params)
                self.results.append(result)

    def generate_report(self):
        """Generate comprehensive validation report."""
        report_file = self.output_dir / "REAL_DATA_VALIDATION.md"

        # Calculate statistics
        total_tests = len(self.results)
        passed = sum(1 for r in self.results if r.identical)
        failed = total_tests - passed

        # Group by config
        by_config = defaultdict(list)
        for r in self.results:
            by_config[r.config].append(r)

        with open(report_file, "w") as f:
            f.write("# Real Genomic Data Validation Report\n\n")
            f.write(f"**Date:** {time.strftime('%Y-%m-%d %H:%M:%S')}\n")
            f.write(f"**Status:** {'✅ PASSED' if failed == 0 else '❌ FAILED'}\n\n")

            f.write("## Overall Summary\n\n")
            f.write(f"- **Total Tests:** {total_tests}\n")
            f.write(f"- **Passed:** {passed}\n")
            f.write(f"- **Failed:** {failed}\n")
            f.write(f"- **Success Rate:** {100 * passed / total_tests:.1f}%\n\n")

            f.write("## Results by Configuration\n\n")
            f.write("| Configuration | Genomes | Passed | Failed | Match Rate |\n")
            f.write("|--------------|---------|--------|--------|------------|\n")

            for config_name in [c[0] for c in TEST_CONFIGS]:
                if config_name not in by_config:
                    continue
                results = by_config[config_name]
                config_passed = sum(1 for r in results if r.identical)
                config_total = len(results)
                match_rate = 100 * config_passed / config_total if config_total > 0 else 0

                f.write(
                    f"| {config_name} | {config_total} | {config_passed} | "
                    f"{config_total - config_passed} | {match_rate:.1f}% |\n"
                )

            # Performance summary
            f.write("\n## Performance Summary\n\n")
            mepcr_times = [r.mepcr_time for r in self.results if r.mepcr_time > 0]
            merpcr_times = [r.merpcr_time for r in self.results if r.merpcr_time > 0]

            if mepcr_times and merpcr_times:
                avg_mepcr = sum(mepcr_times) / len(mepcr_times)
                avg_merpcr = sum(merpcr_times) / len(merpcr_times)
                ratio = avg_merpcr / avg_mepcr

                f.write("| Tool | Avg Time | Total Time | Ratio |\n")
                f.write("|------|----------|------------|-------|\n")
                f.write(f"| me-PCR | {avg_mepcr:.2f}s | {sum(mepcr_times):.1f}s | 1.0x |\n")
                f.write(
                    f"| merPCR | {avg_merpcr:.2f}s | {sum(merpcr_times):.1f}s | {ratio:.2f}x |\n"
                )

            # Detailed failures (if any)
            failures = [r for r in self.results if not r.identical]
            if failures:
                f.write("\n## ⚠️ Differences Found\n\n")
                for result in failures:
                    f.write(f"### {result.genome} - {result.config}\n\n")

                    if result.missing_in_merpcr:
                        f.write(f"**Missing in merPCR:** {len(result.missing_in_merpcr)}\n")
                        for hit in result.missing_in_merpcr[:5]:  # Show first 5
                            f.write(f"- `{hit.to_line()}`\n")

                    if result.extra_in_merpcr:
                        f.write(f"**Extra in merPCR:** {len(result.extra_in_merpcr)}\n")
                        for hit in result.extra_in_merpcr[:5]:  # Show first 5
                            f.write(f"- `{hit.to_line()}`\n")

                    if result.position_diffs:
                        f.write(f"**Position differences:** {len(result.position_diffs)}\n")
                        for mepcr_hit, merpcr_hit in result.position_diffs[:5]:
                            f.write(f"- me-PCR: `{mepcr_hit.to_line()}`\n")
                            f.write(f"- merPCR: `{merpcr_hit.to_line()}`\n")

                    f.write("\n")
            else:
                f.write("\n## ✅ Perfect Match\n\n")
                f.write("All outputs match identically between me-PCR and merPCR!\n\n")

            f.write("\n## Conclusion\n\n")
            if failed == 0:
                f.write(
                    "✅ **merPCR produces identical results to me-PCR on real genomic data.**\n\n"
                )
                f.write("The implementation is verified and ready for production use.\n")
            else:
                f.write(f"❌ **Found {failed} differences that need investigation.**\n\n")
                f.write("See detailed differences above.\n")

        print(f"\n{'='*70}")
        print(f"Validation Complete!")
        print(f"{'='*70}")
        print(f"Report generated: {report_file}")
        print(f"Results: {passed}/{total_tests} passed ({100*passed/total_tests:.1f}%)")
        print(f"{'='*70}\n")

        return report_file


def main():
    parser = argparse.ArgumentParser(description="Validate merPCR against me-PCR on real data")
    parser.add_argument(
        "--quick",
        action="store_true",
        help="Quick mode: test only default parameters",
    )
    parser.add_argument(
        "--output-dir",
        default="/tmp/merpcr_real_validation",
        help="Output directory for results",
    )

    args = parser.parse_args()

    validator = Validator(args.output_dir, quick=args.quick)

    try:
        validator.run_validation()
        report_file = validator.generate_report()

        # Return success if all tests passed
        failed = sum(1 for r in validator.results if not r.identical)
        sys.exit(0 if failed == 0 else 1)

    except KeyboardInterrupt:
        print("\n\nValidation interrupted by user")
        sys.exit(1)
    except Exception as e:
        print(f"\n\nERROR: {e}")
        import traceback

        traceback.print_exc()
        sys.exit(1)


if __name__ == "__main__":
    main()
