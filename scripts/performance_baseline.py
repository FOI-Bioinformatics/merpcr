#!/usr/bin/env python3
"""
Performance Baseline Management System for merPCR

This script manages performance baselines for CI/CD pipelines, allowing for
systematic performance regression detection across different platforms and
Python versions.

Usage:
    python scripts/performance_baseline.py establish [--platform] [--python-version]
    python scripts/performance_baseline.py compare [--platform] [--python-version] [--threshold]
    python scripts/performance_baseline.py report [--format json|text|html]
"""

import json
import time
import platform
import sys
import os
import argparse
import statistics
from pathlib import Path
from typing import Dict, List, Optional, Tuple, Any
import subprocess
import tempfile

class PerformanceBaseline:
    """Manages performance baselines for merPCR operations."""
    
    def __init__(self, baseline_dir: str = ".benchmarks"):
        self.baseline_dir = Path(baseline_dir)
        self.baseline_dir.mkdir(exist_ok=True)
        
        # System information
        self.platform_info = {
            "system": platform.system(),
            "machine": platform.machine(),
            "processor": platform.processor(),
            "python_version": f"{sys.version_info.major}.{sys.version_info.minor}",
            "python_implementation": platform.python_implementation(),
        }
        
        self.baseline_file = self.baseline_dir / f"baseline_{self.get_platform_key()}.json"
    
    def get_platform_key(self) -> str:
        """Generate a unique key for the current platform."""
        return f"{self.platform_info['system']}-{self.platform_info['machine']}-py{self.platform_info['python_version']}"
    
    def create_performance_datasets(self) -> Dict[str, Path]:
        """Create datasets for performance testing."""
        test_data_dir = Path("tests/data/performance")
        test_data_dir.mkdir(parents=True, exist_ok=True)
        
        datasets = {}
        
        # Small dataset for quick tests
        small_sts = test_data_dir / "small_performance.sts"
        if not small_sts.exists():
            with open(small_sts, 'w') as f:
                for i in range(100):
                    f.write(f"SMALL_STS_{i:03d}\t{'ATCG' * 5}\t{'GCTA' * 5}\t{200 + i % 50}\tSmall test {i}\n")
        datasets["small_sts"] = small_sts
        
        # Medium dataset
        medium_sts = test_data_dir / "medium_performance.sts"
        if not medium_sts.exists():
            with open(medium_sts, 'w') as f:
                for i in range(1000):
                    primer1 = 'ATCG' * (4 + i % 3)  # Variable length primers
                    primer2 = 'GCTA' * (4 + i % 3)
                    f.write(f"MED_STS_{i:04d}\t{primer1}\t{primer2}\t{200 + i % 100}\tMedium test {i}\n")
        datasets["medium_sts"] = medium_sts
        
        # Large dataset for stress testing
        large_sts = test_data_dir / "large_performance.sts"
        if not large_sts.exists():
            with open(large_sts, 'w') as f:
                for i in range(5000):
                    primer1 = 'ATCGATCGATCG' + 'ATCG' * (i % 5)
                    primer2 = 'GCTAGCTAGCTA' + 'GCTA' * (i % 5)
                    f.write(f"LARGE_STS_{i:05d}\t{primer1}\t{primer2}\t{150 + i % 200}\tLarge test {i}\n")
        datasets["large_sts"] = large_sts
        
        # Genomic sequence for search testing
        sequence_fa = test_data_dir / "performance_sequence.fa"
        if not sequence_fa.exists():
            with open(sequence_fa, 'w') as f:
                f.write(">performance_test_chromosome Synthetic sequence for performance testing\n")
                # 500KB sequence with embedded STS sites
                sequence = ""
                for i in range(6250):  # 6250 * 80 = 500KB
                    if i % 100 == 0:  # Embed some STS sites
                        sequence += "ATCGATCGATCGATCG" + "N" * 48 + "GCTAGCTAGCTAGCTA"  # 80 bp
                    else:
                        sequence += "ATCG" * 20  # 80 bp of random-ish sequence
                    f.write(sequence[-80:] + "\n")
                    sequence = sequence[-80:]  # Keep only last chunk for memory efficiency
        datasets["sequence_fa"] = sequence_fa
        
        return datasets
    
    def measure_operation(self, operation_name: str, operation_func, *args, **kwargs) -> Dict[str, Any]:
        """Measure the performance of a single operation."""
        measurements = []
        memory_measurements = []
        
        # Run operation multiple times for statistical significance
        for run in range(5):
            try:
                import psutil
                process = psutil.Process()
                
                # Memory before
                memory_before = process.memory_info().rss / 1024 / 1024  # MB
                
                # Time the operation
                start_time = time.perf_counter()
                result = operation_func(*args, **kwargs)
                end_time = time.perf_counter()
                
                # Memory after
                memory_after = process.memory_info().rss / 1024 / 1024  # MB
                memory_used = memory_after - memory_before
                
                execution_time = end_time - start_time
                measurements.append(execution_time)
                memory_measurements.append(memory_used)
                
                print(f"  Run {run + 1}: {execution_time:.4f}s, Memory: {memory_used:.2f}MB")
                
            except Exception as e:
                print(f"  Run {run + 1} failed: {e}")
                continue
        
        if not measurements:
            return None
        
        return {
            "operation": operation_name,
            "mean_time": statistics.mean(measurements),
            "median_time": statistics.median(measurements),
            "stdev_time": statistics.stdev(measurements) if len(measurements) > 1 else 0,
            "min_time": min(measurements),
            "max_time": max(measurements),
            "mean_memory": statistics.mean(memory_measurements),
            "max_memory": max(memory_measurements),
            "measurements": measurements,
            "memory_measurements": memory_measurements,
            "platform": self.platform_info.copy(),
            "timestamp": time.time(),
        }
    
    def establish_baseline(self) -> Dict[str, Any]:
        """Establish performance baseline for current platform."""
        print(f"Establishing performance baseline for {self.get_platform_key()}")
        
        # Create test datasets
        datasets = self.create_performance_datasets()
        
        # Import merPCR
        try:
            from merpcr import MerPCR
        except ImportError:
            print("Error: merPCR not installed. Run 'pip install -e .' first.")
            return None
        
        baseline_results = {}
        
        # Test 1: STS file loading performance
        print("Measuring STS file loading performance...")
        
        def load_sts_small():
            engine = MerPCR()
            return engine.load_sts_file(str(datasets["small_sts"]))
        
        def load_sts_medium():
            engine = MerPCR()
            return engine.load_sts_file(str(datasets["medium_sts"]))
        
        def load_sts_large():
            engine = MerPCR()
            return engine.load_sts_file(str(datasets["large_sts"]))
        
        baseline_results["load_sts_small"] = self.measure_operation("load_sts_small", load_sts_small)
        baseline_results["load_sts_medium"] = self.measure_operation("load_sts_medium", load_sts_medium)
        baseline_results["load_sts_large"] = self.measure_operation("load_sts_large", load_sts_large)
        
        # Test 2: FASTA file loading performance
        print("Measuring FASTA file loading performance...")
        
        def load_fasta():
            engine = MerPCR()
            return engine.load_fasta_file(str(datasets["sequence_fa"]))
        
        baseline_results["load_fasta"] = self.measure_operation("load_fasta", load_fasta)
        
        # Test 3: Search performance
        print("Measuring search performance...")
        
        def search_small():
            engine = MerPCR(wordsize=8, margin=50)
            engine.load_sts_file(str(datasets["small_sts"]))
            records = engine.load_fasta_file(str(datasets["sequence_fa"]))
            return engine.search(records[:1])  # Search first sequence only
        
        def search_medium():
            engine = MerPCR(wordsize=8, margin=50)
            engine.load_sts_file(str(datasets["medium_sts"]))
            records = engine.load_fasta_file(str(datasets["sequence_fa"]))
            return engine.search(records[:1])
        
        baseline_results["search_small"] = self.measure_operation("search_small", search_small)
        baseline_results["search_medium"] = self.measure_operation("search_medium", search_medium)
        
        # Test 4: Multi-threaded performance
        print("Measuring multi-threaded performance...")
        
        def search_multithreaded():
            engine = MerPCR(wordsize=8, margin=50, threads=2)
            engine.load_sts_file(str(datasets["medium_sts"]))
            records = engine.load_fasta_file(str(datasets["sequence_fa"]))
            return engine.search(records[:1])
        
        baseline_results["search_multithreaded"] = self.measure_operation("search_multithreaded", search_multithreaded)
        
        # Save baseline
        with open(self.baseline_file, 'w') as f:
            json.dump(baseline_results, f, indent=2)
        
        print(f"Baseline established and saved to {self.baseline_file}")
        return baseline_results
    
    def compare_with_baseline(self, threshold_percent: float = 20.0) -> Dict[str, Any]:
        """Compare current performance with established baseline."""
        if not self.baseline_file.exists():
            print(f"No baseline found for {self.get_platform_key()}. Establishing new baseline.")
            return self.establish_baseline()
        
        print(f"Comparing performance against baseline for {self.get_platform_key()}")
        
        # Load existing baseline
        with open(self.baseline_file, 'r') as f:
            baseline = json.load(f)
        
        # Measure current performance
        current_results = self.establish_baseline()
        
        if not current_results:
            return None
        
        # Compare results
        comparison_results = {
            "baseline_timestamp": baseline.get("search_small", {}).get("timestamp"),
            "current_timestamp": current_results["search_small"]["timestamp"],
            "platform": self.platform_info,
            "comparisons": {},
            "summary": {
                "regressions": 0,
                "improvements": 0,
                "stable": 0,
                "failed_tests": 0
            }
        }
        
        for operation in current_results:
            if operation not in baseline:
                continue
            
            baseline_time = baseline[operation]["mean_time"]
            current_time = current_results[operation]["mean_time"]
            
            change_percent = ((current_time - baseline_time) / baseline_time) * 100
            
            status = "stable"
            if change_percent > threshold_percent:
                status = "regression"
                comparison_results["summary"]["regressions"] += 1
            elif change_percent < -threshold_percent:
                status = "improvement"
                comparison_results["summary"]["improvements"] += 1
            else:
                comparison_results["summary"]["stable"] += 1
            
            comparison_results["comparisons"][operation] = {
                "baseline_time": baseline_time,
                "current_time": current_time,
                "change_percent": change_percent,
                "change_absolute": current_time - baseline_time,
                "status": status,
                "baseline_memory": baseline[operation]["mean_memory"],
                "current_memory": current_results[operation]["mean_memory"],
                "memory_change_percent": ((current_results[operation]["mean_memory"] - 
                                         baseline[operation]["mean_memory"]) / 
                                        baseline[operation]["mean_memory"]) * 100 if baseline[operation]["mean_memory"] > 0 else 0
            }
        
        # Save comparison results
        comparison_file = self.baseline_dir / f"comparison_{self.get_platform_key()}_{int(time.time())}.json"
        with open(comparison_file, 'w') as f:
            json.dump(comparison_results, f, indent=2)
        
        return comparison_results
    
    def generate_report(self, format_type: str = "text", output_file: Optional[str] = None):
        """Generate a performance report."""
        if not self.baseline_file.exists():
            print("No baseline data available. Run 'establish' first.")
            return
        
        # Find latest comparison
        comparison_files = list(self.baseline_dir.glob(f"comparison_{self.get_platform_key()}_*.json"))
        if not comparison_files:
            print("No comparison data available. Run 'compare' first.")
            return
        
        latest_comparison_file = max(comparison_files, key=lambda x: x.stat().st_mtime)
        
        with open(latest_comparison_file, 'r') as f:
            comparison_data = json.load(f)
        
        if format_type == "text":
            report = self._generate_text_report(comparison_data)
        elif format_type == "json":
            report = json.dumps(comparison_data, indent=2)
        elif format_type == "html":
            report = self._generate_html_report(comparison_data)
        else:
            print(f"Unknown format: {format_type}")
            return
        
        if output_file:
            with open(output_file, 'w') as f:
                f.write(report)
            print(f"Report saved to {output_file}")
        else:
            print(report)
    
    def _generate_text_report(self, comparison_data: Dict[str, Any]) -> str:
        """Generate a text format performance report."""
        lines = []
        lines.append("=== merPCR Performance Analysis Report ===")
        lines.append(f"Platform: {comparison_data['platform']['system']} {comparison_data['platform']['machine']}")
        lines.append(f"Python: {comparison_data['platform']['python_version']} ({comparison_data['platform']['python_implementation']})")
        lines.append(f"Analysis Date: {time.ctime(comparison_data['current_timestamp'])}")
        lines.append("")
        
        # Summary
        summary = comparison_data['summary']
        lines.append("Performance Summary:")
        lines.append(f"  Regressions: {summary['regressions']}")
        lines.append(f"  Improvements: {summary['improvements']}")
        lines.append(f"  Stable: {summary['stable']}")
        lines.append(f"  Failed Tests: {summary['failed_tests']}")
        lines.append("")
        
        # Detailed results
        lines.append("Detailed Analysis:")
        for operation, data in comparison_data['comparisons'].items():
            status_symbol = {
                'regression': '🔴',
                'improvement': '🟢',
                'stable': '⚪'
            }.get(data['status'], '❓')
            
            lines.append(f"{status_symbol} {operation}")
            lines.append(f"    Time: {data['baseline_time']:.4f}s → {data['current_time']:.4f}s ({data['change_percent']:+.1f}%)")
            lines.append(f"    Memory: {data['baseline_memory']:.2f}MB → {data['current_memory']:.2f}MB ({data['memory_change_percent']:+.1f}%)")
            lines.append("")
        
        return "\n".join(lines)
    
    def _generate_html_report(self, comparison_data: Dict[str, Any]) -> str:
        """Generate an HTML format performance report."""
        html = f"""
<!DOCTYPE html>
<html>
<head>
    <title>merPCR Performance Report</title>
    <style>
        body {{ font-family: Arial, sans-serif; margin: 20px; }}
        .regression {{ color: red; }}
        .improvement {{ color: green; }}
        .stable {{ color: gray; }}
        table {{ border-collapse: collapse; width: 100%; }}
        th, td {{ border: 1px solid #ddd; padding: 8px; text-align: left; }}
        th {{ background-color: #f2f2f2; }}
    </style>
</head>
<body>
    <h1>merPCR Performance Analysis Report</h1>
    <p><strong>Platform:</strong> {comparison_data['platform']['system']} {comparison_data['platform']['machine']}</p>
    <p><strong>Python:</strong> {comparison_data['platform']['python_version']} ({comparison_data['platform']['python_implementation']})</p>
    <p><strong>Analysis Date:</strong> {time.ctime(comparison_data['current_timestamp'])}</p>
    
    <h2>Performance Summary</h2>
    <ul>
        <li>Regressions: {comparison_data['summary']['regressions']}</li>
        <li>Improvements: {comparison_data['summary']['improvements']}</li>
        <li>Stable: {comparison_data['summary']['stable']}</li>
        <li>Failed Tests: {comparison_data['summary']['failed_tests']}</li>
    </ul>
    
    <h2>Detailed Results</h2>
    <table>
        <tr>
            <th>Operation</th>
            <th>Baseline Time (s)</th>
            <th>Current Time (s)</th>
            <th>Change (%)</th>
            <th>Memory Change (%)</th>
            <th>Status</th>
        </tr>
        """
        
        for operation, data in comparison_data['comparisons'].items():
            status_class = data['status']
            html += f"""
        <tr class="{status_class}">
            <td>{operation}</td>
            <td>{data['baseline_time']:.4f}</td>
            <td>{data['current_time']:.4f}</td>
            <td>{data['change_percent']:+.1f}%</td>
            <td>{data['memory_change_percent']:+.1f}%</td>
            <td>{data['status'].title()}</td>
        </tr>
            """
        
        html += """
    </table>
</body>
</html>
        """
        
        return html


def main():
    parser = argparse.ArgumentParser(description="merPCR Performance Baseline Management")
    parser.add_argument("action", choices=["establish", "compare", "report"], 
                       help="Action to perform")
    parser.add_argument("--threshold", type=float, default=20.0,
                       help="Performance regression threshold percentage (default: 20.0)")
    parser.add_argument("--format", choices=["text", "json", "html"], default="text",
                       help="Report format (default: text)")
    parser.add_argument("--output", type=str,
                       help="Output file for report")
    
    args = parser.parse_args()
    
    baseline_manager = PerformanceBaseline()
    
    if args.action == "establish":
        baseline_manager.establish_baseline()
    elif args.action == "compare":
        results = baseline_manager.compare_with_baseline(args.threshold)
        if results:
            print(f"\nPerformance comparison completed:")
            print(f"  Regressions: {results['summary']['regressions']}")
            print(f"  Improvements: {results['summary']['improvements']}")
            print(f"  Stable: {results['summary']['stable']}")
            
            if results['summary']['regressions'] > 0:
                print("\n⚠️ Performance regressions detected!")
                sys.exit(1)
            else:
                print("\n✅ No significant performance regressions detected")
    elif args.action == "report":
        baseline_manager.generate_report(args.format, args.output)


if __name__ == "__main__":
    main()