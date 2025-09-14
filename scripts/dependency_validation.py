#!/usr/bin/env python3
"""
Comprehensive Cross-Platform Dependency Validation for merPCR

This script validates dependencies across different platforms and Python versions,
ensuring consistent behavior and identifying potential compatibility issues.

Usage:
    python scripts/dependency_validation.py [--check-imports] [--check-versions] [--check-conflicts]
"""

import argparse
import importlib
import platform
import subprocess
import sys
import json
import pkg_resources
from pathlib import Path
from typing import Dict, List, Optional, Tuple, Any
import tempfile
import os

class DependencyValidator:
    """Validates dependencies across platforms and Python versions."""
    
    def __init__(self):
        self.platform_info = {
            "system": platform.system(),
            "machine": platform.machine(),
            "processor": platform.processor(),
            "python_version": f"{sys.version_info.major}.{sys.version_info.minor}",
            "python_implementation": platform.python_implementation(),
            "architecture": platform.architecture()[0],
        }
        
        self.validation_results = {
            "platform": self.platform_info,
            "imports": {},
            "versions": {},
            "conflicts": {},
            "compatibility": {},
            "summary": {
                "total_dependencies": 0,
                "successful_imports": 0,
                "failed_imports": 0,
                "version_mismatches": 0,
                "conflicts_found": 0,
                "platform_issues": 0
            }
        }
    
    def get_installed_packages(self) -> Dict[str, str]:
        """Get all installed packages and their versions."""
        try:
            installed_packages = {}
            for dist in pkg_resources.working_set:
                installed_packages[dist.project_name.lower()] = dist.version
            return installed_packages
        except Exception as e:
            print(f"Error getting installed packages: {e}")
            return {}
    
    def check_core_dependencies(self) -> bool:
        """Check that core merPCR dependencies are available and importable."""
        print("Validating core dependencies...")
        
        # Core Python modules that merPCR depends on
        core_modules = [
            "os", "sys", "pathlib", "tempfile", "threading", "multiprocessing",
            "json", "time", "argparse", "logging", "collections", "itertools",
            "functools", "operator", "warnings", "traceback", "concurrent.futures"
        ]
        
        # Optional dependencies
        optional_modules = [
            "pytest", "psutil", "coverage", "hypothesis"
        ]
        
        success_count = 0
        total_count = len(core_modules) + len(optional_modules)
        
        # Test core modules
        for module in core_modules:
            try:
                importlib.import_module(module)
                self.validation_results["imports"][module] = {"status": "success", "type": "core"}
                success_count += 1
                print(f"  ✓ {module}")
            except ImportError as e:
                self.validation_results["imports"][module] = {
                    "status": "failed", 
                    "type": "core", 
                    "error": str(e)
                }
                print(f"  ✗ {module}: {e}")
        
        # Test optional modules
        for module in optional_modules:
            try:
                importlib.import_module(module)
                self.validation_results["imports"][module] = {"status": "success", "type": "optional"}
                success_count += 1
                print(f"  ✓ {module} (optional)")
            except ImportError as e:
                self.validation_results["imports"][module] = {
                    "status": "failed",
                    "type": "optional", 
                    "error": str(e)
                }
                print(f"  ? {module} (optional): {e}")
        
        self.validation_results["summary"]["successful_imports"] = success_count
        self.validation_results["summary"]["failed_imports"] = total_count - success_count
        
        # Core modules must all be available
        core_failures = [m for m in core_modules if 
                        self.validation_results["imports"][m]["status"] == "failed"]
        
        if core_failures:
            print(f"  CRITICAL: Core modules failed to import: {core_failures}")
            return False
        else:
            print(f"  SUCCESS: All core dependencies validated ({success_count}/{total_count} total)")
            return True
    
    def check_merpcr_import(self) -> bool:
        """Test importing merPCR and basic functionality."""
        print("Validating merPCR import and basic functionality...")
        
        try:
            # Test basic import
            import merpcr
            print(f"  ✓ merpcr import successful")
            print(f"  ✓ merpcr version: {merpcr.__version__}")
            
            # Test core class import
            from merpcr import MerPCR
            print(f"  ✓ MerPCR class import successful")
            
            # Test basic instantiation
            engine = MerPCR()
            print(f"  ✓ MerPCR instantiation successful")
            
            # Test method availability
            required_methods = ["load_sts_file", "load_fasta_file", "search"]
            for method in required_methods:
                if hasattr(engine, method):
                    print(f"  ✓ Method {method} available")
                else:
                    print(f"  ✗ Method {method} missing")
                    return False
            
            # Test CLI module
            try:
                import merpcr.__main__
                print(f"  ✓ CLI module available")
            except ImportError:
                print(f"  ? CLI module not available (may be expected)")
            
            self.validation_results["imports"]["merpcr"] = {
                "status": "success", 
                "version": merpcr.__version__,
                "methods": required_methods
            }
            
            return True
            
        except ImportError as e:
            print(f"  ✗ merPCR import failed: {e}")
            self.validation_results["imports"]["merpcr"] = {
                "status": "failed", 
                "error": str(e)
            }
            return False
        except Exception as e:
            print(f"  ✗ merPCR functionality test failed: {e}")
            self.validation_results["imports"]["merpcr"] = {
                "status": "partial", 
                "error": str(e)
            }
            return False
    
    def check_version_consistency(self) -> bool:
        """Check for version consistency across dependencies."""
        print("Checking dependency version consistency...")
        
        installed_packages = self.get_installed_packages()
        self.validation_results["summary"]["total_dependencies"] = len(installed_packages)
        
        # Define expected version ranges for key dependencies
        version_requirements = {
            "pytest": {"min": "6.0", "max": "8.0"},
            "coverage": {"min": "5.0", "max": "8.0"},
            "psutil": {"min": "5.0", "max": "6.0"},
        }
        
        issues_found = 0
        
        for package, requirements in version_requirements.items():
            if package in installed_packages:
                installed_version = installed_packages[package]
                
                try:
                    from packaging import version
                    installed_ver = version.parse(installed_version)
                    min_ver = version.parse(requirements["min"])
                    max_ver = version.parse(requirements["max"])
                    
                    if installed_ver < min_ver or installed_ver >= max_ver:
                        print(f"  ⚠️ {package}: version {installed_version} outside recommended range [{requirements['min']}, {requirements['max']})")
                        issues_found += 1
                        self.validation_results["versions"][package] = {
                            "installed": installed_version,
                            "required_min": requirements["min"],
                            "required_max": requirements["max"],
                            "status": "out_of_range"
                        }
                    else:
                        print(f"  ✓ {package}: version {installed_version} OK")
                        self.validation_results["versions"][package] = {
                            "installed": installed_version,
                            "status": "ok"
                        }
                except ImportError:
                    print(f"  ? {package}: cannot verify version (packaging module not available)")
                    self.validation_results["versions"][package] = {
                        "installed": installed_version,
                        "status": "cannot_verify"
                    }
            else:
                if package == "pytest":  # Required for testing
                    print(f"  ⚠️ {package}: not installed (required for testing)")
                    issues_found += 1
                else:
                    print(f"  ? {package}: not installed (optional)")
        
        self.validation_results["summary"]["version_mismatches"] = issues_found
        return issues_found == 0
    
    def check_dependency_conflicts(self) -> bool:
        """Check for dependency conflicts."""
        print("Checking for dependency conflicts...")
        
        try:
            # Use pip check to find dependency conflicts
            result = subprocess.run([sys.executable, "-m", "pip", "check"], 
                                  capture_output=True, text=True, timeout=30)
            
            if result.returncode == 0:
                print("  ✓ No dependency conflicts detected")
                self.validation_results["conflicts"]["pip_check"] = {
                    "status": "no_conflicts",
                    "output": result.stdout
                }
                return True
            else:
                print("  ⚠️ Dependency conflicts detected:")
                conflicts = result.stdout.strip().split('\n') if result.stdout else []
                for conflict in conflicts:
                    if conflict.strip():
                        print(f"    {conflict}")
                
                self.validation_results["conflicts"]["pip_check"] = {
                    "status": "conflicts_found",
                    "conflicts": conflicts,
                    "output": result.stdout
                }
                self.validation_results["summary"]["conflicts_found"] = len(conflicts)
                return False
                
        except subprocess.TimeoutExpired:
            print("  ⚠️ pip check timed out")
            self.validation_results["conflicts"]["pip_check"] = {
                "status": "timeout"
            }
            return False
        except Exception as e:
            print(f"  ⚠️ Error checking conflicts: {e}")
            self.validation_results["conflicts"]["pip_check"] = {
                "status": "error",
                "error": str(e)
            }
            return False
    
    def check_platform_compatibility(self) -> bool:
        """Check platform-specific compatibility issues."""
        print(f"Checking platform compatibility for {self.platform_info['system']}...")
        
        platform_issues = 0
        
        # Check Python version compatibility
        if sys.version_info < (3, 11):
            print(f"  ⚠️ Python {sys.version_info.major}.{sys.version_info.minor} is not supported (requires 3.11+)")
            platform_issues += 1
        else:
            print(f"  ✓ Python {sys.version_info.major}.{sys.version_info.minor} version OK")
        
        # Check platform-specific modules
        platform_modules = {
            "Windows": ["winsound"],
            "Darwin": [],  # macOS
            "Linux": []
        }
        
        current_platform = self.platform_info['system']
        if current_platform in platform_modules:
            for module in platform_modules[current_platform]:
                try:
                    importlib.import_module(module)
                    print(f"  ✓ Platform module {module} available")
                except ImportError:
                    print(f"  ? Platform module {module} not available (may be expected)")
        
        # Test file system operations
        try:
            with tempfile.NamedTemporaryFile(mode='w', suffix='.test', delete=False) as tmp:
                tmp.write("test")
                tmp_path = tmp.name
            
            # Test file operations
            test_path = Path(tmp_path)
            if test_path.exists():
                test_path.unlink()
                print("  ✓ File system operations working")
            else:
                print("  ⚠️ File system operations may have issues")
                platform_issues += 1
                
        except Exception as e:
            print(f"  ⚠️ File system test failed: {e}")
            platform_issues += 1
        
        # Test multiprocessing
        try:
            import multiprocessing
            cpu_count = multiprocessing.cpu_count()
            print(f"  ✓ Multiprocessing available ({cpu_count} CPUs)")
            
            # Test basic multiprocessing
            def test_worker():
                return os.getpid()
            
            with multiprocessing.Pool(1) as pool:
                result = pool.apply(test_worker)
                if result:
                    print("  ✓ Multiprocessing test successful")
                else:
                    print("  ⚠️ Multiprocessing test failed")
                    platform_issues += 1
                    
        except Exception as e:
            print(f"  ⚠️ Multiprocessing test failed: {e}")
            platform_issues += 1
        
        self.validation_results["summary"]["platform_issues"] = platform_issues
        self.validation_results["compatibility"] = {
            "python_version_ok": sys.version_info >= (3, 11),
            "filesystem_ok": platform_issues == 0,
            "multiprocessing_ok": platform_issues == 0,
            "platform_issues": platform_issues
        }
        
        return platform_issues == 0
    
    def run_comprehensive_validation(self) -> bool:
        """Run all validation checks."""
        print("=== merPCR Dependency Validation ===")
        print(f"Platform: {self.platform_info['system']} {self.platform_info['machine']}")
        print(f"Python: {self.platform_info['python_version']} ({self.platform_info['python_implementation']})")
        print("")
        
        all_passed = True
        
        # Run all validation checks
        checks = [
            ("Core Dependencies", self.check_core_dependencies),
            ("merPCR Import", self.check_merpcr_import),
            ("Version Consistency", self.check_version_consistency),
            ("Dependency Conflicts", self.check_dependency_conflicts),
            ("Platform Compatibility", self.check_platform_compatibility)
        ]
        
        results = {}
        for check_name, check_func in checks:
            print(f"\n--- {check_name} ---")
            try:
                passed = check_func()
                results[check_name] = passed
                if not passed:
                    all_passed = False
            except Exception as e:
                print(f"  ERROR: {check_name} check failed: {e}")
                results[check_name] = False
                all_passed = False
        
        # Summary
        print("\n=== Validation Summary ===")
        for check_name, passed in results.items():
            status = "✅ PASS" if passed else "❌ FAIL"
            print(f"{status} {check_name}")
        
        if all_passed:
            print("\n🎉 All dependency validations passed!")
        else:
            print(f"\n⚠️ Some validations failed. Review the details above.")
        
        return all_passed
    
    def save_results(self, output_file: str = "dependency_validation_results.json"):
        """Save validation results to a JSON file."""
        with open(output_file, 'w') as f:
            json.dump(self.validation_results, f, indent=2)
        print(f"\nValidation results saved to {output_file}")


def main():
    parser = argparse.ArgumentParser(description="merPCR Dependency Validation")
    parser.add_argument("--check-imports", action="store_true",
                       help="Check import availability")
    parser.add_argument("--check-versions", action="store_true", 
                       help="Check version consistency")
    parser.add_argument("--check-conflicts", action="store_true",
                       help="Check for dependency conflicts")
    parser.add_argument("--output", type=str, default="dependency_validation_results.json",
                       help="Output file for results")
    parser.add_argument("--comprehensive", action="store_true", default=True,
                       help="Run all validation checks (default)")
    
    args = parser.parse_args()
    
    validator = DependencyValidator()
    
    if args.comprehensive or not any([args.check_imports, args.check_versions, args.check_conflicts]):
        # Run comprehensive validation
        success = validator.run_comprehensive_validation()
    else:
        # Run specific checks
        success = True
        if args.check_imports:
            success &= validator.check_core_dependencies()
            success &= validator.check_merpcr_import()
        if args.check_versions:
            success &= validator.check_version_consistency()
        if args.check_conflicts:
            success &= validator.check_dependency_conflicts()
    
    # Save results
    validator.save_results(args.output)
    
    # Exit with appropriate code
    sys.exit(0 if success else 1)


if __name__ == "__main__":
    main()