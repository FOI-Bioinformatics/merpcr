# Changelog

All notable changes to this project will be documented in this file.

## [1.0.0] - 2024-08-26

### Added
- Complete refactoring of merPCR into a proper Python package structure
- Modular architecture with separate core, IO, and CLI modules  
- Comprehensive test suite with real me-PCR test data
- Modern Python packaging with pyproject.toml
- Type hints throughout the codebase
- Proper logging configuration
- **GitHub Actions CI/CD Pipeline**:
  - Matrix testing across Python 3.9-3.12 on Ubuntu, macOS, and Windows
  - Automated code quality checks with flake8, black, and isort
  - Security scanning with bandit, safety, pip-audit, and CodeQL
  - Test coverage reporting with codecov integration
  - Automated performance benchmarking
  - Automatic code formatting on pull requests
- **Enhanced Testing Infrastructure**:
  - 76 tests across 5 categories: unit, integration, CLI, performance, and I/O
  - pytest configuration with custom test markers
  - Performance benchmarks with CI-optimized settings
  - Memory efficiency and scalability testing
  - Cross-platform compatibility validation
- **Professional Documentation**:
  - Comprehensive API reference documentation
  - Detailed user guide with parameter tuning
  - Step-by-step tutorials and real-world examples
  - Performance optimization guide
- **Development Tools**:
  - tox for multi-environment testing
  - Coverage reporting with .coveragerc
  - Makefile for convenient development commands
  - GitHub issue and PR templates

### Changed
- Restructured project into `src/merpcr/` package layout
- Split monolithic `merPCR.py` into logical modules:
  - `core/engine.py` - Main search engine
  - `core/models.py` - Data models
  - `io/fasta.py` - FASTA file handling  
  - `cli.py` - Command-line interface
- Moved tests to `tests/` directory with proper data organization
- Updated import statements for new package structure

### Fixed  
- Duplicate STS processing bug causing incorrect record counts
- Boolean return type issues in sequence comparison
- Reverse strand search algorithm accuracy
- All tests now pass with restructured codebase
- **Test Suite Corrections**:
  - Fixed ambiguous base reverse complement calculations in test expectations
  - Corrected FASTA loader test assertions to match actual implementation behavior
  - Fixed file handling tests to properly expect FileNotFoundError for missing files
  - Updated sequence filtering tests to account for all IUPAC nucleotide codes
  - Optimized performance tests for CI environments to prevent timeouts

### Removed
- Outdated development files and directories
- Legacy setup.py and requirements.txt
- Duplicate and unused code sections