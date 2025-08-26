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

### Removed
- Outdated development files and directories
- Legacy setup.py and requirements.txt
- Duplicate and unused code sections