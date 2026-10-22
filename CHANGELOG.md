# Change Log

This document records significant modifications and enhancements to the merPCR software package following semantic versioning principles.

## [1.0.2] - Structured results from the Python API

- Repository layout: ad-hoc comparison and benchmark scripts moved from the
  root to scripts/; the redundant scripts/merpcr launcher removed (use the
  installed merpcr command or python -m merpcr); pytest and coverage
  configuration moved into pyproject.toml; make clean removes all local
  artifacts.
- New `MerPCR.find_hits(records)` returns the hits as `SearchHit` objects
  (sequence label, 1-based start and end, STS record, strand, product size,
  and `format()` for the output line). `search()` now calls it and writes the
  lines; the last results are kept in `MerPCR.hits`.

## [1.0.1] - Audit against me-PCR 1.0.6

Every behavior below was checked by running the me-PCR 1.0.6 binary on the same
input. A differential test harness (scripts/compare_with_mepcr.py) and the
regression tests in tests/test_mepcr_fidelity.py document the corrections.

### Corrections for me-PCR compatibility
- Output: the alias column is omitted when the STS line has no fifth field, and
  text after the fourth tab is reported verbatim (including further tabs).
- FASTA parsing: uses the me-PCR nucleotide alphabet. X and U are dropped from
  sequences; the gap character "-" is kept and counts as a position.
- IUPAC mode: matching is one-directional, as in me-PCR. A primer code matches
  the sequence bases it expands to, but an ambiguity code in the sequence does
  not match a plain primer base. Only primers containing ambiguity codes are
  compared this way.
- Hashing: U is treated as ambiguous when choosing the hash word.
- Hit order: hits are reported in scan order, chunk by chunk, with records that
  share a hash value checked most-recent-first, instead of being sorted by
  start position.
- Threading: duplicate hits in chunk overlaps were not removed because
  absolute positions were compared with a chunk-relative limit. Chunk sizes now
  use the same rounding as me-PCR and the threading decision uses the FASTA
  file size.
- STS parsing: lines longer than the -S limit are an error; an STS whose first
  primer has no hash is skipped entirely; the secondary loader in
  merpcr.io.sts and merpcr.core.utils.hash_value now share the engine logic.
- A defline consisting only of ">" gives an empty label instead of an error.
- The legacy "-margin N" option of me-PCR is accepted.

### Tests and documentation
- New regression tests (tests/test_mepcr_fidelity.py) and a differential test
  harness against the me-PCR binary (scripts/compare_with_mepcr.py).
- The threading performance test no longer asserts on wall-clock time, and the
  memory test uses a random sequence instead of a periodic one that matched
  every marker at every position.
- Documentation consolidated into README.md and five files under docs/.
  CI reduced to a test/lint/build workflow and a release workflow.

### Requirements
- Python 3.13 or later is required. Formatting and import order are now
  enforced in CI with black and isort.

### Deliberate differences from me-PCR 1.0.6
- PCR size ranges ("200-220") and "-" are accepted as described in the me-PCR
  source; the 1.0.6 binary rejects them with a format error.
- STS lines following a comment line are reported correctly; me-PCR 1.0.6
  reads the wrong line when reporting such hits.

## [1.0.0] - Production Release

### Implementation Enhancements
- Architectural restructuring into modular Python package design
- Implementation of separation of concerns through distinct core, I/O, and interface modules
- Development of comprehensive validation framework utilizing authentic me-PCR reference datasets
- Adoption of contemporary Python packaging standards via pyproject.toml
- Integration of comprehensive type annotation system
- Implementation of structured logging infrastructure with configurable verbosity levels
- **Continuous Integration Framework**:
  - Multi-platform validation across Python 3.8-3.12 on Linux, macOS, and Windows environments
  - Automated code quality assessment using industry-standard linting tools
  - Comprehensive security vulnerability scanning including static analysis and dependency auditing
  - Integrated test coverage reporting with threshold enforcement
  - Automated performance regression testing
  - Standardized code formatting through automated tooling
- **Comprehensive Testing Infrastructure**:
  - 221 validation tests spanning unit, integration, CLI, performance, and I/O categories
  - Advanced testing methodologies including property-based testing, stress testing, and error injection
  - pytest framework with custom markers for categorical test execution
  - Performance characterization with environment-optimized benchmarking
  - Cross-platform compatibility verification
- **Scientific Documentation Standards**:
  - Complete API reference with usage specifications
  - Detailed parameter optimization guidelines
  - Methodological examples with real-world genomic datasets
  - Performance tuning recommendations for various computational environments
- **Development Infrastructure**:
  - Multi-environment testing via tox automation
  - Code coverage analysis with configurable reporting thresholds
  - Makefile automation for streamlined development workflows
  - Standardized issue reporting and contribution templates

### Architectural Modifications
- Restructured codebase according to standard Python package conventions (`src/merpcr/` layout)
- Decomposed monolithic implementation into functionally distinct modules:
  - `core/engine.py` - Primary computational algorithms and search logic
  - `core/models.py` - Data structure definitions and type specifications
  - `core/utils.py` - Utility functions for sequence manipulation and mathematical operations
  - `io/fasta.py` - FASTA format parsing and validation
  - `io/sts.py` - STS marker file processing and hash table construction
  - `cli.py` - Command-line interface with me-PCR compatibility layer
- Organized testing framework within dedicated `tests/` directory with structured data management
- Refactored import dependencies to support modular architecture

### Resolved Issues
- Corrected STS record duplication artifacts affecting search result accuracy
- Resolved type system inconsistencies in sequence comparison algorithms
- Enhanced reverse complement search algorithm precision
- Achieved complete test suite validation with refactored architecture
- **Validation Framework Corrections**:
  - Corrected IUPAC ambiguity code complement calculations in test expectations
  - Aligned FASTA parsing test assertions with implementation specifications
  - Standardized error handling expectations for file system operations
  - Updated nucleotide filtering validation to encompass complete IUPAC specification
  - Optimized performance benchmarks for continuous integration environment constraints
  - Implemented comprehensive me-PCR output format compatibility validation

### Deprecated Components
- Eliminated obsolete development artifacts and legacy directory structures
- Removed deprecated packaging configurations (setup.py, requirements.txt) in favor of modern standards
- Eliminated redundant code segments and unused import dependencies
- Cleaned legacy testing data and superseded validation approaches