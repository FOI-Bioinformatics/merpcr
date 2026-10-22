# CLAUDE.md - merPCR Production System Guide

This file provides comprehensive guidance to Claude Code (claude.ai/code) when working with the merPCR codebase - a production-ready Python reimplementation of the me-PCR (Multithreaded Electronic PCR) program.

## System Status: PRODUCTION READY ✅

- **Version**: 1.0.0
- **Test Coverage**: **94% for engine.py** (critical component), comprehensive edge case testing
- **Test Suite**: 277 comprehensive tests across 19 test files
- **Architecture**: Modular, typed, thread-safe
- **Compatibility**: me-PCR output compatibility verified on 42 real genomes and by differential fuzzing against the me-PCR 1.0.6 binary (scripts/compare_with_mepcr.py)
- **Performance**: Optimized multithreading for large datasets
- **Critical Fixes**: 3 algorithmic differences fixed (see VERIFICATION.md), plus the 2026-10 binary audit (alias column, FASTA alphabet, one-directional IUPAC matching, scan-order reporting, overlap deduplication; see CHANGELOG.md)

## Critical Compatibility Fixes ⚠️

Three critical algorithmic differences were identified and fixed to achieve 100% compatibility with me-PCR:

### 1. Hash Computation (backward vs forward search)
- **Issue:** merPCR searched from position 0 forward; me-PCR searches backward from end
- **Impact:** Different hash values → different search results
- **Status:** ✅ FIXED in engine.py:_hash_value()
- **Tests:** tests/test_hash_computation_equivalence.py (8 tests)

### 2. PCR Range Margin Adjustment
- **Issue:** PCR ranges like "200-220" didn't adjust margin properly
- **Impact:** Incorrect search window for size ranges
- **Status:** ✅ FIXED in engine.py:_parse_pcr_size()
- **Tests:** tests/test_pcr_range_handling.py (9 tests)

### 3. Forward Strand Primer Matching (MOST CRITICAL)
- **Issue:** primer2 stored as-is instead of reverse complement for forward direction
- **Impact:** **ALL forward strand (+) hits were missed** - only found reverse strand
- **Example:** F. tularensis had 7 hits before fix, 20 hits after (13 missing forward hits)
- **Status:** ✅ FIXED in engine.py:load_sts_file()
- **Tests:** tests/test_forward_strand_matching.py (8 tests)

Full details: [docs/VERIFICATION.md](docs/VERIFICATION.md). Documentation files: docs/USER_GUIDE.md, docs/API.md, docs/VERIFICATION.md, docs/PERFORMANCE.md, docs/DEVELOPMENT.md (CI and release process).

---

## Development Commands

### Core Testing
- `make test` or `pytest` - Run all 277 tests
- `make test-unit` - Unit tests only
- `make test-integration` - Integration tests
- `make test-performance` - Performance benchmarks
- `make coverage` - Generate coverage report (94% for engine.py achieved!)
- `pytest -m "not slow"` - Skip slow tests during development

### Advanced Testing
- `pytest tests/test_property_based.py` - Hypothesis-based robustness testing
- `pytest tests/test_threading_stress.py` - Concurrency stress tests
- `pytest tests/test_error_injection.py` - Fault tolerance testing
- `pytest --timeout=300` - Set test timeouts for CI/CD

### Code Quality & Formatting
- `make lint` - Run flake8, black --check, isort --check
- `make format` - Auto-format with black (100-char line length)
- `mypy src/` - Type checking (strict mode enabled)
- `ruff check src/` - Additional linting

### Development Workflow
- `make dev-install` - Install in development mode
- `make clean` - Clean all artifacts
- `make build` - Build distribution packages
- `make upload` - Upload to PyPI

### Compatibility Testing
- `python test_compatibility.py` - Verify me-PCR output compatibility
- `python scripts/compare_with_mepcr.py --mepcr /path/to/me-PCR -n 50` - Differential fuzzing against the me-PCR binary (add `--pure-python` to test the fallback path)
- Deliberate deviations from me-PCR 1.0.6: size ranges and "-" are accepted (the binary rejects them); STS lines after a comment are reported correctly (the binary misreports them)
- Test both argument formats: `-M 50` and `M=50`

## Production Architecture

### Overview
merPCR is a high-performance bioinformatics tool for electronic PCR (e-PCR) analysis, searching genomic sequences for Sequence-Tagged Sites (STS) markers using primer pairs. It achieves full compatibility with the original me-PCR while providing modern Python architecture.

### Core Package Structure
```
src/merpcr/
├── __init__.py           # Package initialization, exports MerPCR, STSRecord, FASTARecord, STSHit
├── __main__.py           # Module entry point (python -m merpcr)
├── cli.py                # Command-line interface with me-PCR compatibility
└── core/
    ├── engine.py         # Main MerPCR class - search algorithm with threading
    ├── models.py         # Data models (STSRecord, FASTARecord, STSHit, ThreadData)
    └── utils.py          # Utility functions (reverse_complement, hash_value, IUPAC)
└── io/
    ├── fasta.py          # FASTA file loader with validation
    └── sts.py            # STS file loader with hash table construction
```

### Key Production Features

#### 1. me-PCR Full Compatibility
- **Argument Compatibility**: Supports both modern (`-M 50`) and legacy (`M=50`) formats
- **Output Compatibility**: Identical output format including alias fields
- **Parameter Compatibility**: All original parameters with validation
- **File Format Support**: Standard STS and FASTA formats with error handling

#### 2. High-Performance Search Engine (`engine.py`)
- **Multithreading**: Automatic threading for files >100KB with ThreadPoolExecutor
- **Hash-Based Lookup**: O(1) STS lookup using 2-bit encoded hash tables
- **Memory Efficient**: Streaming file processing for large datasets
- **IUPAC Support**: Full IUPAC ambiguity code handling
- **Bidirectional Search**: Forward and reverse complement primer matching

#### 3. Robust Data Models (`models.py`)
- **STSRecord**: STS marker with primer sequences, PCR size, aliases
- **FASTARecord**: Genomic sequence with metadata extraction
- **STSHit**: Search result with position and match details
- **ThreadData**: Thread-safe data containers for parallel processing

#### 4. Production I/O Handling
- **FASTA Loader**: Handles multi-sequence files, validates nucleotide characters
- **STS Loader**: Tab-delimited parsing with range support, error reporting
- **File Validation**: Size checks, format validation, graceful error handling
- **Large File Support**: Efficient processing of GB-scale genomic files

#### 5. Comprehensive CLI (`cli.py`)
- **Dual Format Support**: `-M 50` (modern) and `M=50` (legacy me-PCR)
- **Parameter Validation**: Range checking with descriptive error messages
- **Logging System**: Debug, info, warning levels with timestamp
- **Output Control**: File or stdout with proper buffering

### Search Algorithm Details

#### Hash-Based Primer Lookup
1. **Hash Computation**: 2-bit encoding (A=0, C=1, G=2, T=3) for word-size k-mers
2. **Collision Handling**: Hash table with chaining for multiple STS per hash
3. **Bidirectional**: Forward primer hash + reverse complement second primer
4. **IUPAC Handling**: Skip ambiguous positions during hash computation

#### Threading Architecture
- **Automatic Scaling**: Thread count based on file size and CPU cores
- **Work Distribution**: Sequence chunks with overlap handling
- **Thread Safety**: Immutable data structures, atomic hit counting
- **Load Balancing**: Dynamic chunk sizing for optimal CPU utilization

#### Sequence Comparison
- **Mismatch Tolerance**: Configurable mismatch count (0-10)
- **3' Protection**: Prevents mismatches in 3'-ward bases
- **Case Insensitive**: Handles mixed case input sequences
- **Length Validation**: Ensures primer compatibility with PCR size

### Production Configuration

#### System Requirements
- **Python**: 3.11+ (tested on 3.11-3.13)
- **Memory**: Scales with dataset size, ~1GB for human genome
- **CPU**: Multithreading benefits from 2+ cores
- **Storage**: Handles files up to several GB

#### Default Parameters (me-PCR Compatible)
- **Word Size**: 11 (range: 3-16)
- **Margin**: 50bp (range: 0-10000)
- **Mismatches**: 0 (range: 0-10)
- **3' Protection**: 1bp (minimum: 0)
- **PCR Size**: 240bp (range: 1-10000)
- **Threads**: 1 (auto-scaling available)
- **IUPAC Mode**: Disabled (0=off, 1=on)

### Test Suite Architecture (94% Coverage for engine.py - Excellent!)

#### Core Test Categories
1. **Unit Tests** (`test_basic.py`, `test_utils_comprehensive.py`)
   - Individual function testing
   - Boundary condition validation
   - Error handling verification

2. **Integration Tests** (`test_comprehensive.py`, `test_io_modules.py`)
   - End-to-end workflow testing
   - File format compatibility
   - Real-world data processing

3. **CLI Tests** (`test_cli.py`, `test_cli_enhanced.py`, `test_module_entry_point.py`)
   - Argument parsing (both formats)
   - Error handling and validation
   - Module entry point testing

4. **Compatibility Tests** (Critical for me-PCR compatibility)
   - **Hash Computation** (`test_hash_computation_equivalence.py`): Backward search validation (8 tests)
   - **PCR Range Handling** (`test_pcr_range_handling.py`): Margin adjustment tests (9 tests)
   - **Forward Strand Matching** (`test_forward_strand_matching.py`): Primer RC storage tests (8 tests)

5. **Edge Case Tests** (New - Comprehensive coverage improvement)
   - **Engine Edge Cases** (`test_engine_edge_cases.py`): 31 tests targeting uncovered code paths
     - Parameter validation edge cases
     - STS loading error paths
     - PCR size parsing edge cases
     - Hash computation edge cases
     - Threading and multiprocessing paths
     - IUPAC mode functionality
     - Match STS edge cases

6. **Advanced Testing** (Production Features)
   - **Property-Based** (`test_property_based.py`): Hypothesis-generated test cases
   - **Threading Stress** (`test_threading_stress.py`): Concurrency validation
   - **Error Injection** (`test_error_injection.py`): Fault tolerance testing
   - **Performance** (`test_performance.py`): Benchmark validation

#### Test Coverage Breakdown (Current)
- `core/engine.py`: **94%** (20/344 lines missed) - **Excellent!** ✅
  - Improved from 62% by adding comprehensive edge case tests
  - Remaining missed lines: 199-200, 222, 227-230, 335-340, 454-455, 493, 533, 540, 605, 625, 635-643, 667
- `cli.py`: 98% (2/120 lines missed)
- `core/utils.py`: 100% (0/38 lines missed)
- `io/fasta.py`: 86% (5/36 lines missed)
- `core/models.py`: 98% (1/41 lines missed)

**Achievement:** 32 percentage point improvement in engine.py coverage (62% → 94%)

### API Reference

#### Primary Class: MerPCR
```python
from merpcr import MerPCR

# Initialize with parameters
engine = MerPCR(
    wordsize=11,           # Hash word size (3-16)
    margin=50,             # Search margin in bp (0-10000)
    mismatches=0,          # Allowed mismatches (0-10)
    three_prime_match=1,   # 3' protection bases (≥0)
    iupac_mode=0,          # IUPAC ambiguity handling (0/1)
    default_pcr_size=240,  # Default PCR size (1-10000)
    threads=1              # Thread count (≥1)
)

# Load data files
success = engine.load_sts_file("primers.sts")
records = engine.load_fasta_file("genome.fa")

# Perform search
hit_count = engine.search(records, output_file="results.txt")
```

#### Data Models
```python
from merpcr import STSRecord, FASTARecord, STSHit

# STS marker definition
sts = STSRecord(
    id="STS_001",
    primer1="ATCGATCGATCG",
    primer2="GCTAGCTAGCTA",
    pcr_size=200,
    alias="Test STS"
)

# Genomic sequence
seq = FASTARecord(
    defline=">chr1 Human chromosome 1",
    sequence="ATCGATCG..."
)

# Search result
hit = STSHit(pos1=1000, pos2=1200, sts=sts)
```

### Performance Characteristics

#### Benchmarks (Typical Hardware)
- **Small Dataset** (<1MB): <1 second, single-threaded
- **Medium Dataset** (10-100MB): 10-60 seconds, multithreaded
- **Large Dataset** (>100MB): Scales linearly with threading
- **Memory Usage**: ~2-3x file size during processing
- **Threading Benefit**: 2-8x speedup on multi-core systems

#### Optimization Features
- **Streaming Processing**: Constant memory usage for large files
- **Hash Table Efficiency**: O(1) average lookup time
- **Thread Pool Management**: Optimal core utilization
- **Lazy Loading**: On-demand sequence processing

### Development Workflows

#### Adding New Features
1. Write failing tests first (TDD approach)
2. Implement feature with type hints
3. Run full test suite (`make test`)
4. Check coverage (`make coverage`)
5. Format code (`make format`)
6. Update documentation

#### Bug Fixing
1. Create reproduction test case
2. Fix implementation
3. Verify fix with stress tests
4. Run compatibility tests
5. Update CHANGELOG.md

#### Performance Optimization
1. Profile with `test_performance.py`
2. Identify bottlenecks
3. Implement optimizations
4. Verify with threading stress tests
5. Benchmark against me-PCR for compatibility

### Production Deployment

#### Installation Options
```bash
# Production installation
pip install merpcr

# Development installation
git clone <repo>
cd merpcr
make dev-install

# Container deployment
docker build -t merpcr .
docker run merpcr input.sts genome.fa
```

#### Usage Examples
```bash
# Modern format
merpcr primers.sts genome.fa -M 50 -N 1 -W 11

# Legacy me-PCR format
merpcr primers.sts genome.fa M=50 N=1 W=11

# High-throughput processing
merpcr large_primers.sts human_genome.fa -T 8 -O results.txt

# Debug mode
merpcr primers.sts genome.fa --debug -Q 0
```