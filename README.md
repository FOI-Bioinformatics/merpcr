# merPCR - Modern Electronic Rapid PCR

Python reimplementation of the me-PCR (Multithreaded Electronic PCR) program originally developed by
Gregory Schuler at NCBI and enhanced by Kevin Murphy at Children's Hospital of Philadelphia.

## Description

merPCR is a tool for searching large DNA sequences for Sequence-Tagged Sites (STS) markers. STSs are defined as two short subsequences (primers) separated by an approximate distance. They are commonly used in genome mapping and chromosome identification.

The original me-PCR was developed as an enhanced version of the e-PCR program, adding multithreading capability and other improvements. This Python reimplementation (merPCR) retains all functionality while providing modern programming practices, a clean modular architecture, and comprehensive testing.

## ✅ Verification with me-PCR

This implementation has been thoroughly tested against the original me-PCR and produces identical results:

- **Algorithm accuracy**: Matches me-PCR output exactly on test data
- **Parameter compatibility**: All command-line options work identically  
- **Performance**: Comparable search performance with improved code maintainability
- **Comprehensive testing**: Includes both unit tests and integration tests with real genomic data

## Features

- **Multithreaded**: Uses Python's concurrent processing for faster analysis on multi-core systems
- **IUPAC Support**: Optional handling of IUPAC ambiguity codes in primers
- **Margin Control**: Adjust the allowed distance between primers
- **Mismatch Tolerance**: Configure the number of mismatches allowed
- **3' End Protection**: Control how many bases at the 3' end must match exactly
- **Modern Python**: Implemented using typing, argparse, and proper error handling

## Installation

### From Source

```bash
# Clone the repository
git clone https://github.com/yourusername/merpcr.git
cd merpcr

# Install in development mode
pip install -e .

# Or install directly
pip install .
```

### For Development

```bash
# Install with development dependencies
pip install -e ".[dev]"
```

## Usage

After installation, you can use merPCR in several ways:

### Command Line

```bash
# Using the installed command
merpcr [options] sts_file fasta_file

# Using the script directly  
./scripts/merpcr [options] sts_file fasta_file

# As a Python module
python -m merpcr [options] sts_file fasta_file
```

### Arguments

- `sts_file`: Tab-delimited file containing STS marker data
- `fasta_file`: FASTA-formatted sequence file

### Options

- `-M, --margin`: Margin (default: 50)
- `-N, --mismatches`: Number of mismatches allowed (default: 0)
- `-W, --wordsize`: Word size (default: 11)
- `-T, --threads`: Number of threads (default: 1)
- `-X, --three-prime-match`: Number of 3'-ward bases in which to disallow mismatches (default: 1)
- `-O, --output`: Output file name (default: stdout)
- `-Q, --quiet`: Quiet flag (0=verbose, 1=quiet)
- `-Z, --default-pcr-size`: Default PCR size when not specified in STS file (default: 240)
- `-I, --iupac`: IUPAC flag (0=don't honor IUPAC ambiguity symbols, 1=honor IUPAC symbols)
- `--debug`: Enable debug logging
- `-v, --version`: Show program's version number and exit

## STS File Format

The STS file should be a tab-delimited text file with the following format:

```
ID	Primer1	Primer2	PCR_Size [Alias]
```

Example:
```
AFM248yg9	GCTAAAAATACACGGATGG	TGCAAGACTGCGTCTC	193	(D17S932)  Chr.17, 63.7 cM
```

- **ID**: Unique identifier for the STS
- **Primer1**: First primer sequence
- **Primer2**: Second primer sequence
- **PCR_Size**: Expected PCR product size in base pairs (can be a range like "100-150")
- **Alias**: Optional additional information

## Output Format

The output is in the format:
```
SequenceLabel	Position1..Position2	STS_ID	(Orientation)
```

Where:
- **SequenceLabel**: Label from the FASTA sequence
- **Position1..Position2**: 1-based positions of the match
- **STS_ID**: ID of the matched STS
- **Orientation**: + for forward strand, - for reverse strand

## Tips for Optimal Use

1. **Word Size (W)**:
   - Larger word sizes (10-12) are faster but use more memory
   - W=11 offers good performance on modern systems
   - Use smaller W values for more thorough searching

2. **Margin (M)**:
   - Larger margins allow more flexibility in finding STSs
   - Values of 50-100 are typical
   - For STSs with unknown size, use larger margins

3. **Threads (T)**:
   - Set to the number of available CPU cores for best performance
   - Small sequences (<100KB) automatically use a single thread

4. **IUPAC Mode (I)**:
   - Enable (I=1) if your STSs contain ambiguity codes like 'N', 'W', etc.

## Example

```bash
merpcr -M 50 -N 1 -W 11 -T 4 -I 1 -O results.txt sts_markers.txt genome.fa
```

This runs merPCR with a margin of 50, allowing 1 mismatch, using word size 11, 4 threads, with IUPAC support, and writing results to results.txt.

### Python API

```python
from merpcr import MerPCR

# Initialize with custom parameters
mer_pcr = MerPCR(
    wordsize=11,
    margin=50,
    mismatches=1,
    threads=4,
    iupac_mode=1
)

# Load data and search
mer_pcr.load_sts_file("sts_markers.txt")
fasta_records = mer_pcr.load_fasta_file("genome.fa")
hit_count = mer_pcr.search(fasta_records, "results.txt")
```

## Testing

merPCR includes comprehensive test suites to ensure accuracy and reliability:

### Quick Test Commands

```bash
# Run all tests
make test

# Run unit tests only  
make test-unit

# Run integration tests
make test-integration

# Run performance benchmarks
make test-performance

# Generate coverage report
make coverage
```

### Using pytest directly

```bash
# All tests
pytest

# Specific test categories
pytest -m unit          # Unit tests only
pytest -m integration   # Integration tests only
pytest -m performance   # Performance tests only
pytest -m cli          # CLI tests only

# Verbose output
pytest -v

# With coverage
pytest --cov=src/merpcr --cov-report=html
```

### Test Categories

- **Unit Tests** (`tests/test_core_models.py`, `tests/test_engine_internals.py`): Test individual components and functions
- **Integration Tests** (`tests/test_comprehensive.py`): End-to-end testing with real data from me-PCR test suite
- **I/O Tests** (`tests/test_io_modules.py`): Test file loading and parsing functionality
- **CLI Tests** (`tests/test_cli.py`): Test command-line interface and parameter handling
- **Performance Tests** (`tests/test_performance.py`): Benchmarking and scalability testing

The comprehensive tests use actual data from the original me-PCR test suite and verify:
- Exact match with me-PCR output format
- Parameter validation and bounds checking
- IUPAC ambiguity code support
- Threading behavior and performance
- Memory efficiency
- Error handling and edge cases

## Troubleshooting

### Common Issues

1. **No hits found**:
   - Try increasing the margin (-M) parameter
   - Try allowing more mismatches (-N)
   - Check if IUPAC support (-I 1) is needed for your data

2. **Performance issues**:
   - Adjust word size (-W) according to available memory
   - Very large sequences may require splitting into smaller files
   - Make sure not to use more threads than you have CPU cores

### Debug Mode

Use the `--debug` flag to get detailed information about the run:

```bash
merpcr --debug sts_file.txt genome.fa
```

## Project Structure

```
merpcr/
├── src/merpcr/          # Main package
│   ├── __init__.py      # Package exports
│   ├── __main__.py      # Module entry point
│   ├── cli.py           # Command-line interface
│   ├── core/            # Core functionality
│   │   ├── engine.py    # Main search engine
│   │   ├── models.py    # Data models
│   │   └── utils.py     # Utility functions
│   └── io/              # Input/output modules
│       ├── fasta.py     # FASTA file handling
│       └── sts.py       # STS file handling
├── tests/               # Test suite
│   ├── data/            # Test data files
│   ├── test_basic.py    # Basic unit tests
│   └── test_comprehensive.py  # Integration tests
├── scripts/             # Entry point scripts
├── docs/                # Documentation
├── pyproject.toml       # Modern Python packaging
└── README.md            # This file
```

## License

This project is licensed under the GNU General Public License v3.0 - see the LICENSE file for details.

## Acknowledgments

- Gregory Schuler - Original e-PCR developer at NCBI
- Kevin Murphy - Developer of me-PCR at Children's Hospital of Philadelphia

## References

- Schuler, G.D. (1997) "Sequence mapping by electronic PCR." Genome Research 7: 541-550.