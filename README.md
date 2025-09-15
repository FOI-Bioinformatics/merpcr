<div align="center">
  <img src="assets/images/merPCR_logo.png" alt="merPCR Logo" width="400">
  
  # merPCR - Modern Electronic PCR Implementation
  
  A Python reimplementation of the me-PCR (Multithreaded Electronic PCR) program, originally developed by Gregory Schuler at NCBI and subsequently enhanced by Kevin Murphy at Children's Hospital of Philadelphia.
</div>

## Overview

merPCR provides computational tools for locating Sequence-Tagged Sites (STS) within large genomic sequences. STS markers consist of two short primer subsequences separated by a defined genomic interval and serve as fundamental components in genome mapping and chromosomal analysis.

This implementation maintains full compatibility with the original me-PCR while incorporating modern software engineering practices, including comprehensive type safety, modular architecture, and extensive validation testing.

## Validation and Compatibility

Extensive validation against the original me-PCR demonstrates complete functional equivalence:

- **Algorithmic fidelity**: Produces identical results across diverse genomic datasets
- **Parameter consistency**: Maintains complete command-line interface compatibility
- **Performance characteristics**: Achieves comparable computational efficiency with enhanced maintainability
- **Quality assurance**: Incorporates comprehensive unit and integration testing with authentic genomic data

## Key Features

- **Concurrent processing**: Utilizes Python's threading capabilities for enhanced performance on multi-core architectures
- **IUPAC compatibility**: Provides optional support for International Union of Pure and Applied Chemistry nucleotide ambiguity codes
- **Configurable search parameters**: Allows adjustment of primer separation distances and mismatch tolerance
- **Primer specificity control**: Implements 3'-end protection mechanisms for enhanced specificity
- **Modern software architecture**: Employs comprehensive type annotation, robust error handling, and modular design principles

## Installation

### From Source

```bash
# Clone the repository
git clone https://github.com/FOI-Bioinformatics/merpcr.git
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

## Input File Specifications

### STS Marker Format

STS input files require tab-delimited text format with the following structure:

```
Identifier	Forward_Primer	Reverse_Primer	Amplicon_Size	[Optional_Annotation]
```

Example entry:
```
AFM248yg9	GCTAAAAATACACGGATGG	TGCAAGACTGCGTCTC	193	(D17S932)  Chr.17, 63.7 cM
```

Field specifications:
- **Identifier**: Unique STS designation
- **Forward_Primer**: 5' to 3' sequence of the forward amplification primer
- **Reverse_Primer**: 5' to 3' sequence of the reverse amplification primer
- **Amplicon_Size**: Expected PCR product length in base pairs (accepts range notation: "150-200")
- **Optional_Annotation**: Additional metadata or mapping information

### Output Format Specification

Results are presented in tab-delimited format:
```
Sequence_Identifier	Genomic_Coordinates	STS_Identifier	Strand_Orientation
```

Field descriptions:
- **Sequence_Identifier**: FASTA header designation for the target sequence
- **Genomic_Coordinates**: 1-indexed positional boundaries of the predicted amplicon
- **STS_Identifier**: Corresponding STS marker designation from input file
- **Strand_Orientation**: Amplification directionality (+ indicates forward strand, - indicates reverse complement)

## Parameter Optimization Guidelines

### Hash Word Size (-W)
- **Performance consideration**: Larger values (10-12) improve computational speed but increase memory requirements
- **Recommended setting**: W=11 provides optimal performance on contemporary hardware
- **Sensitivity adjustment**: Smaller values enhance detection sensitivity at computational cost

### Search Margin (-M)
- **Flexibility control**: Increased values accommodate greater variance in amplicon size
- **Typical range**: Values between 50-100 base pairs are commonly employed
- **Unknown amplicon sizes**: Larger margins recommended for STS markers with undefined size ranges

### Thread Utilization (-T)
- **Performance optimization**: Configure to match available CPU core count for maximum throughput
- **Automatic optimization**: Sequences below 100KB automatically utilize single-thread processing

### IUPAC Ambiguity Support (-I)
- **Enable when necessary**: Activate (I=1) for primer sequences containing degenerate nucleotides (N, W, R, Y, etc.)
- **Default setting**: Standard nucleotides (A, T, G, C) require no special handling

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

## Quality Assurance and Validation

merPCR incorporates extensive testing frameworks to ensure computational accuracy and software reliability:

### Test Suite Execution

```bash
# Complete validation suite
make test

# Component-specific testing
make test-unit          # Individual function validation
make test-integration   # End-to-end workflow testing
make test-performance   # Computational benchmarking

# Code coverage analysis
make coverage
```

### Direct pytest Interface

```bash
# Comprehensive testing
pytest

# Targeted test categories
pytest -m unit          # Isolated component testing
pytest -m integration   # System integration validation
pytest -m performance   # Performance characterization
pytest -m cli          # Command-line interface testing

# Detailed reporting
pytest -v --cov=src/merpcr --cov-report=html
```

### Validation Framework Categories

- **Unit Testing**: Validates individual algorithmic components and data structures
- **Integration Testing**: Confirms end-to-end functionality using authentic genomic datasets
- **I/O Validation**: Verifies file parsing and data serialization accuracy
- **Interface Testing**: Ensures command-line parameter processing and error handling
- **Performance Benchmarking**: Characterizes computational efficiency and scalability

**Validation Standards**: All tests utilize reference datasets from the original me-PCR implementation, confirming:
- Exact algorithmic equivalence with legacy output formats
- Comprehensive parameter boundary validation
- IUPAC ambiguity code processing accuracy
- Multi-threading stability and performance characteristics
- Memory utilization efficiency
- Robust error handling across edge cases

## Diagnostic Procedures and Troubleshooting

### Common Analytical Challenges

**Low Detection Sensitivity**:
- Increase search margin parameters (-M) to accommodate amplicon size variance
- Adjust mismatch tolerance (-N) to account for sequence polymorphisms
- Verify IUPAC support activation (-I 1) for degenerate primer sequences

**Performance Optimization**:
- Configure hash word size (-W) based on available system memory
- Consider sequence partitioning for exceptionally large genomic files
- Ensure thread count (-T) does not exceed available CPU cores

### Diagnostic Mode

Detailed execution logging is available through debug mode:

```bash
merpcr --debug input_markers.sts target_sequence.fa
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

We acknowledge the foundational contributions of:
- Gregory D. Schuler (NCBI) - Original e-PCR algorithm development
- Kevin Murphy (Children's Hospital of Philadelphia) - me-PCR multithreading enhancements

## References

1. Schuler, G.D. (1997) "Sequence mapping by electronic PCR." *Genome Research* **7**: 541-550. doi:10.1101/gr.7.5.541

2. Altschul, S.F., Gish, W., Miller, W., Myers, E.W., and Lipman, D.J. (1990) "Basic local alignment search tool." *Journal of Molecular Biology* **215**: 403-410. doi:10.1016/S0022-2836(05)80360-2
