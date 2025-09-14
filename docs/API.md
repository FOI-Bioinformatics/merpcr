# API Technical Specification

This document provides comprehensive technical documentation for the merPCR software package interface.

## Primary Classes

### MerPCR

Principal computational engine class implementing STS marker detection algorithms for genomic sequence analysis.

#### Constructor

```python
MerPCR(
    wordsize: int = 11,
    margin: int = 50,
    mismatches: int = 0,
    three_prime_match: int = 1,
    threads: int = 1,
    quiet: int = 1,
    default_pcr_size: int = 240,
    iupac_mode: int = 0
)
```

**Parameter Specifications:**
- `wordsize` (int): Hash word length for initial sequence matching (valid range: 3-16). Default: 11
- `margin` (int): Maximum permitted deviation from expected PCR amplicon size (valid range: 0-10000). Default: 50
- `mismatches` (int): Permitted mismatch count in primer-template alignment (valid range: 0-10). Default: 0
- `three_prime_match` (int): Required exact match length at primer 3' terminus. Default: 1
- `threads` (int): Concurrent processing thread count for computational parallelization. Default: 1
- `quiet` (int): Logging verbosity control (0=detailed output, 1=minimal output). Default: 1
- `default_pcr_size` (int): Default amplicon size for undefined PCR products (valid range: 50-10000). Default: 240
- `iupac_mode` (int): IUPAC ambiguity nucleotide code processing (0=disabled, 1=enabled). Default: 0

**Exception Handling:**
- `ValueError`: Raised when input parameters exceed defined validation ranges

#### Methods

##### load_sts_file(filename: str) -> bool

Imports STS marker definitions from tab-delimited input files.

**Input Parameters:**
- `filename` (str): Filesystem path to STS marker definition file

**Return Values:**
- `bool`: Success indicator (True=successful import, False=import failure)

**File Format Specification:**
Tab-delimited format: `Identifier\tForward_Primer\tReverse_Primer\tAmplicon_Size\t[Optional_Annotation]`

##### load_fasta_file(filename: str) -> List[FASTARecord]

Imports genomic sequences from FASTA-formatted input files.

**Input Parameters:**
- `filename` (str): Filesystem path to FASTA sequence file

**Return Values:**
- `List[FASTARecord]`: Collection of parsed sequence records with associated metadata

##### search(fasta_records: List[FASTARecord], output_file: Optional[str] = None) -> int

Executes STS marker detection algorithms against target genomic sequences.

**Input Parameters:**
- `fasta_records` (List[FASTARecord]): Target genomic sequences for analysis
- `output_file` (Optional[str]): Output destination specification (None directs to stdout)

**Return Values:**
- `int`: Total count of detected STS marker matches

## Data Structure Specifications

### STSRecord

Data structure representing Sequence-Tagged Site marker definitions.

```python
@dataclass
class STSRecord:
    id: str                    # Unique STS identifier
    primer1: str              # Forward primer sequence (5' to 3')
    primer2: str              # Reverse primer sequence (5' to 3')
    pcr_size: int            # Predicted PCR amplicon size
    alias: str = ""          # Optional descriptive annotation
    offset: int = 0          # Positional offset within input file
    hash_offset: int = 0     # Hash computation positional offset
    direct: str = '+'        # Amplification strand orientation ('+' forward, '-' reverse)
    ambig_primer: int = 0    # IUPAC ambiguity presence indicator
```

### FASTARecord

Represents a FASTA sequence record.

```python
@dataclass
class FASTARecord:
    defline: str             # FASTA header line
    sequence: str            # DNA sequence
    label: str = ""          # Extracted sequence label
    
    def __post_init__(self):
        # Automatically extracts label from defline if not provided
```

### STSHit

Represents a found STS match.

```python
@dataclass
class STSHit:
    pos1: int                # Start position of match
    pos2: int                # End position of match  
    sts: STSRecord           # The STS that matched
```

### ThreadData

Data structure for multithreaded processing.

```python
@dataclass
class ThreadData:
    thread_id: int           # Thread identifier
    sequence: str            # Sequence segment
    offset: int              # Global offset
    length: int              # Segment length
    hits: List[STSHit] = field(default_factory=list)  # Found hits
```

## I/O Modules

### FASTALoader

Static methods for loading FASTA files.

#### load_file(filename: str) -> List[FASTARecord]

Load sequences from a FASTA file with proper error handling.

**Parameters:**
- `filename` (str): Path to FASTA file

**Returns:**
- `List[FASTARecord]`: List of parsed sequences

**Features:**
- Handles multiline sequences
- Filters invalid characters
- Preserves IUPAC ambiguity codes
- Robust error handling

## Internal Methods

The following methods are available but primarily intended for internal use:

### MerPCR Internal Methods

#### _hash_value(sequence: str) -> Tuple[int, int]

Compute hash value for sequence matching.

**Returns:**
- `Tuple[int, int]`: (offset, hash_value) or (-1, 0) if invalid

#### _compare_seqs(seq1: str, seq2: str, direction: str) -> bool

Compare two sequences allowing for mismatches and 3' protection.

**Parameters:**
- `seq1` (str): First sequence
- `seq2` (str): Second sequence  
- `direction` (str): Strand direction ('+' or '-')

**Returns:**
- `bool`: True if sequences match within tolerance

#### _reverse_complement(sequence: str) -> str

Generate reverse complement of a DNA sequence.

**Parameters:**
- `sequence` (str): Input DNA sequence

**Returns:**
- `str`: Reverse complement sequence

## Usage Examples

### Basic Usage

```python
from merpcr import MerPCR

# Create engine with default parameters
mer_pcr = MerPCR()

# Load data
success = mer_pcr.load_sts_file("markers.sts")
if not success:
    print("Failed to load STS file")
    exit(1)

sequences = mer_pcr.load_fasta_file("genome.fa")
if not sequences:
    print("No sequences loaded")
    exit(1)

# Search for matches
hit_count = mer_pcr.search(sequences, "results.txt")
print(f"Found {hit_count} hits")
```

### Advanced Usage

```python
from merpcr import MerPCR, FASTARecord

# Create engine with custom parameters
mer_pcr = MerPCR(
    wordsize=11,
    margin=100,
    mismatches=2,
    three_prime_match=2,
    threads=4,
    iupac_mode=1,
    quiet=0  # Verbose output
)

# Manual FASTA record creation
sequence = FASTARecord(
    defline=">chromosome1",
    sequence="ATCGATCGATCGATCG"
)

# Search with programmatically created sequences
mer_pcr.load_sts_file("markers.sts")
hits = mer_pcr.search([sequence])
```

### Error Handling

```python
from merpcr import MerPCR
import logging

# Enable logging to see detailed error messages
logging.basicConfig(level=logging.INFO)

try:
    # Parameters will be validated
    mer_pcr = MerPCR(wordsize=20)  # Will raise ValueError
except ValueError as e:
    print(f"Invalid parameter: {e}")

# File loading with error checking
mer_pcr = MerPCR()
if not mer_pcr.load_sts_file("nonexistent.sts"):
    print("STS file not found or invalid")

sequences = mer_pcr.load_fasta_file("nonexistent.fa")
if not sequences:
    print("FASTA file not found or empty")
```

## Performance Considerations

### Thread Usage

- For sequences < 100KB, single-threaded processing is used regardless of thread setting
- Optimal thread count is typically equal to the number of CPU cores
- Memory usage increases with thread count and word size

### Word Size Selection

- Larger word sizes (11-12) provide faster initial matching but use more memory
- Smaller word sizes (7-10) are more thorough but slower
- Word size must be ≤ shortest primer length

### Memory Usage

- Memory usage scales with: `O(4^wordsize + sequence_length + num_sts)`
- For large genomes, consider processing chromosome by chromosome
- Monitor memory usage with performance tests

## Compatibility

### me-PCR Compatibility

This implementation is designed to be fully compatible with the original me-PCR:

- Identical command-line interface
- Same output format
- Equivalent search algorithms
- Matching parameter validation

### Python Version Requirements

- Python 3.9+
- Standard library only (no external dependencies for core functionality)
- Optional: psutil for memory monitoring in tests