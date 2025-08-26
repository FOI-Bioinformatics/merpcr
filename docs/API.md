# API Reference

This document provides detailed API documentation for the merPCR package.

## Core Classes

### MerPCR

The main search engine class that handles STS marker searching in genomic sequences.

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

**Parameters:**
- `wordsize` (int): Size of the hash word for initial matching (3-16). Default: 11
- `margin` (int): Maximum allowed distance deviation from expected PCR size (0-10000). Default: 50
- `mismatches` (int): Number of mismatches allowed in primer matching (0-10). Default: 0
- `three_prime_match` (int): Number of bases at 3' end that must match exactly. Default: 1
- `threads` (int): Number of threads to use for processing. Default: 1
- `quiet` (int): Verbosity level (0=verbose, 1=quiet). Default: 1
- `default_pcr_size` (int): Default PCR product size when not specified (50-10000). Default: 240
- `iupac_mode` (int): Whether to honor IUPAC ambiguity codes (0=no, 1=yes). Default: 0

**Raises:**
- `ValueError`: If any parameter is outside valid range

#### Methods

##### load_sts_file(filename: str) -> bool

Load STS markers from a tab-delimited file.

**Parameters:**
- `filename` (str): Path to the STS file

**Returns:**
- `bool`: True if loading succeeded, False otherwise

**File Format:**
Each line should contain: `ID\tPrimer1\tPrimer2\tPCR_Size\t[Alias]`

##### load_fasta_file(filename: str) -> List[FASTARecord]

Load sequences from a FASTA file.

**Parameters:**
- `filename` (str): Path to the FASTA file

**Returns:**
- `List[FASTARecord]`: List of loaded FASTA records

##### search(fasta_records: List[FASTARecord], output_file: Optional[str] = None) -> int

Search for STS matches in the given sequences.

**Parameters:**
- `fasta_records` (List[FASTARecord]): Sequences to search
- `output_file` (Optional[str]): Output file path (None for stdout)

**Returns:**
- `int`: Total number of hits found

## Data Models

### STSRecord

Represents a Sequence-Tagged Site marker.

```python
@dataclass
class STSRecord:
    id: str                    # Unique identifier
    primer1: str              # First primer sequence
    primer2: str              # Second primer sequence
    pcr_size: int            # Expected PCR product size
    alias: str = ""          # Optional alias/description
    offset: int = 0          # Offset in sequence
    hash_offset: int = 0     # Hash computation offset
    direct: str = '+'        # Strand direction ('+' or '-')
    ambig_primer: int = 0    # IUPAC ambiguity flag
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