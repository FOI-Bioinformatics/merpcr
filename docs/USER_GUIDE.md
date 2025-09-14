# merPCR User Guide

This guide provides detailed instructions for using merPCR effectively for genomic STS marker searching.

## Table of Contents

1. [Getting Started](#getting-started)
2. [Understanding STS Markers](#understanding-sts-markers)
3. [Preparing Your Data](#preparing-your-data)
4. [Command Line Usage](#command-line-usage)
5. [Parameter Tuning](#parameter-tuning)
6. [Interpreting Results](#interpreting-results)
7. [Common Workflows](#common-workflows)
8. [Troubleshooting](#troubleshooting)
9. [Performance Optimization](#performance-optimization)

## Getting Started

### Quick Start

1. **Install merPCR**:
   ```bash
   pip install -e .
   ```

2. **Prepare your data**:
   - STS marker file (tab-delimited)
   - FASTA sequence file

3. **Run a basic search**:
   ```bash
   merpcr markers.sts genome.fa
   ```

### Verifying Installation

Test your installation with the provided test data:

```bash
# Run basic functionality test
make test-unit

# Run comprehensive integration tests
make test-integration
```

## Understanding STS Markers

### What are STS Markers?

Sequence-Tagged Sites (STS) are short DNA sequences that occur only once in the genome and can be easily detected by PCR. They consist of:

- **Two primers**: Short DNA sequences (typically 18-24 nucleotides)
- **Expected distance**: The anticipated size of the PCR product
- **Unique identifier**: A name or ID for the marker

### STS Applications

- **Genome mapping**: Physical mapping of chromosomes
- **Genetic linkage**: Tracking inheritance patterns
- **Comparative genomics**: Cross-species sequence analysis
- **Quality control**: Verifying sequence assemblies

## Preparing Your Data

### STS Marker File Specification

Construct tab-delimited input files using the following format:

```
Identifier	Forward_Primer	Reverse_Primer	Amplicon_Size	[Annotation]
```

**Representative Examples**:
```
D1S2893	TGAGTCAGATGTTTGATTTTG	ATGCCACATCAACTTATACTG	126	1p36.33
D1S468	GAATGAACAGAGATGATGCCT	CACACACACACACACCACAC	142-148	1p36.22
AFM248yg9	GCTAAAAATACACGGATGG	TGCAAGACTGCGTCTC	193	(D17S932) Chr.17, 63.7 cM
```

**Field Specifications**:
- **Identifier**: Unique STS designation (alphanumeric characters and underscores)
- **Forward_Primer**: 5' to 3' sequence of the forward amplification primer (IUPAC nucleotide codes supported)
- **Reverse_Primer**: 5' to 3' sequence of the reverse amplification primer
- **Amplicon_Size**: Predicted PCR product size in base pairs
  - Single value: `200`
  - Range notation: `150-250` (system calculates mean value: 200)
- **Annotation**: Optional descriptive metadata (supports spaces and special characters)

**Format Requirements**:
- Utilize tab characters for field separation (not spaces)
- Lines prefixed with `#` are interpreted as comments
- Empty lines are ignored during processing
- Primer sequences must be specified in 5' to 3' orientation
- Amplicon size excludes primer sequences themselves

### FASTA File Format

Standard FASTA format is supported:

```
>sequence_identifier Description
ATCGATCGATCGATCG
GGCCTTAAGGCCTTAA
>another_sequence
AAAATTTTCCCCGGGG
```

**Requirements**:
- Standard FASTA headers starting with `>`
- DNA sequences using A, T, C, G
- IUPAC ambiguity codes (N, R, Y, etc.) are supported when IUPAC mode is enabled
- Sequences can span multiple lines
- No maximum sequence length

## Command Line Usage

### Basic Syntax

```bash
merpcr [OPTIONS] sts_file fasta_file
```

### Essential Options

#### Search Parameters
- `-M, --margin INT`: Distance tolerance around expected PCR size (default: 50)
  - Example: If PCR size is 200 and margin is 50, matches from 150-250 bp are accepted

- `-N, --mismatches INT`: Number of mismatches allowed in primers (default: 0)
  - Each primer can have up to N mismatches
  - Mismatches in 3' region may be restricted (see `-X`)

- `-W, --wordsize INT`: Hash word size for initial matching (3-16, default: 11)
  - Larger values: faster but more memory usage
  - Must be ≤ length of shortest primer

#### Primer Matching
- `-X, --three-prime-match INT`: 3' bases requiring exact match (default: 1)
  - Forward primer: protects rightmost X bases
  - Reverse primer: protects leftmost X bases (5' of reverse complement)
  - Set to 0 to disable 3' protection

- `-I, --iupac INT`: Handle IUPAC ambiguity codes (0=no, 1=yes, default: 0)
  - When enabled: N matches any base, R matches A/G, etc.
  - Useful for degenerate primers

#### Output Options
- `-O, --output FILE`: Write results to file instead of stdout
- `-Q, --quiet INT`: Verbosity level
  - 0: Verbose (shows progress, statistics)
  - 1: Quiet (results only, default)

#### Performance Options
- `-T, --threads INT`: Number of processing threads (default: 1)
  - Set to number of CPU cores for best performance
  - Automatically limited to 1 for small sequences (<100KB)

- `-Z, --default-pcr-size INT`: Default PCR size when not specified (default: 240)

### Advanced Options

- `--debug`: Enable detailed logging and debugging output
- `-v, --version`: Show version information
- `--help`: Show all options and usage examples

## Parameter Tuning

### Choosing the Right Parameters

#### Word Size (-W)
```bash
# Conservative search (slower, more thorough)
merpcr -W 8 markers.sts genome.fa

# Standard search (balanced)
merpcr -W 11 markers.sts genome.fa

# Fast search (faster, uses more memory)
merpcr -W 14 markers.sts genome.fa
```

**Guidelines**:
- W=8-10: Use for highly divergent sequences or short primers
- W=11-12: Good balance for most applications
- W=13-16: Use for very similar sequences or when speed is critical

#### Margin (-M)
```bash
# Strict size matching
merpcr -M 10 markers.sts genome.fa

# Standard tolerance  
merpcr -M 50 markers.sts genome.fa

# Flexible size matching
merpcr -M 200 markers.sts genome.fa
```

**Guidelines**:
- Small margin (10-25): When PCR sizes are well-characterized
- Medium margin (50-100): Standard for most applications
- Large margin (100+): For poorly characterized STSs or divergent species

#### Mismatches (-N)
```bash
# Exact matching only
merpcr -N 0 markers.sts genome.fa

# Allow some variation
merpcr -N 1 markers.sts genome.fa

# Highly permissive
merpcr -N 3 markers.sts genome.fa
```

**Guidelines**:
- N=0: High-quality sequences, exact matches needed
- N=1-2: Standard for most genomic searches
- N=3+: Highly divergent sequences or cross-species analysis

### Common Parameter Combinations

#### High-Specificity Search
```bash
# Strict parameters for high-confidence results
merpcr -M 25 -N 0 -X 2 -W 11 markers.sts genome.fa
```

#### Sensitive Search
```bash
# Permissive parameters to find more hits
merpcr -M 100 -N 2 -X 0 -W 8 -I 1 markers.sts genome.fa
```

#### Performance-Optimized
```bash
# Fast search with multiple threads
merpcr -W 12 -T 8 -Q 1 markers.sts genome.fa
```

## Interpreting Results

### Output Format

Each hit is reported in this format:
```
SequenceLabel	Position1..Position2	STS_ID	(Direction)
```

**Example Output**:
```
chr1	1000..1200	D1S2893	(+)
chr1	5000..5150	D1S468	(-)
chr17	30000..30193	AFM248yg9	(+)
```

**Field Descriptions**:
- **SequenceLabel**: Name from FASTA header (first word after `>`)
- **Position1**: Start position of the match (1-based)
- **Position2**: End position of the match (1-based)  
- **STS_ID**: Identifier from the STS file
- **Direction**: 
  - `(+)`: Forward strand match
  - `(-)`: Reverse strand match

### Understanding Positions

Positions represent the **outermost boundaries** of the PCR product:

```
Sequence: ATCG[PRIMER1]NNNNNNNN[PRIMER2]GCTA
          1234 5     12      20     27 30
```

For a match with primers at positions 5-12 and 20-27:
- **Position1**: 5 (start of first primer)
- **Position2**: 27 (end of second primer)
- **PCR Product Length**: 27 - 5 + 1 = 23 bp

### Strand Interpretation

#### Forward Strand (+)
- Primer1 found at lower position
- Primer2 found at higher position  
- Natural PCR orientation

#### Reverse Strand (-)
- Primer1 (reverse complement) found at higher position
- Primer2 (reverse complement) found at lower position
- PCR would amplify the reverse complement

## Common Workflows

### Workflow 1: Basic Genomic Screening

```bash
# 1. Prepare standard parameters
merpcr -M 50 -N 1 -W 11 -T 4 -O results.txt markers.sts genome.fa

# 2. Review results
head -20 results.txt

# 3. Count total hits
wc -l results.txt
```

### Workflow 2: Cross-Species Analysis

```bash
# Use permissive parameters for divergent sequences
merpcr -M 100 -N 3 -X 0 -I 1 -W 9 -O cross_species.txt \
       human_markers.sts mouse_genome.fa
```

### Workflow 3: Quality Control

```bash
# Strict search to verify expected markers
merpcr -M 25 -N 0 -X 2 -O qc_results.txt \
       known_markers.sts assembly.fa

# Compare with expected results
diff expected_results.txt qc_results.txt
```

### Workflow 4: Large-Scale Analysis

```bash
# Process multiple chromosomes efficiently
for chr in chr{1..22} chrX chrY; do
    merpcr -M 50 -N 1 -T 8 -O ${chr}_results.txt \
           markers.sts ${chr}.fa
done

# Combine results
cat chr*_results.txt > all_results.txt
```

## Troubleshooting

### No Results Found

**Possible Causes and Solutions**:

1. **Margin too small**:
   ```bash
   # Try larger margin
   merpcr -M 200 markers.sts genome.fa
   ```

2. **No mismatches allowed**:
   ```bash
   # Allow some variation
   merpcr -N 2 markers.sts genome.fa
   ```

3. **3' end protection too strict**:
   ```bash
   # Reduce or disable 3' protection
   merpcr -X 0 markers.sts genome.fa
   ```

4. **IUPAC codes in primers**:
   ```bash
   # Enable IUPAC mode
   merpcr -I 1 markers.sts genome.fa
   ```

### Too Many Results

**Solutions**:

1. **Tighten parameters**:
   ```bash
   # More stringent search
   merpcr -M 25 -N 0 -X 2 markers.sts genome.fa
   ```

2. **Check for repetitive sequences**:
   - Review primer specificity
   - Use longer or more specific primers
   - Filter results by position

### Performance Issues

**Memory Problems**:
```bash
# Reduce word size to save memory
merpcr -W 9 markers.sts genome.fa

# Process smaller sequence files
```

**Slow Performance**:
```bash
# Increase word size (if memory allows)
merpcr -W 13 -T 8 markers.sts genome.fa

# Use multiple threads
merpcr -T $(nproc) markers.sts genome.fa
```

### File Format Errors

**STS File Issues**:
- Check for consistent tab-delimited format
- Verify all required fields are present
- Remove or comment out problematic lines with `#`

**FASTA File Issues**:
- Ensure proper FASTA format with `>` headers
- Check for non-DNA characters in sequences
- Verify file encoding (should be ASCII/UTF-8)

### Debugging

Enable debug mode for detailed information:

```bash
merpcr --debug -Q 0 markers.sts genome.fa 2> debug.log
```

This will show:
- Parameter validation
- File loading progress  
- STS processing details
- Search statistics
- Error messages

## Performance Optimization

### System Requirements

**Minimum**:
- 4 GB RAM
- Single CPU core
- 1 GB free disk space

**Recommended**:
- 16+ GB RAM for large genomes
- 4+ CPU cores for multithreading
- SSD storage for large sequence files

### Memory Usage Estimation

Approximate memory usage:
```
Memory ≈ 4^W * 8 bytes + Sequence_Size + STS_Count * 200 bytes
```

**Examples**:
- W=11, 1GB sequence, 1000 STSs: ~4.2 GB RAM
- W=13, 100MB sequence, 100 STSs: ~1.1 GB RAM

### Performance Tips

1. **Optimize word size**:
   ```bash
   # Test different word sizes
   time merpcr -W 10 markers.sts genome.fa
   time merpcr -W 11 markers.sts genome.fa
   time merpcr -W 12 markers.sts genome.fa
   ```

2. **Use multiple threads appropriately**:
   ```bash
   # Check CPU cores
   nproc
   
   # Use all cores for large sequences
   merpcr -T $(nproc) markers.sts large_genome.fa
   ```

3. **Process in chunks**:
   ```bash
   # Split large files
   split -l 1000000 huge_genome.fa genome_chunk_
   
   # Process each chunk
   for chunk in genome_chunk_*; do
       merpcr markers.sts $chunk > ${chunk}_results.txt
   done
   ```

4. **Optimize I/O**:
   ```bash
   # Use fast storage for temporary files
   export TMPDIR=/fast/ssd/tmp
   
   # Redirect output efficiently  
   merpcr markers.sts genome.fa > /fast/ssd/results.txt
   ```

### Benchmarking

Use the built-in performance tests:

```bash
# Run performance benchmarks
make test-performance

# Or manually
SKIP_PERFORMANCE_TESTS= pytest tests/test_performance.py -v
```

This will test:
- Large sequence processing speed
- Memory usage patterns
- Threading efficiency
- Scalability characteristics

The results help you choose optimal parameters for your specific system and data.