# Examples and Tutorials

This document provides practical examples and step-by-step tutorials for using merPCR in various scenarios.

## Table of Contents

1. [Basic Examples](#basic-examples)
2. [Tutorial 1: First-Time Setup](#tutorial-1-first-time-setup)
3. [Tutorial 2: Human Genome Analysis](#tutorial-2-human-genome-analysis)
4. [Tutorial 3: Cross-Species Comparison](#tutorial-3-cross-species-comparison)
5. [Tutorial 4: High-Throughput Analysis](#tutorial-4-high-throughput-analysis)
6. [Advanced Examples](#advanced-examples)
7. [Python API Examples](#python-api-examples)
8. [Real-World Case Studies](#real-world-case-studies)

## Basic Examples

### Example 1: Simple Search

**Scenario**: Search for a few known STSs in a single chromosome.

**Data**: 
- `markers.sts`: 5 STS markers
- `chr21.fa`: Human chromosome 21

**Command**:
```bash
merpcr markers.sts chr21.fa
```

**Expected Output**:
```
chr21	14400000..14400150	D21S11	(+)
chr21	25600000..25600200	D21S1435	(-)
```

### Example 2: Permissive Search

**Scenario**: Search with relaxed parameters to find more matches.

**Command**:
```bash
merpcr -M 100 -N 2 -X 0 markers.sts chr21.fa
```

**Use Case**: When working with:
- Divergent species
- Low-quality sequence data
- Degenerate primers

### Example 3: Performance-Optimized

**Scenario**: Fast search of large genome with multiple cores.

**Command**:
```bash
merpcr -W 12 -T 8 -Q 1 -O results.txt markers.sts whole_genome.fa
```

**Benefits**:
- Faster execution with 8 threads
- Larger word size for speed
- Quiet mode for clean output

## Tutorial 1: First-Time Setup

This tutorial walks you through your first merPCR analysis from start to finish.

### Step 1: Install merPCR

```bash
# Clone the repository
git clone https://github.com/yourusername/merpcr.git
cd merpcr

# Install in development mode
pip install -e .

# Verify installation
merpcr --version
```

### Step 2: Prepare Test Data

Create a simple STS file (`test_markers.sts`):
```
# Simple test markers
TEST001	ATCGATCGATCGATCG	CGATCGATCGATCGAT	150	Test marker 1
TEST002	GGCCTTAAGGCCTTAA	TTAAGGCCTTAAGGCC	200	Test marker 2
TEST003	AAAACCCCGGGGTTTT	AAAACCCCGGGGTTTT	100-300	Variable size marker
```

Create a test FASTA file (`test_sequence.fa`):
```
>test_chromosome Test sequence for merPCR
ATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG
NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN
ATCGATCGATCGATCGCGATCGATCGATCGATTTAAGGCCTTAAGGCCAAAACCCCGGGG
TTTTAAAACCCCGGGGTTTTCGATCGATCGATCGATCG
```

### Step 3: Run Your First Search

```bash
# Basic search
merpcr test_markers.sts test_sequence.fa

# With verbose output
merpcr -Q 0 test_markers.sts test_sequence.fa
```

### Step 4: Understand the Output

Expected results might look like:
```
test_chromosome	1..150	TEST001	(+)
test_chromosome	150..350	TEST003	(+)
```

### Step 5: Experiment with Parameters

```bash
# Allow more mismatches
merpcr -N 1 test_markers.sts test_sequence.fa

# Increase margin
merpcr -M 100 test_markers.sts test_sequence.fa

# Enable debug mode
merpcr --debug -Q 0 test_markers.sts test_sequence.fa
```

## Tutorial 2: Human Genome Analysis

This tutorial demonstrates a realistic analysis using human genomic data.

### Scenario

You have a set of STSs from the Human Genome Database and want to map them to chromosome 22.

### Step 1: Download Data

```bash
# Create data directory
mkdir -p data/tutorial2

# Download chromosome 22 (example URLs - replace with actual sources)
wget -O data/tutorial2/chr22.fa.gz \
  "http://hgdownload.cse.ucsc.edu/goldenPath/hg38/chromosomes/chr22.fa.gz"
gunzip data/tutorial2/chr22.fa.gz

# Prepare STS file (example data)
cat > data/tutorial2/chr22_markers.sts << 'EOF'
D22S303	TGAGTCAGATGTTTGATTTTG	ATGCCACATCAACTTATACTG	126	22q11.1
D22S420	GAATGAACAGAGATGATGCCT	CACACACACACACACCACAC	142	22q11.21
WI-18708	GCTAAAAATACACGGATGG	TGCAAGACTGCGTCTC	193	22q12.1
D22S685	CTCTCTCTCTCTCTCTCTGA	GAGAGAGAGAGAGAGAGAGC	160-180	22q12.3
SHGC-145891	AAGGAGGGAAAAGTAGAAGC	TCTGTTCTGTTTCTGTGTCT	201	22q13.1
EOF
```

### Step 2: Basic Analysis

```bash
# Standard search
merpcr -M 50 -N 1 -W 11 -O data/tutorial2/basic_results.txt \
       data/tutorial2/chr22_markers.sts data/tutorial2/chr22.fa

# View results
head -10 data/tutorial2/basic_results.txt
wc -l data/tutorial2/basic_results.txt
```

### Step 3: Quality Control

```bash
# Strict search for high-confidence hits
merpcr -M 25 -N 0 -X 2 -O data/tutorial2/strict_results.txt \
       data/tutorial2/chr22_markers.sts data/tutorial2/chr22.fa

# Permissive search to find potential additional sites
merpcr -M 100 -N 2 -X 0 -O data/tutorial2/permissive_results.txt \
       data/tutorial2/chr22_markers.sts data/tutorial2/chr22.fa

# Compare results
echo "Strict hits: $(wc -l < data/tutorial2/strict_results.txt)"
echo "Basic hits: $(wc -l < data/tutorial2/basic_results.txt)"
echo "Permissive hits: $(wc -l < data/tutorial2/permissive_results.txt)"
```

### Step 4: Analysis and Interpretation

```bash
# Extract unique STSs found
cut -f3 data/tutorial2/basic_results.txt | sort | uniq -c | sort -nr

# Check for multiple hits per STS (potential duplications)
cut -f3 data/tutorial2/basic_results.txt | sort | uniq -d

# Analyze size distribution
awk '{print $3, ($2 - $1 + 1)}' data/tutorial2/basic_results.txt | \
    awk '{split($1,pos,".."); print $2, pos[2]-pos[1]+1}' | \
    sort -k2 -n
```

## Tutorial 3: Cross-Species Comparison

This tutorial shows how to use merPCR for comparative genomics.

### Scenario

You have human STS markers and want to find homologous regions in mouse genome.

### Step 1: Prepare Data

```bash
mkdir -p data/tutorial3

# Human STSs (example)
cat > data/tutorial3/human_markers.sts << 'EOF'
HUM_001	ATCGATCGATCGATCG	CGATCGATCGATCGAT	150	Human specific
HUM_002	GGCCTTAAGGCCTTAA	TTAAGGCCTTAAGGCC	200	Conserved region
HUM_003	AAAACCCCGGGGTTTT	TTTTGGGGCCCCAAAA	180	Variable region
CONS_001	GAATTCGAATTCGAAT	AATTCGAATTCGAATT	160	Highly conserved
EOF

# Mouse sequence (example - replace with real data)
cat > data/tutorial3/mouse_chr1.fa << 'EOF'
>mouse_chromosome1 Mus musculus chromosome 1
ATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG
NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN
ATCGATCGATCGATCGCGATCGATCGATCGATTTAAGGCCTTAAGGCCAAAACCCCGGGG
TTTTAAAACCCCGGGGTTTTCGATCGATCGATCGATCGAATTCGAATTCGAATTCGAAT
TCGAATTCGAATTCGAATTCGAATTCGAATTCGCGATCGATCGATCGAT
EOF
```

### Step 2: Cross-Species Search

```bash
# Permissive search allowing for species differences
merpcr -M 150 -N 3 -X 0 -I 1 -W 9 -Q 0 \
       -O data/tutorial3/cross_species_results.txt \
       data/tutorial3/human_markers.sts data/tutorial3/mouse_chr1.fa
```

### Step 3: Analyze Conservation

```bash
# Check which markers were found
echo "Markers found in mouse:"
cut -f3 data/tutorial3/cross_species_results.txt | sort | uniq

# Compare hit sizes with expected (indicates conservation level)
awk '{
    split($2, pos, "\\.\\."); 
    hit_size = pos[2] - pos[1] + 1;
    print $3, hit_size
}' data/tutorial3/cross_species_results.txt > data/tutorial3/hit_sizes.txt

# Join with original expected sizes for comparison
# (This requires more complex shell scripting or manual analysis)
```

### Step 4: Visualization

```bash
# Create simple report
cat > data/tutorial3/conservation_report.txt << 'EOF'
Cross-Species STS Analysis Report
=================================

Search Parameters:
- Margin: 150 bp
- Mismatches allowed: 3
- 3' protection: disabled
- IUPAC mode: enabled

Results:
EOF

echo "Total hits found: $(wc -l < data/tutorial3/cross_species_results.txt)" >> \
     data/tutorial3/conservation_report.txt

echo -e "\nHit details:" >> data/tutorial3/conservation_report.txt
cat data/tutorial3/cross_species_results.txt >> data/tutorial3/conservation_report.txt
```

## Tutorial 4: High-Throughput Analysis

This tutorial covers batch processing of multiple files and large-scale analysis.

### Scenario

You need to search 1000 STSs across all human chromosomes.

### Step 1: Organize Data

```bash
mkdir -p data/tutorial4/{chromosomes,results}

# Example structure:
# data/tutorial4/
# ├── chromosomes/
# │   ├── chr1.fa
# │   ├── chr2.fa
# │   └── ...
# ├── markers_set1.sts
# ├── markers_set2.sts
# └── results/
```

### Step 2: Batch Processing Script

Create `batch_analysis.sh`:

```bash
#!/bin/bash

# Configuration
STS_DIR="data/tutorial4"
CHR_DIR="data/tutorial4/chromosomes"  
RES_DIR="data/tutorial4/results"
THREADS=8

# Parameters
MARGIN=75
MISMATCHES=1
WORDSIZE=11

echo "Starting batch STS analysis..."
echo "Using ${THREADS} threads per job"

# Process each STS file against all chromosomes
for sts_file in "${STS_DIR}"/markers_set*.sts; do
    sts_name=$(basename "$sts_file" .sts)
    echo "Processing $sts_name..."
    
    # Process each chromosome
    for chr_file in "${CHR_DIR}"/chr*.fa; do
        chr_name=$(basename "$chr_file" .fa)
        output_file="${RES_DIR}/${sts_name}_${chr_name}.txt"
        
        echo "  Searching $chr_name..."
        
        # Run merPCR with timeout (in case of issues)
        timeout 1800 merpcr \
            -M $MARGIN \
            -N $MISMATCHES \
            -W $WORDSIZE \
            -T $THREADS \
            -Q 1 \
            -O "$output_file" \
            "$sts_file" "$chr_file"
        
        # Check for success
        if [ $? -eq 0 ]; then
            hits=$(wc -l < "$output_file")
            echo "    Found $hits hits"
        else
            echo "    ERROR: Processing failed"
            rm -f "$output_file"
        fi
    done
done

echo "Batch processing complete."
```

### Step 3: Parallel Processing

For even faster processing, use GNU parallel:

```bash
# Install parallel if not available
# sudo apt-get install parallel  # Ubuntu/Debian
# brew install parallel          # macOS

# Create parallel processing script
create_parallel_jobs() {
    for sts_file in data/tutorial4/markers_set*.sts; do
        for chr_file in data/tutorial4/chromosomes/chr*.fa; do
            sts_name=$(basename "$sts_file" .sts)
            chr_name=$(basename "$chr_file" .fa)
            output_file="data/tutorial4/results/${sts_name}_${chr_name}.txt"
            
            echo "merpcr -M 75 -N 1 -W 11 -T 2 -Q 1 -O $output_file $sts_file $chr_file"
        done
    done > parallel_jobs.txt
}

# Generate job list
create_parallel_jobs

# Run in parallel (adjust -j based on your system)
parallel -j 4 < parallel_jobs.txt
```

### Step 4: Results Aggregation

```bash
# Combine all results
cat data/tutorial4/results/*.txt > data/tutorial4/all_results.txt

# Generate summary statistics
echo "=== STS Analysis Summary ===" > data/tutorial4/summary.txt
echo "Total hits found: $(wc -l < data/tutorial4/all_results.txt)" >> data/tutorial4/summary.txt
echo "" >> data/tutorial4/summary.txt

echo "Hits per chromosome:" >> data/tutorial4/summary.txt
awk '{print $1}' data/tutorial4/all_results.txt | sort | uniq -c | sort -nr >> data/tutorial4/summary.txt
echo "" >> data/tutorial4/summary.txt

echo "Hits per STS:" >> data/tutorial4/summary.txt  
awk '{print $3}' data/tutorial4/all_results.txt | sort | uniq -c | sort -nr | head -20 >> data/tutorial4/summary.txt

# Find STSs with multiple hits (potential duplications)
echo "" >> data/tutorial4/summary.txt
echo "STSs with multiple hits:" >> data/tutorial4/summary.txt
awk '{print $3}' data/tutorial4/all_results.txt | sort | uniq -c | awk '$1 > 1 {print $0}' >> data/tutorial4/summary.txt
```

## Advanced Examples

### Example 1: Parameter Optimization

Find optimal parameters for your specific dataset:

```bash
# Test different parameter combinations
test_parameters() {
    local sts_file=$1
    local fasta_file=$2
    local output_prefix=$3
    
    for margin in 25 50 100; do
        for mismatches in 0 1 2; do
            for wordsize in 9 11 13; do
                output_file="${output_prefix}_M${margin}_N${mismatches}_W${wordsize}.txt"
                
                echo "Testing M=$margin N=$mismatches W=$wordsize"
                time merpcr -M $margin -N $mismatches -W $wordsize \
                           -O "$output_file" "$sts_file" "$fasta_file"
                
                hits=$(wc -l < "$output_file")
                echo "  Hits found: $hits"
            done
        done
    done
}

# Usage
test_parameters markers.sts genome.fa results/param_test
```

### Example 2: Quality Filtering

Filter results based on various criteria:

```bash
# Filter by PCR product size (only 100-300 bp products)
awk '{
    split($2, pos, "\\.\\."); 
    size = pos[2] - pos[1] + 1;
    if (size >= 100 && size <= 300) print $0
}' results.txt > filtered_results.txt

# Filter by chromosome (exclude random/unknown sequences)
grep -E '^(chr[1-9]|chr1[0-9]|chr2[0-2]|chrX|chrY)' results.txt > main_chr_results.txt

# Remove duplicate hits (same STS, same position)
sort results.txt | uniq > unique_results.txt

# Filter by strand (only forward strand hits)
grep '(+)' results.txt > forward_only_results.txt
```

### Example 3: Custom Output Formatting

Create custom output formats:

```bash
# Convert to BED format
awk 'BEGIN{OFS="\t"} {
    split($2, pos, "\\.\\."); 
    print $1, pos[1]-1, pos[2], $3, ".", ($4 == "(+)" ? "+" : "-")
}' results.txt > results.bed

# Create GFF3 format
awk 'BEGIN{OFS="\t"; print "##gff-version 3"} {
    split($2, pos, "\\.\\.");
    print $1, "merPCR", "STS", pos[1], pos[2], ".", ($4 == "(+)" ? "+" : "-"), ".", "ID="$3
}' results.txt > results.gff3

# Create summary table
awk '{
    split($2, pos, "\\.\\."); 
    size = pos[2] - pos[1] + 1;
    print $1, pos[1], pos[2], $3, size, $4
}' results.txt | column -t > formatted_results.txt
```

## Python API Examples

### Example 1: Basic Python Usage

```python
#!/usr/bin/env python3
"""
Basic merPCR usage example
"""

from merpcr import MerPCR
import sys

def basic_search():
    # Initialize with custom parameters
    mer_pcr = MerPCR(
        wordsize=11,
        margin=50,
        mismatches=1,
        threads=4,
        quiet=0  # Verbose output
    )
    
    # Load data
    print("Loading STS file...")
    if not mer_pcr.load_sts_file("markers.sts"):
        print("ERROR: Failed to load STS file")
        return False
    
    print("Loading FASTA file...")
    sequences = mer_pcr.load_fasta_file("genome.fa")
    if not sequences:
        print("ERROR: Failed to load FASTA file")
        return False
    
    print(f"Loaded {len(sequences)} sequences")
    
    # Perform search
    print("Searching for STS matches...")
    hit_count = mer_pcr.search(sequences, "python_results.txt")
    
    print(f"Found {hit_count} total hits")
    return True

if __name__ == "__main__":
    basic_search()
```

### Example 2: Parameter Comparison

```python
#!/usr/bin/env python3
"""
Compare different parameter settings
"""

from merpcr import MerPCR
import time

def compare_parameters():
    # Test different parameter combinations
    test_configs = [
        {"wordsize": 9, "margin": 25, "mismatches": 0, "name": "conservative"},
        {"wordsize": 11, "margin": 50, "mismatches": 1, "name": "standard"},  
        {"wordsize": 13, "margin": 100, "mismatches": 2, "name": "permissive"}
    ]
    
    sequences = None  # Load once, reuse
    results = {}
    
    for config in test_configs:
        print(f"Testing {config['name']} configuration...")
        
        # Create engine with test parameters
        mer_pcr = MerPCR(
            wordsize=config["wordsize"],
            margin=config["margin"], 
            mismatches=config["mismatches"],
            quiet=1
        )
        
        # Load data (only once)
        if sequences is None:
            mer_pcr.load_sts_file("markers.sts")
            sequences = mer_pcr.load_fasta_file("genome.fa")
        else:
            mer_pcr.load_sts_file("markers.sts")
        
        # Time the search
        start_time = time.time()
        hit_count = mer_pcr.search(sequences)
        elapsed_time = time.time() - start_time
        
        results[config["name"]] = {
            "hits": hit_count,
            "time": elapsed_time,
            "config": config
        }
        
        print(f"  Hits: {hit_count}, Time: {elapsed_time:.2f}s")
    
    # Print summary
    print("\n=== Parameter Comparison Summary ===")
    for name, result in results.items():
        print(f"{name:12s}: {result['hits']:4d} hits in {result['time']:5.2f}s")

if __name__ == "__main__":
    compare_parameters()
```

### Example 3: Batch Processing

```python
#!/usr/bin/env python3
"""
Batch process multiple files
"""

from merpcr import MerPCR
from pathlib import Path
import json

def batch_process():
    # Configuration
    config = {
        "wordsize": 11,
        "margin": 50,
        "mismatches": 1,
        "threads": 4
    }
    
    # File paths
    sts_files = list(Path("data").glob("*.sts"))
    fasta_files = list(Path("data").glob("*.fa"))
    
    results_summary = {}
    
    # Process each combination
    for sts_file in sts_files:
        for fasta_file in fasta_files:
            print(f"Processing {sts_file.name} vs {fasta_file.name}")
            
            # Create engine
            mer_pcr = MerPCR(**config)
            
            # Load data
            if not mer_pcr.load_sts_file(str(sts_file)):
                print(f"  ERROR: Could not load {sts_file}")
                continue
            
            sequences = mer_pcr.load_fasta_file(str(fasta_file))
            if not sequences:
                print(f"  ERROR: Could not load {fasta_file}")
                continue
            
            # Search
            output_file = f"results/{sts_file.stem}_{fasta_file.stem}.txt"
            hit_count = mer_pcr.search(sequences, output_file)
            
            # Store results
            key = f"{sts_file.stem}_{fasta_file.stem}"
            results_summary[key] = {
                "sts_file": str(sts_file),
                "fasta_file": str(fasta_file),
                "hits": hit_count,
                "output": output_file
            }
            
            print(f"  Found {hit_count} hits -> {output_file}")
    
    # Save summary
    with open("batch_results_summary.json", "w") as f:
        json.dump(results_summary, f, indent=2)
    
    print(f"\nBatch processing complete. Processed {len(results_summary)} combinations.")

if __name__ == "__main__":
    batch_process()
```

## Real-World Case Studies

### Case Study 1: Genome Assembly Validation

**Problem**: Validate that a new genome assembly contains expected STS markers.

**Approach**:
```bash
# Use strict parameters for high confidence
merpcr -M 25 -N 0 -X 2 -W 11 -O validation_results.txt \
       known_markers.sts new_assembly.fa

# Check coverage
total_markers=$(wc -l < known_markers.sts)
found_markers=$(cut -f3 validation_results.txt | sort | uniq | wc -l)
coverage=$(echo "scale=2; $found_markers * 100 / $total_markers" | bc)

echo "Assembly validation: $found_markers/$total_markers markers found (${coverage}%)"
```

**Interpretation**:
- High coverage (>95%): Assembly likely complete
- Low coverage (<80%): Potential assembly gaps or errors
- Multiple hits per STS: Potential duplications or repeats

### Case Study 2: Comparative Genomics

**Problem**: Identify conserved genomic regions between human and mouse.

**Approach**:
```bash
# Permissive search allowing for species divergence
merpcr -M 200 -N 3 -X 0 -I 1 -W 9 \
       human_conserved_stss.sts mouse_genome.fa > conserved_regions.txt

# Analyze conservation patterns
awk '{
    split($2, pos, "\\.\\."); 
    size = pos[2] - pos[1] + 1;
    print $1, pos[1], pos[2], $3, size
}' conserved_regions.txt | sort -k1,1 -k2,2n > sorted_conserved.txt
```

**Results Analysis**:
- Clustered hits suggest conserved synteny blocks
- Size variations indicate evolutionary insertions/deletions
- Missing STSs highlight species-specific regions

### Case Study 3: Medical Genetics

**Problem**: Map disease-associated STSs to patient genome sequences.

**Approach**:
```bash
# High-sensitivity search for clinical variants
merpcr -M 100 -N 2 -I 1 -W 10 -T 8 \
       disease_associated_stss.sts patient_genome.fa > patient_hits.txt

# Check for missing markers (potential deletions)
comm -23 <(cut -f1 disease_associated_stss.sts | sort) \
         <(cut -f3 patient_hits.txt | sort) > missing_markers.txt

echo "Missing disease markers:"
cat missing_markers.txt
```

**Clinical Interpretation**:
- Present markers: Normal genomic structure
- Missing markers: Potential deletions requiring further investigation
- Multiple hits: Duplications that might affect dosage

### Case Study 4: Phylogenetic Analysis

**Problem**: Use STS presence/absence for evolutionary analysis.

**Approach**:
```python
#!/usr/bin/env python3
"""
Create presence/absence matrix for phylogenetic analysis
"""

from merpcr import MerPCR
import pandas as pd

def phylogenetic_matrix():
    species_genomes = [
        "human_genome.fa",
        "chimp_genome.fa", 
        "gorilla_genome.fa",
        "orangutan_genome.fa"
    ]
    
    sts_file = "phylogenetic_markers.sts"
    
    # Get STS list
    with open(sts_file) as f:
        sts_ids = [line.split('\t')[0] for line in f if not line.startswith('#')]
    
    # Create matrix
    matrix = pd.DataFrame(0, index=sts_ids, columns=species_genomes)
    
    # Process each species
    for genome in species_genomes:
        print(f"Processing {genome}...")
        
        mer_pcr = MerPCR(margin=100, mismatches=2, quiet=1)
        mer_pcr.load_sts_file(sts_file)
        sequences = mer_pcr.load_fasta_file(genome)
        
        # Get results
        mer_pcr.search(sequences, f"temp_{genome}_results.txt")
        
        # Mark present STSs
        with open(f"temp_{genome}_results.txt") as f:
            for line in f:
                sts_id = line.split('\t')[2]
                matrix.loc[sts_id, genome] = 1
    
    # Save matrix
    matrix.to_csv("phylogenetic_matrix.csv")
    print("Phylogenetic matrix saved to phylogenetic_matrix.csv")

if __name__ == "__main__":
    phylogenetic_matrix()
```

This matrix can then be used with phylogenetic software like PHYLIP or RAxML.

## Tips and Best Practices

### Performance Tips

1. **Start small**: Test with a subset of data before full analysis
2. **Monitor resources**: Use `htop` or `top` to monitor CPU and memory usage
3. **Optimize parameters**: Use performance tests to find optimal settings
4. **Parallelize wisely**: Don't use more threads than CPU cores

### Data Quality Tips

1. **Validate input files**: Check for proper formatting before large runs
2. **Use appropriate parameters**: Match parameters to your data quality
3. **Filter results**: Remove low-quality hits based on your criteria
4. **Cross-validate**: Compare results with other tools when possible

### Troubleshooting Tips

1. **Enable debug mode**: Use `--debug` for detailed error information
2. **Check file paths**: Ensure all input files exist and are readable
3. **Verify permissions**: Make sure you can write to output directories
4. **Test incrementally**: Start with simple cases and build complexity

### Documentation Tips

1. **Record parameters**: Always document the parameters used for reproducibility
2. **Version control**: Keep track of merPCR version and input file versions  
3. **Archive results**: Save intermediate results for future reference
4. **Share protocols**: Document your workflows for team members