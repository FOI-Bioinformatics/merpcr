# Migrating from me-PCR to merPCR

This guide helps users transition from the original me-PCR (C++) implementation to merPCR (Python).

## TL;DR - Quick Migration

merPCR is **100% compatible** with me-PCR v1.0.6:
- ✅ **Same command-line arguments** - no changes needed
- ✅ **Identical output** - byte-for-byte identical results
- ✅ **Drop-in replacement** - works with existing scripts

```bash
# Before (me-PCR)
me-PCR primers.sts genome.fa M=50 N=1 O=results.txt

# After (merPCR) - same command works!
merpcr primers.sts genome.fa M=50 N=1 O=results.txt

# Or use modern format
merpcr primers.sts genome.fa -M 50 -N 1 -O results.txt
```

---

## What's the Same

### 1. Command-Line Interface
merPCR supports **all** me-PCR arguments in **both** formats:

| Parameter | Legacy Format | Modern Format | Description |
|-----------|--------------|---------------|-------------|
| Margin | `M=50` | `-M 50` | Search margin in bp |
| Mismatches | `N=1` | `-N 1` | Allowed mismatches |
| Word size | `W=11` | `-W 11` | Hash word size |
| Threads | `T=4` | `-T 4` | Number of threads |
| Output | `O=file.txt` | `-O file.txt` | Output file |
| PCR size | `Z=240` | `-Z 240` | Default PCR size |
| IUPAC | `I=1` | `-I 1` | IUPAC ambiguity support |
| 3' match | `X=1` | `-X 1` | 3' protection |
| Quiet | `Q=1` | `-Q 1` | Quiet mode |

**Mix and match:** You can use both formats in the same command!
```bash
merpcr primers.sts genome.fa M=50 -N 1 --debug
```

### 2. Input File Formats
- **STS files:** Identical tab-delimited format
- **FASTA files:** Standard FASTA format
- **No conversion needed** - use your existing files

### 3. Output Format
- **Exact same output** - byte-for-byte identical
- **Same field order:** seq_label, positions, STS_ID, alias, direction
- **Same delimiters:** Tab-separated
- **Same coordinates:** 1-based positions

### 4. Behavior
- **Same search algorithm** - identical results
- **Same threading** - automatic thread scaling
- **Same file handling** - stdout and file output

---

## What's Different (Better!)

### 1. Installation is Easier

**me-PCR (C++):**
```bash
# Need to compile from source
make
# Manual path setup
export PATH=$PATH:/path/to/me-PCR/src
```

**merPCR (Python):**
```bash
# One command installation
pip install merpcr

# Or from source
pip install -e .
```

### 2. Error Messages are Clearer

**me-PCR:**
```
Error: bad STS file format
```

**merPCR:**
```
ERROR - Bad STS file format at line 42.
Expected at least 4 fields (ID, primer1, primer2, size).
Found only 3 fields.
```

### 3. Logging is Better

**merPCR adds:**
- Structured logging with timestamps
- Debug mode (`--debug`) for troubleshooting
- Progress information
- Clear warnings

```bash
# Enable detailed logging
merpcr --debug primers.sts genome.fa

# Output:
# 2025-10-13 18:45:32 INFO - Reading STS file: primers.sts
# 2025-10-13 18:45:32 INFO - Loaded 7 STS records in 0.01 seconds
# 2025-10-13 18:45:32 INFO - Processing sequence: chr1 (1862444 bp)
# 2025-10-13 18:45:33 INFO - Total hits found: 20
```

### 4. Python API Available

**New capability** - Use merPCR programmatically:

```python
from merpcr import MerPCR

# Initialize
engine = MerPCR(wordsize=11, margin=50, mismatches=1, threads=4)

# Load data
engine.load_sts_file("primers.sts")
records = engine.load_fasta_file("genome.fa")

# Search
hits = engine.search(records, "results.txt")

# Access results
for record in engine.sts_records:
    print(f"STS: {record.id}, Primers: {record.primer1}, {record.primer2}")
```

### 5. Modern Development Practices

- **Type hints** - Better IDE support
- **Unit tests** - 324 comprehensive tests
- **Documentation** - Extensive API docs
- **CI/CD** - Automated testing
- **Version control** - Active development

---

## Verifying Compatibility

### Quick Test

Run both tools on the same data and compare:

```bash
# Run me-PCR
me-PCR primers.sts genome.fa O=mepcr_output.txt

# Run merPCR
merpcr primers.sts genome.fa -O merpcr_output.txt

# Compare outputs (should be identical)
diff mepcr_output.txt merpcr_output.txt
# No output = files are identical ✓

# Or check MD5 checksums
md5 mepcr_output.txt merpcr_output.txt
# Should show same hash
```

### Comprehensive Validation

merPCR has been validated against me-PCR with:
- **15/15 compatibility tests** passed
- **42 bacterial genomes** (90MB data)
- **Multiple parameter combinations**
- **Real-world STS markers**

See [VERIFICATION.md](VERIFICATION.md) for detailed results.

---

## Performance Comparison

### Speed
- **Small files (<1MB):** Nearly identical speed
- **Medium files (10-100MB):** merPCR ~20% slower (Python overhead)
- **Large files (>100MB):** Similar with multithreading

**Example (42 genomes, 90MB total):**
- me-PCR: 0.02s average per genome
- merPCR: 0.02s average per genome
- **No practical difference for most use cases**

### Memory
- Similar memory usage
- Both scale well with file size
- Both handle GB-scale genomes efficiently

### Threading
- Both scale well with CPU cores
- Automatic thread optimization
- Similar speedup characteristics

---

## Migration Strategies

### Strategy 1: Drop-In Replacement (Recommended)

Simply replace `me-PCR` with `merpcr`:

```bash
# Update your PATH or create alias
alias me-PCR='merpcr'

# Or symbolic link
ln -s $(which merpcr) /usr/local/bin/me-PCR

# Existing scripts work unchanged
```

### Strategy 2: Gradual Migration

1. **Install merPCR** alongside me-PCR
2. **Test on sample data** to verify identical results
3. **Update scripts** one at a time
4. **Remove me-PCR** when comfortable

### Strategy 3: Parallel Running

```bash
# Run both and compare
me-PCR "$@" O=mepcr_out.txt
merpcr "$@" -O merpcr_out.txt
diff mepcr_out.txt merpcr_out.txt || echo "Outputs differ!"
```

---

## Common Migration Scenarios

### Scenario 1: Shell Scripts

**Before:**
```bash
#!/bin/bash
for genome in genomes/*.fa; do
    me-PCR primers.sts "$genome" M=50 N=1 O="results/$(basename $genome .fa).txt"
done
```

**After:**
```bash
#!/bin/bash
for genome in genomes/*.fa; do
    merpcr primers.sts "$genome" -M 50 -N 1 -O "results/$(basename $genome .fa).txt"
done
```

### Scenario 2: Makefile

**Before:**
```makefile
results/%.txt: data/%.fa primers.sts
	me-PCR primers.sts $< M=50 O=$@
```

**After:**
```makefile
results/%.txt: data/%.fa primers.sts
	merpcr primers.sts $< -M 50 -O $@
```

### Scenario 3: Python Subprocess

**Before:**
```python
import subprocess
subprocess.run([
    "me-PCR", "primers.sts", "genome.fa",
    "M=50", "O=results.txt"
])
```

**After (Option 1 - Command line):**
```python
import subprocess
subprocess.run([
    "merpcr", "primers.sts", "genome.fa",
    "-M", "50", "-O", "results.txt"
])
```

**After (Option 2 - Python API):**
```python
from merpcr import MerPCR

engine = MerPCR(margin=50)
engine.load_sts_file("primers.sts")
records = engine.load_fasta_file("genome.fa")
engine.search(records, "results.txt")
```

---

## Troubleshooting Migration Issues

### Issue: "Command not found: merpcr"

**Solution:**
```bash
# Check installation
pip list | grep merpcr

# If not found, install
pip install merpcr

# Or install from source
cd merpcr && pip install -e .

# Verify
which merpcr
merpcr --version
```

### Issue: "Different results between me-PCR and merPCR"

**Solution:**
1. Check versions:
   ```bash
   me-PCR --version  # Should be v1.0.6
   merpcr --version   # Should be v1.0.0+
   ```

2. Verify exact same parameters:
   ```bash
   # Use --debug to see what merPCR is doing
   merpcr --debug primers.sts genome.fa [parameters]
   ```

3. Check file formats:
   ```bash
   # Verify STS file format
   head -5 primers.sts
   # Should be tab-delimited with 4-5 fields
   ```

4. Report issue:
   - Run with `--debug` flag
   - Include both outputs
   - Share command used
   - Report at: https://github.com/FOI-Bioinformatics/merpcr/issues

### Issue: "merPCR is slower"

**Explanation:** Python has ~20% overhead vs C++ for CPU-intensive tasks. This is expected.

**Solutions:**
1. Use multithreading:
   ```bash
   merpcr -T 4 primers.sts genome.fa  # Use 4 threads
   ```

2. For high-throughput pipelines, me-PCR's speed advantage may matter. Consider:
   - Using merPCR for development/testing (better errors, logging)
   - Using me-PCR for production batches (if speed critical)

---

## Benefits of Migrating

### For Users
- ✅ **Easier installation** - `pip install merpcr`
- ✅ **Better error messages** - know what went wrong
- ✅ **Better logging** - understand what's happening
- ✅ **Python API** - integrate into Python workflows
- ✅ **Active development** - bug fixes and improvements

### For Developers
- ✅ **Modern codebase** - Python, type hints, tests
- ✅ **Easy to extend** - add features without C++ complexity
- ✅ **Well documented** - API docs, examples, guides
- ✅ **CI/CD pipeline** - automated testing and releases
- ✅ **Open development** - contribute via pull requests

### For Projects
- ✅ **Maintainable** - easier to fix bugs and add features
- ✅ **Portable** - runs anywhere Python runs
- ✅ **Testable** - comprehensive test suite
- ✅ **Documented** - extensive documentation
- ✅ **Supported** - active development and support

---

## FAQ

### Q: Do I need to rewrite my scripts?
**A:** No! merPCR is a drop-in replacement. Your existing scripts will work without changes.

### Q: Will the output be exactly the same?
**A:** Yes! merPCR produces byte-for-byte identical output to me-PCR v1.0.6.

### Q: Is merPCR slower?
**A:** Slightly (~20%) for single-threaded. With multithreading, the difference is negligible for most use cases.

### Q: Can I use both me-PCR and merPCR?
**A:** Yes! They can coexist. Install both and use whichever is appropriate.

### Q: How do I get help?
**A:**
- Documentation: [docs/USER_GUIDE.md](USER_GUIDE.md)
- Examples: [docs/EXAMPLES.md](EXAMPLES.md)
- Issues: https://github.com/FOI-Bioinformatics/merpcr/issues

### Q: What if I find a bug?
**A:** Please report it! Include:
- Command used
- Input files (or minimal example)
- Expected vs actual output
- Debug log (`--debug` flag)

### Q: Can I contribute?
**A:** Yes! merPCR is open source. See the main [README.md](../README.md) for contribution guidelines.

---

## Getting Help

### Resources
- **Documentation:** [docs/README.md](README.md)
- **User Guide:** [docs/USER_GUIDE.md](USER_GUIDE.md)
- **API Reference:** [docs/API.md](API.md)
- **Examples:** [docs/EXAMPLES.md](EXAMPLES.md)
- **Verification:** [docs/VERIFICATION.md](VERIFICATION.md)

### Support
- **GitHub Issues:** https://github.com/FOI-Bioinformatics/merpcr/issues
- **Email:** (if you want to provide support email)

---

## Summary

✅ **Migration is easy** - merPCR is a drop-in replacement for me-PCR
✅ **No changes needed** - existing scripts work as-is
✅ **Identical results** - verified on real genomic data
✅ **Better experience** - improved errors, logging, and API
✅ **Well supported** - active development and comprehensive docs

**Recommended approach:** Install merPCR, test on your data, then gradually migrate. Both tools can coexist during transition.

---

*Last updated: October 2025*
*Verified against: me-PCR v1.0.6, merPCR v1.0.0*
