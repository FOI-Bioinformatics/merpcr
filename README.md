<div align="center">
  <img src="assets/images/merPCR_logo.png" alt="merPCR Logo" width="400">

  # merPCR - Modern Electronic PCR Implementation

  [![Tests](https://github.com/FOI-Bioinformatics/merpcr/workflows/tests/badge.svg)](https://github.com/FOI-Bioinformatics/merpcr/actions)
  [![codecov](https://codecov.io/gh/FOI-Bioinformatics/merpcr/branch/main/graph/badge.svg)](https://codecov.io/gh/FOI-Bioinformatics/merpcr)
  [![Python 3.8+](https://img.shields.io/badge/python-3.8+-blue.svg)](https://www.python.org/downloads/)
  [![License: GPL v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)

  **100% Compatible Python reimplementation of me-PCR**

  [📖 Documentation](docs/README.md) | [🚀 Quick Start](#quick-start) | [✅ Verification](docs/VERIFICATION.md) | [🔄 Migration Guide](docs/MIGRATION.md)
</div>

---

## Overview

merPCR locates Sequence-Tagged Sites (STS) within genomic sequences using computational PCR. It's a **drop-in replacement** for the original me-PCR (Multithreaded Electronic PCR), producing **identical results** while offering modern Python architecture, better error messages, and comprehensive documentation.

**Key Highlights:**
- ✅ **100% Compatible** - Verified byte-for-byte identical output to me-PCR v1.0.6
- 🚀 **Drop-In Replacement** - Same command-line interface, no changes needed
- 🐍 **Python API** - Use programmatically in your Python workflows
- 📚 **Well Documented** - Extensive guides, examples, and API reference
- 🧪 **Thoroughly Tested** - 277 tests with 94% coverage for engine.py on real genomic data

## Compatibility & Validation

merPCR has been **extensively validated** against me-PCR:
- **Real Genomic Data:** Tested on 42 Francisellaceae genomes (90MB)
- **Compatibility Tests:** 15/15 passed (100%)
- **Output Identity:** MD5 checksums match exactly
- **Critical Fixes:** Three algorithmic differences identified and fixed

Three critical compatibility fixes were implemented:
1. **Hash computation** - Backward vs forward search
2. **PCR range margins** - Proper margin adjustment for size ranges
3. **Forward strand matching** - Reverse complement handling for primer2

📄 **Full verification details:** [docs/VERIFICATION.md](docs/VERIFICATION.md)

---

## Quick Start

### Installation

```bash
# From PyPI (when available)
pip install merpcr

# From source
git clone https://github.com/FOI-Bioinformatics/merpcr.git
cd merpcr
pip install -e .
```

### Basic Usage

```bash
# Use merPCR exactly like me-PCR
merpcr primers.sts genome.fa

# With parameters (both formats supported)
merpcr primers.sts genome.fa -M 50 -N 1 -O results.txt
merpcr primers.sts genome.fa M=50 N=1 O=results.txt  # Legacy format works too!

# Multiple parameters, mixed formats
merpcr primers.sts genome.fa -M 50 -N 1 -T 4 --debug
```

### Python API

```python
from merpcr import MerPCR

# Initialize with parameters
engine = MerPCR(wordsize=11, margin=50, mismatches=1, threads=4)

# Load data and search
engine.load_sts_file("primers.sts")
records = engine.load_fasta_file("genome.fa")
hit_count = engine.search(records, "results.txt")

print(f"Found {hit_count} hits")
```

📖 **More examples:** [docs/EXAMPLES.md](docs/EXAMPLES.md)

---

## Key Features

### Computational Features
- **Multithreaded Processing** - Automatic thread scaling for large files
- **Hash-Based Search** - O(1) STS lookup with 2-bit encoding
- **IUPAC Support** - Optional ambiguity code handling
- **Flexible Parameters** - Configurable margins, mismatches, and word sizes
- **3' Protection** - Prevents mismatches in primer 3' regions

### Software Features
- **Modern Architecture** - Type-safe Python with comprehensive error handling
- **Better Diagnostics** - Clear error messages with context
- **Debug Mode** - Detailed logging for troubleshooting
- **Extensive Testing** - 277 tests covering edge cases and real data
- **CI/CD Pipeline** - Automated testing on multiple platforms

---

## Documentation

### For Users
- 📖 **[User Guide](docs/USER_GUIDE.md)** - Comprehensive usage guide
- 💡 **[Examples](docs/EXAMPLES.md)** - Real-world use cases and tutorials
- 🔄 **[Migration Guide](docs/MIGRATION.md)** - Migrating from me-PCR

### For Developers
- 🔧 **[API Reference](docs/API.md)** - Complete API documentation
- ✅ **[Verification Report](docs/VERIFICATION.md)** - Compatibility testing details
- ⚙️ **[CI/CD Documentation](docs/CI_CD.md)** - GitHub Actions workflows

### Quick Reference

**Input Format (STS file):**
```
STS_ID	Forward_Primer	Reverse_Primer	PCR_Size	[Optional_Alias]
```

**Output Format:**
```
Sequence_ID	pos1..pos2	STS_ID	Alias	(+/-)
```

**Common Parameters:**
- `-M, --margin` - Search margin in bp (default: 50)
- `-N, --mismatches` - Allowed mismatches (default: 0)
- `-W, --wordsize` - Hash word size (default: 11)
- `-T, --threads` - Number of threads (default: 1)
- `-O, --output` - Output file (default: stdout)
- `--debug` - Enable debug logging

📖 **Full parameter list:** [docs/USER_GUIDE.md#parameters](docs/USER_GUIDE.md#parameters)

---

## Testing

```bash
# Run all tests
make test

# Run specific test categories
make test-unit          # Unit tests only
make test-integration   # Integration tests
make test-performance   # Performance benchmarks

# Generate coverage report
make coverage

# Run compatibility tests
python test_compatibility.py
```

**Current Status:**
- 277 tests (all passing)
- 94% code coverage for engine.py (critical component)
- Real genomic data validation on 42 genomes

---

## Performance

merPCR performs comparably to me-PCR for most use cases:

| Dataset Size | me-PCR | merPCR | Notes |
|-------------|--------|---------|-------|
| Small (<1MB) | ~0.5s | ~0.6s | Negligible difference |
| Medium (10-100MB) | ~10s | ~12s | ~20% Python overhead |
| Large (>100MB) | Scales linearly | Scales linearly | Multithreading mitigates overhead |

**Real Genomic Data (42 genomes, 90MB total):**
- me-PCR: 0.02s average per genome
- merPCR: 0.02s average per genome
- **No practical difference**

🚀 **For high-throughput:** Use multithreading (`-T` parameter) to maximize performance.

---

## Project Structure

```
merpcr/
├── src/merpcr/          # Main package
│   ├── cli.py           # Command-line interface
│   ├── core/            # Core functionality
│   │   ├── engine.py    # Search engine
│   │   ├── models.py    # Data models
│   │   └── utils.py     # Utilities
│   └── io/              # Input/output
│       ├── fasta.py     # FASTA handling
│       └── sts.py       # STS handling
├── tests/               # Comprehensive test suite (277 tests)
├── docs/                # Documentation
│   ├── USER_GUIDE.md    # Usage documentation
│   ├── API.md           # API reference
│   ├── EXAMPLES.md      # Practical examples
│   ├── VERIFICATION.md  # Compatibility verification
│   ├── MIGRATION.md     # Migration from me-PCR
│   └── CI_CD.md         # CI/CD documentation
├── pyproject.toml       # Modern Python packaging
└── Makefile             # Development commands
```

---

## Contributing

Contributions are welcome! Please:
1. Fork the repository
2. Create a feature branch
3. Add tests for new features
4. Ensure all tests pass (`make test`)
5. Format code (`make format`)
6. Submit a pull request

**Development Setup:**
```bash
git clone https://github.com/FOI-Bioinformatics/merpcr.git
cd merpcr
make dev-install  # Install with development dependencies
make test         # Verify installation
```

---

## License

This project is licensed under the GNU General Public License v3.0 - see the [LICENSE](LICENSE) file for details.

---

## Acknowledgments

merPCR builds upon the pioneering work of:
- **Gregory D. Schuler** (NCBI) - Original e-PCR algorithm development
- **Kevin Murphy** (Children's Hospital of Philadelphia) - me-PCR multithreading enhancements

## References

1. Schuler, G.D. (1997) "Sequence mapping by electronic PCR." *Genome Research* **7**: 541-550. doi:10.1101/gr.7.5.541

2. Altschul, S.F., Gish, W., Miller, W., Myers, E.W., and Lipman, D.J. (1990) "Basic local alignment search tool." *Journal of Molecular Biology* **215**: 403-410. doi:10.1016/S0022-2836(05)80360-2

---

<div align="center">
  <p><strong>merPCR</strong> - Production-ready electronic PCR in Python</p>
  <p>Maintained with ❤️ by the merPCR team</p>
</div>
