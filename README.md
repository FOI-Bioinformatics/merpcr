<div align="center">
  <img src="assets/images/merPCR_logo.png" alt="merPCR: modern electronic PCR implementation" width="600">

  [![CI](https://github.com/FOI-Bioinformatics/merpcr/actions/workflows/ci.yml/badge.svg)](https://github.com/FOI-Bioinformatics/merpcr/actions/workflows/ci.yml)
  [![Python 3.13+](https://img.shields.io/badge/python-3.13+-blue.svg)](https://www.python.org/downloads/)
  [![License: GPL v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)
</div>

merPCR is a Python reimplementation of me-PCR (Multithreaded Electronic PCR; Schuler 1997, Murphy 2008). It locates sequence-tagged sites (STS) in nucleotide sequences by searching for a primer pair at the expected product size, and it reproduces the output of me-PCR 1.0.6 for the inputs that program accepts. The command-line interface accepts both the original `M=50` option syntax and conventional `-M 50` options, and the search engine is also available as a Python API.

## Installation

```bash
pip install merpcr            # from PyPI
pip install -e .              # from a clone of this repository
```

Python 3.13 or later and NumPy are required. A Cython extension is built during installation when a C compiler is available; without it the pure-Python implementation is used, with identical results.

## Usage

```bash
merpcr primers.sts genome.fa                    # defaults, results on stdout
merpcr primers.sts genome.fa -M 100 -N 1 -T 4   # margin, mismatches, threads
merpcr primers.sts genome.fa M=100 N=1 T=4      # me-PCR option syntax
merpcr primers.sts genome.fa -O hits.txt        # write results to a file
```

| Option | me-PCR form | Default | Meaning |
|--------|-------------|---------|---------|
| `-M` | `M=` | 50 | Allowed deviation from the expected product size (bp) |
| `-N` | `N=` | 0 | Mismatches allowed per primer |
| `-W` | `W=` | 11 | Hash word size (3 to 16) |
| `-X` | `X=` | 1 | Number of 3' bases in which no mismatch is allowed |
| `-I` | `I=` | 0 | Interpret IUPAC ambiguity codes in primers (0 or 1) |
| `-Z` | `Z=` | 240 | Product size used when the STS file gives 0 |
| `-T` | `T=` | 1 | Number of worker processes |
| `-S` | `S=` | 1022 | Maximum STS line length |
| `-O` | `O=` | stdout | Output file |
| `-Q` | `Q=` | 1 | 0 prints progress messages, 1 is quiet |

The STS file is tab-delimited with the fields identifier, forward primer, reverse primer, expected product size and an optional annotation. The output has one line per hit:

```
L78833	75823..76023	AFM248yg9	(D17S932)  Chr.17, 63.7 cM	(-)
```

The fields are sequence label, 1-based start and end of the product, STS identifier, the annotation when present, and the strand. The [user guide](docs/USER_GUIDE.md) describes the file formats, the parameters and the search algorithm in detail.

### Python API

```python
from merpcr import MerPCR

engine = MerPCR(margin=50, mismatches=1, threads=4)
engine.load_sts_file("primers.sts")
records = engine.load_fasta_file("genome.fa")
hits = engine.find_hits(records)          # list of SearchHit objects
engine.search(records, output_file="hits.txt")  # or write the text output
```

See the [API reference](docs/API.md).

## Compatibility with me-PCR

Output equivalence with the me-PCR 1.0.6 binary has been established in three ways: on the me-PCR test data, on 42 bacterial genomes with a common primer set, and by differential testing on randomly generated genomes with planted primer pairs under random parameter settings (`scripts/compare_with_mepcr.py`). The [verification report](docs/VERIFICATION.md) lists the algorithmic details that had to be matched, including the choice of hash word, the handling of size ranges, reverse-complement orientation, the one-directional interpretation of IUPAC codes and the order in which hits are reported.

Two deliberate differences remain. Size ranges such as `150-250` and the value `-` are accepted as described in the me-PCR source, although the 1.0.6 binary rejects them. STS lines that follow a comment line are reported correctly, whereas me-PCR 1.0.6 reports the wrong line.

## Performance

With the Cython extension merPCR processes bacterial genomes of 2 to 3 Mb in about 0.3 s per genome with a single worker, roughly 2.6 times faster than the pure-Python path. Multiple workers split each sequence into overlapping chunks, as me-PCR does. Measurements and build notes are given in [docs/PERFORMANCE.md](docs/PERFORMANCE.md).

## Development

```bash
make test        # test suite
make lint        # flake8, black and isort checks
make format      # apply black and isort
python scripts/compare_with_mepcr.py --mepcr /path/to/me-PCR -n 50
```

Contribution guidelines, the repository layout and the continuous-integration setup are described in [docs/DEVELOPMENT.md](docs/DEVELOPMENT.md).

## References

- Schuler GD (1997) Sequence mapping by electronic PCR. Genome Research 7:541-550. doi:10.1101/gr.7.5.541
- Murphy K, Raj T, Winters RS, White PS (2004) me-PCR: a refined ultrafast algorithm for identifying sequence-defined genomic elements. Bioinformatics 20:588-590. doi:10.1093/bioinformatics/btg449

The original e-PCR algorithm was written by Gregory D. Schuler (NCBI) and extended with multithreading by Kevin Murphy (Children's Hospital of Philadelphia).

## License

GNU General Public License v3.0. See [LICENSE](LICENSE).
