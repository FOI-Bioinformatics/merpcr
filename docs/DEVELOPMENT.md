# Development

## Layout

```
src/merpcr/
  cli.py              command-line interface, me-PCR option conversion
  core/engine.py      MerPCR class: STS loading, hashing, search, reporting
  core/engine_fast.pyx  Cython implementation of the scan and comparison
  core/models.py      STSRecord, FASTARecord, STSHit, ThreadData
  core/utils.py       reverse complement and hash helpers
  io/fasta.py         FASTA reader
  io/sts.py           STSLoader wrapper around the engine's parser
tests/                pytest suite; tests/data holds the me-PCR test files
scripts/              comparison and benchmark scripts (see below)
docs/                 this documentation
```

`scripts/compare_with_mepcr.py` is the differential test against a me-PCR binary and takes the binary path as an argument. `scripts/test_compatibility.py`, `scripts/validate_real_data.py` and `scripts/benchmark_cython.py` compare against a local me-PCR build and local genome data; their paths are hard-coded and need to be edited before use.

Tool configuration (pytest, coverage, black, isort, mypy) lives in `pyproject.toml`; `tox.ini` holds the tox environments and the flake8 settings.

## Setup

```bash
git clone https://github.com/FOI-Bioinformatics/merpcr.git
cd merpcr
pip install -e ".[dev]"
```

The editable install compiles the Cython extension. After changing `engine_fast.pyx`, run `python setup.py build_ext --inplace`.

## Tests

```bash
make test                       # whole suite
pytest -m "not slow"            # skip long-running tests
pytest tests/test_mepcr_fidelity.py   # behaviors verified against the binary
make coverage                   # HTML report in htmlcov/
```

Markers: `unit`, `integration`, `cli`, `performance`, `slow`. Performance tests can be skipped by setting `SKIP_PERFORMANCE_TESTS=1`.

When changing anything in the search path, run the differential test against a me-PCR binary and include both code paths:

```bash
python scripts/compare_with_mepcr.py --mepcr /path/to/me-PCR -n 50
python scripts/compare_with_mepcr.py --mepcr /path/to/me-PCR -n 20 --pure-python
```

## Style

Black with a line length of 100, isort, flake8. `make format` applies the formatters and `make lint` checks them. Keep comments and documentation plain and factual.

## Continuous integration

Two workflows in `.github/workflows`:

- `ci.yml` runs on pushes and pull requests to `main`. The `test` job installs the package and runs pytest with Python 3.13 on Ubuntu. The `lint` job runs black, isort and flake8 checks. The `build` job builds the source distribution and wheel and checks them with twine.
- `release.yml` runs when a tag `v*` is pushed. It builds the distributions, publishes them to PyPI through trusted publishing (the `pypi` environment must be configured on GitHub), and creates a GitHub release with the files attached.

The differential test against me-PCR is not run in CI because it needs the compiled me-PCR binary.

## Releasing

1. Update the version in `pyproject.toml` and `src/merpcr/__init__.py`, and add a section to `CHANGELOG.md`.
2. Commit, tag `vX.Y.Z` and push the tag.
