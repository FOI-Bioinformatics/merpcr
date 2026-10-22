# Performance

## Implementation

The scan over the sequence and the primer comparison are implemented twice: in `src/merpcr/core/engine.py` (pure Python, NumPy lookup tables) and in `src/merpcr/core/engine_fast.pyx` (Cython). The Cython module is compiled during `pip install` when a C compiler is available and is used automatically; otherwise the engine falls back to the Python implementation. Both paths are tested against each other and against me-PCR and produce identical output.

The Cython path avoids interpreter overhead in the rolling-hash loop and compares primers through byte memoryviews. Compiler flags are `-O3`, `-march=native` and `-ffast-math`; bounds checking and negative indexing are disabled.

To rebuild the extension after editing the `.pyx` file:

```bash
python setup.py build_ext --inplace
```

To check which path is active:

```bash
merpcr primers.sts genome.fa --debug 2>&1 | grep -i cython
```

## Measurements

Single worker, default parameters, Apple M-series laptop, wall-clock time for the complete run.

| Genome | Size | Pure Python | Cython | Ratio |
|--------|------|-------------|--------|-------|
| Francisella tularensis | 1.8 Mb | 0.51 s | 0.24 s | 2.1 |
| Fangia hongkongensis | 2.8 Mb | 0.80 s | 0.27 s | 2.9 |
| Cysteiniphilum litorale | 3.1 Mb | 0.86 s | 0.30 s | 2.8 |

The ratio grows with sequence length because process start-up and file parsing are a fixed cost. Run time is dominated by the scan and scales linearly with sequence length; the number of markers has little influence until hash collisions become frequent.

`benchmark_cython.py` reproduces the comparison by running merPCR with and without the extension on a set of genomes; its data paths must be edited for a different machine.

## Multiple workers

`-T n` splits each sequence into `n` overlapping chunks that are searched in separate processes. The overlap is `max product size + margin - 1` bases, and duplicated hits in the overlap are removed. Files below 100 kB are processed with a single worker regardless of `-T`. Each worker receives a copy of the hash table, so memory grows with the number of workers; for bacterial genomes a single worker is usually adequate.

## Memory

The hash table is a Python dictionary keyed by the hash values that occur in the primer set, so memory does not grow with `4^W` as it does in me-PCR. Memory use is dominated by the sequence itself: roughly two to three times the FASTA file size during a run.
