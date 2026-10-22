# Verification against me-PCR

This report records how output equivalence between merPCR and me-PCR 1.0.6 was established and which details of the original implementation had to be reproduced.

## Reference

me-PCR 1.0.6 (Murphy, Children's Hospital of Philadelphia), compiled from the distributed C++ source. The relevant code is in `stsmatch.cpp` (STS parsing, hashing, matching, threading and reporting) and `fasta-io.cpp` (sequence parsing).

## Methods

1. The me-PCR test data (`test/test.sts`, `test/test.fa`) and the parameter variations in `scripts/test_compatibility.py`.
2. Forty-two Francisellaceae genome assemblies (about 90 Mb in total) with a shared primer set, compared line by line under several parameter settings (`scripts/validate_real_data.py`; the data are not part of the repository).
3. Differential testing with `scripts/compare_with_mepcr.py`: random genomes of 3 kb to 400 kb with planted primer pairs in both orientations, with and without mismatches, IUPAC codes in primers, N blocks, lower-case input, multi-record files and every form of the annotation column, searched under random `M`, `N`, `W`, `X`, `T`, `Z` and `I` values. Standard output and exit codes are compared byte for byte. At the time of writing, 140 of 140 cases over three seeds (1849 hits) were identical for both the Cython and the pure-Python code path.
4. Regression tests in `tests/test_mepcr_fidelity.py` for every behavior listed below, with expected values taken from the binary.

## Behaviors that had to be matched

The first three were found by reading the source during the initial implementation; the remainder were found in the audit against the binary (October 2026).

| Area | me-PCR behavior | Consequence of getting it wrong |
|------|-----------------|---------------------------------|
| Hash word | Taken from the 3' end of the primer; moves 5'-wards past ambiguous bases | Different candidate positions, different sensitivity to mismatches |
| Size ranges | Mean of the range, margin widened by `high - mean + 1` | Products at the ends of the range missed |
| Strand orientation | Forward record stores primer1 and the reverse complement of primer2 | All `(+)` hits missed |
| Annotation column | Omitted when the line has no fifth field; text after the fourth tab reported verbatim | Extra empty column; truncated annotations |
| FASTA alphabet | `GATCNBDHKMRSVWY-`; X and U dropped, `-` kept | Coordinates shifted in sequences containing X or gaps |
| IUPAC matching | Primer code matches the bases it expands to; a sequence code never matches a plain primer base; only primers containing codes use this path | Spurious hits across N in the genome |
| U in primers | Not hashable and never equal to T | Spurious hits in IUPAC mode |
| Hit order | Scan order per chunk; records sharing a hash checked most-recent-first | Different line order |
| Chunk overlap | Hits ending inside the leading overlap of a later chunk are dropped (chunk-relative coordinates) | Duplicate hits with `T > 1` |
| Chunk size | Rounded up; threading decided from the FASTA file size | Different chunk boundaries |
| STS line length | Error above the `S=` limit | Long lines accepted |
| Unhashable forward primer | Whole marker skipped | Extra `(-)` record |
| Empty defline | Empty label | Exception |
| `-margin N` | Accepted | Rejected |

## Deliberate differences

- Size ranges (`150-250`) and `-` in the size column are accepted. me-PCR 1.0.6 parses the size with `strtol` and rejects any trailing text, so the range code in the same function is unreachable and the binary exits with "size should be integer".
- Markers that follow a comment or blank line are reported correctly. me-PCR 1.0.6 does not advance its file offset for such lines and reports the wrong line, or nothing, for the next marker.
- Diagnostics are written through Python logging rather than to stderr with fixed wording.

Neither difference affects output for inputs that me-PCR accepts and reports correctly.

## Reproducing the comparison

```bash
# me-PCR test data
me-PCR test/test.sts test/test.fa > a.txt
merpcr test/test.sts test/test.fa > b.txt
diff a.txt b.txt

# randomized comparison
python scripts/compare_with_mepcr.py --mepcr /path/to/me-PCR -n 80 --seed 2
python scripts/compare_with_mepcr.py --mepcr /path/to/me-PCR -n 20 --pure-python
```

Failing cases are kept in a temporary directory whose path is printed.
