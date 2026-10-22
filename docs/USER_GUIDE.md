# User Guide

This guide describes the input and output formats of merPCR, the search parameters, the search algorithm, and how to move from me-PCR.

## Input files

### STS file

The STS file is tab-delimited. Each line describes one marker:

```
identifier<TAB>forward primer<TAB>reverse primer<TAB>product size[<TAB>annotation]
```

```
D1S2893	TGAGTCAGATGTTTGATTTTG	ATGCCACATCAACTTATACTG	126	1p36.33
D1S468	GAATGAACAGAGATGATGCCT	CACACACACACACACCACAC	142-148	1p36.22
AFM248yg9	GCTAAAAATACACGGATGG	TGCAAGACTGCGTCTC	193	(D17S932) Chr.17, 63.7 cM
```

- Both primers are given 5' to 3'. Lower-case letters are accepted. IUPAC ambiguity codes are allowed and are interpreted when `-I 1` is set.
- The product size is the length of the amplicon including both primers. A value of `0` or `-` means unknown, and the default size (`-Z`, 240 bp) is used. A range such as `142-148` is replaced by its mean, and the margin for that marker is widened by `high - mean + 1` so that the whole range is covered.
- If the sum of the primer lengths exceeds the product size, the size is raised to that sum and a warning is logged.
- A primer shorter than the word size cannot be searched for; the marker is skipped with a warning. A marker whose forward primer contains no window of `W` unambiguous bases is skipped as well.
- The annotation is everything after the fourth tab, including further tabs, and is reproduced verbatim in the output.
- Lines starting with `#` and blank lines are ignored. Lines longer than the `-S` limit (1022 characters by default) are an error, as in me-PCR.

### FASTA file

Standard FASTA with one or more records. The sequence label is the text between `>` and the first space or tab. Sequence lines are upper-cased and filtered with the me-PCR nucleotide alphabet `GATCNBDHKMRSVWY-`: other characters, including `X` and `U`, are removed, while the gap character `-` is kept and occupies a position.

## Output

One line per hit, tab-delimited:

```
label<TAB>start..end<TAB>identifier[<TAB>annotation]<TAB>(strand)
```

Positions are 1-based and inclusive and span the whole product from the first base of one primer to the last base of the other. The annotation column is present only if the STS line had one. The strand is `(+)` when the forward primer matches the given strand and the reverse primer matches the opposite strand, and `(-)` for the reverse arrangement.

Hits are reported per sequence in the order the scan finds them. Because the scan triggers on the 3'-most word of a primer, a marker with a shorter primer can be reported before a marker that starts earlier, and markers sharing a hash word are checked most-recent-first. This matches me-PCR; sort the output by position if a positional order is needed.

## Parameters

| Option | Default | Range | Description |
|--------|---------|-------|-------------|
| `-M`, `--margin` | 50 | 0 to 10000 | Products from `size - M` to `size + M` are accepted |
| `-N`, `--mismatches` | 0 | 0 to 10 | Mismatches allowed in each primer |
| `-W`, `--wordsize` | 11 | 3 to 16 | Length of the hash word taken from the 3' end of each primer |
| `-X`, `--three-prime-match` | 1 | 0 or more | Number of 3'-terminal bases in which no mismatch is allowed |
| `-I`, `--iupac` | 0 | 0 or 1 | Interpret IUPAC codes in primers |
| `-Z`, `--default-pcr-size` | 240 | 1 to 10000 | Product size used for markers with size 0 or `-` |
| `-T`, `--threads` | 1 | 1 or more | Worker processes per sequence |
| `-S`, `--max-sts-line-length` | 1022 | 1 or more | Maximum length of an STS line |
| `-O`, `--output` | stdout | | Output file |
| `-Q`, `--quiet` | 1 | 0 or 1 | 0 prints progress and statistics to stderr |
| `--debug` | | | Debug logging |

Every option can also be written in the me-PCR form `M=50`, and the legacy `-margin 50` is accepted.

### Choosing parameters

- The word size controls sensitivity to mismatches near the 3' end of the forward primer: the `W` 3'-most bases of the primer that triggers the hash lookup must match exactly, independent of `-N`. Smaller values find more candidates and cost more time.
- With `-N 0` the `-X` setting has no effect. With `-N > 0`, `-X 1` (the default) still requires the terminal base of each primer to match, which mirrors the requirement for polymerase extension.
- In IUPAC mode a code in a primer matches each base it stands for (R matches A or G, N matches any base). An ambiguity code in the sequence does not match a plain primer base, so an N in the genome counts as a mismatch. Only primers that contain codes are compared this way; others are compared exactly.
- The margin is applied symmetrically around the expected size. me-PCR tests the expected size first and then increasing deviations, and all positions that match are reported, so one marker can produce several hits with the same start.

## Search algorithm

1. Every primer is reduced to a 2-bit hash of its 3'-most `W` unambiguous bases. For each marker two records are stored: forward primer followed by the reverse complement of the reverse primer (strand `+`), and reverse primer followed by the reverse complement of the forward primer (strand `-`).
2. The sequence is scanned once with a rolling hash. Windows containing an ambiguous base are skipped. Each hash hit is verified by comparing the whole first primer, allowing `N` mismatches outside the protected 3' region.
3. The second primer is searched at the expected distance, then at distances `size - 1`, `size + 1`, `size - 2`, ... up to the margin, clipped so that the primers do not overlap and the product stays within the sequence.
4. With `T > 1` the sequence is cut into `T` chunks that overlap by `max product size + margin - 1` bases. Hits whose end lies within the leading overlap of a chunk other than the first are duplicates and are dropped. Files smaller than 100 kB are processed with one worker.

## Moving from me-PCR

Replace the executable name; the arguments and the output are unchanged:

```bash
me-PCR markers.sts genome.fa M=50 N=1 > hits.txt
merpcr markers.sts genome.fa M=50 N=1 > hits.txt
```

Differences that may be noticed:

- Diagnostics go through Python logging and are more descriptive. Progress messages appear with `-Q 0` or `--debug`.
- Size ranges and `-` in the size column are accepted. me-PCR 1.0.6 exits with a format error for both.
- Markers that follow a comment line are reported correctly. me-PCR 1.0.6 reads the wrong line when reporting them.
- Threads are worker processes. Results do not depend on `-T`.

To verify equivalence on your own data, run both programs and compare the outputs with `diff`, or use `scripts/compare_with_mepcr.py` for a randomized comparison.

## Troubleshooting

- No hits: increase `-M`, allow mismatches with `-N`, set `-X 0`, or enable `-I 1` if primers contain ambiguity codes. Check that the STS file is tab-delimited.
- Too many hits: lower `-M`, set `-N 0`, raise `-X`, or review primer specificity.
- Format error on an STS line: the line has fewer than four tab-separated fields or exceeds the `-S` limit.
- Slow runs: use `-T` with the number of available cores on large sequences, and make sure the Cython extension is in use (`--debug` logs whether it was loaded).
