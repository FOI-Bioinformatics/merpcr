# API Reference

The public interface is the `MerPCR` class and the data classes exported from the `merpcr` package.

```python
from merpcr import MerPCR, STSRecord, FASTARecord, STSHit
```

## MerPCR

```python
MerPCR(
    wordsize=11,            # hash word size, 3 to 16
    margin=50,              # size tolerance in bp, 0 to 10000
    mismatches=0,           # mismatches per primer, 0 to 10
    three_prime_match=1,    # protected 3' bases, 0 or more
    iupac_mode=0,           # 1 interprets IUPAC codes in primers
    default_pcr_size=240,   # size for markers given as 0 or "-", 1 to 10000
    threads=1,              # worker processes per sequence
    max_sts_line_length=1022,
)
```

Parameters outside their range raise `ValueError`.

### load_sts_file(filename) -> bool

Reads a tab-delimited STS file and builds the hash table. Returns `False`, after logging an error, when the file is empty, a line has fewer than four fields, or a line exceeds `max_sts_line_length`. Markers with unusable primers are skipped with a warning. After loading, `sts_records` holds two `STSRecord` objects per usable marker (strands `+` and `-`), `sts_table` maps hash values to lists of records, and `max_pcr_size` is the largest product size seen.

### load_fasta_file(filename) -> list[FASTARecord]

Reads all records from a FASTA file. Sequences are upper-cased and restricted to the alphabet `GATCNBDHKMRSVWY-`. The file size is remembered and used to decide whether `threads` may exceed one, as in me-PCR.

### find_hits(fasta_records) -> list[SearchHit]

Searches every record and returns the hits in output order as `SearchHit` objects (see below). The list is also stored in `hits`.

### search(fasta_records, output_file=None) -> int

Calls `find_hits` and writes one line per hit to `output_file`, or to standard output when the argument is `None` or `"stdout"`. Returns the number of hits, which is also stored in `total_hits`.

## Data classes

### STSRecord

| Field | Type | Meaning |
|-------|------|---------|
| `id` | str | Identifier from the first column |
| `primer1` | str | Primer searched on the given strand |
| `primer2` | str | Reverse complement of the partner primer, searched downstream |
| `pcr_size` | int | Expected product size after range and minimum adjustments |
| `alias` | str or None | Text after the fourth tab; `None` when the line had no fifth field |
| `offset` | int | Line number in the STS file |
| `hash_offset` | int | Position of the hash word within `primer1` |
| `direct` | str | `"+"` or `"-"` |
| `ambig_primer` | int | Bit 1: `primer1` contains IUPAC codes; bit 2: `primer2` does (set only in IUPAC mode) |
| `margin` | int | Margin for this marker, including the range adjustment |

### FASTARecord

`defline`, `sequence` and `label`. The label defaults to the first word of the defline after `>`; it is empty if the defline has no text.

### SearchHit

Returned by `find_hits`. Fields `seq_label`, `start`, `end` (1-based, inclusive) and `sts`; properties `sts_id`, `alias`, `strand` and `size`; `format()` returns the output line.

### STSHit

Internal representation with 0-based positions `pos1`, `pos2` and the matching `STSRecord`.

## Lower-level functions

- `merpcr.core.utils.hash_value(primer, wordsize)` returns `(offset, hash)` for the 3'-most unambiguous word, or `(-1, 0)`.
- `merpcr.core.utils.reverse_complement(sequence)` handles IUPAC codes and maps unknown characters to `N`.
- `merpcr.io.fasta.FASTALoader.load_file(filename)` is the FASTA reader used by `MerPCR`.
- `merpcr.io.sts.STSLoader(wordsize, margin, default_pcr_size, iupac_mode=0).load_file(filename)` returns `(records, table, max_pcr_size)` using the engine's parser.

Methods of `MerPCR` whose names start with an underscore (`_hash_value`, `_compare_seqs`, `_match_sts`, `_process_thread`) implement the algorithm and may change between releases.

## Example

```python
import logging
from merpcr import MerPCR

logging.basicConfig(level=logging.INFO)

engine = MerPCR(margin=100, mismatches=1, iupac_mode=1, threads=4)
if not engine.load_sts_file("markers.sts"):
    raise SystemExit("could not read markers.sts")

records = engine.load_fasta_file("assembly.fa")
for hit in engine.find_hits(records):
    print(hit.seq_label, hit.start, hit.end, hit.sts_id, hit.strand, hit.size)

n = engine.search(records, output_file="hits.txt")  # same hits, written as text
```
