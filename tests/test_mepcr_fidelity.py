"""
Regression tests for behaviors verified against the me-PCR 1.0.6 binary.

Each test documents a point where merPCR previously diverged from me-PCR.
The expected values were taken from running me-PCR on the same input.
"""

import io
import os
import random
import tempfile
from contextlib import redirect_stdout

import pytest

from merpcr import MerPCR
from merpcr.cli import convert_mepcr_arguments
from merpcr.core.models import FASTARecord, STSRecord
from merpcr.core.utils import hash_value
from merpcr.io.fasta import FASTALoader
from merpcr.io.sts import STSLoader

P1 = "ACGTTGCAGGCTAGCTTACG"
P2 = "GGATCCTTAGCGATCGGTAC"


def rc(seq):
    return seq[::-1].translate(str.maketrans("ACGT", "TGCA"))


def planted_sequence(seed=1, length=1000, p1=P1, p2=P2, start=100, size=200):
    rng = random.Random(seed)
    seq = [rng.choice("ACGT") for _ in range(length)]
    seq[start : start + len(p1)] = p1
    seq[start + size - len(p2) : start + size] = rc(p2)
    return "".join(seq)


@pytest.fixture
def tmpdir_path():
    with tempfile.TemporaryDirectory() as d:
        yield d


def write(path, text):
    with open(path, "w") as fh:
        fh.write(text)
    return path


def run_search(sts_text, fasta_text, tmp, **kwargs):
    sts = write(os.path.join(tmp, "t.sts"), sts_text)
    fa = write(os.path.join(tmp, "t.fa"), fasta_text)
    engine = MerPCR(**kwargs)
    assert engine.load_sts_file(sts)
    records = engine.load_fasta_file(fa)
    buf = io.StringIO()
    with redirect_stdout(buf):
        engine.search(records)
    return buf.getvalue().splitlines()


class TestOutputFormat:
    def test_no_alias_field_is_omitted(self, tmpdir_path):
        # me-PCR: "seq1\t101..300\tS1\t(+)"
        lines = run_search(
            f"S1\t{P1}\t{P2}\t200\n", ">seq1\n" + planted_sequence() + "\n", tmpdir_path
        )
        assert lines == ["seq1\t101..300\tS1\t(+)"]

    def test_empty_alias_after_trailing_tab(self, tmpdir_path):
        lines = run_search(
            f"S1\t{P1}\t{P2}\t200\t\n", ">seq1\n" + planted_sequence() + "\n", tmpdir_path
        )
        assert lines == ["seq1\t101..300\tS1\t\t(+)"]

    def test_alias_keeps_everything_after_fourth_tab(self, tmpdir_path):
        lines = run_search(
            f"S1\t{P1}\t{P2}\t200\tALIAS1\textra\n",
            ">seq1\n" + planted_sequence() + "\n",
            tmpdir_path,
        )
        assert lines == ["seq1\t101..300\tS1\tALIAS1\textra\t(+)"]

    def test_alias_default_is_none(self):
        sts = STSRecord(id="x", primer1=P1, primer2=P2, pcr_size=200)
        assert sts.alias is None


class TestFastaCharacters:
    def test_x_is_dropped_and_gap_is_kept(self, tmpdir_path):
        seq = planted_sequence()
        with_x = seq[:50] + "X" + seq[50:]
        with_gap = seq[:50] + "-" + seq[50:]
        assert run_search(f"S1\t{P1}\t{P2}\t200\n", ">sx\n" + with_x + "\n", tmpdir_path) == [
            "sx\t101..300\tS1\t(+)"
        ]
        assert run_search(f"S1\t{P1}\t{P2}\t200\n", ">sd\n" + with_gap + "\n", tmpdir_path) == [
            "sd\t102..301\tS1\t(+)"
        ]

    def test_loader_alphabet(self, tmpdir_path):
        fa = write(os.path.join(tmpdir_path, "a.fa"), ">a\nacgtXU-nrykmswbdhv*0\n")
        records = FASTALoader.load_file(fa)
        assert records[0].sequence == "ACGT-NRYKMSWBDHV"

    def test_empty_defline_gives_empty_label(self):
        assert FASTARecord(defline=">", sequence="ACGT").label == ""
        assert FASTARecord(defline=">  name rest", sequence="ACGT").label == "name"


class TestIupacSemantics:
    def test_sequence_ambiguity_does_not_match_plain_primer(self, tmpdir_path):
        seq = list(planted_sequence())
        seq[105] = "N"
        fasta = ">sn\n" + "".join(seq) + "\n"
        assert run_search(f"S1\t{P1}\t{P2}\t200\n", fasta, tmpdir_path, iupac_mode=1) == []
        # With one mismatch allowed the N counts as a mismatch, as in me-PCR
        assert run_search(
            f"S1\t{P1}\t{P2}\t200\n", fasta, tmpdir_path, iupac_mode=1, mismatches=1
        ) == ["sn\t101..300\tS1\t(+)"]

    def test_primer_ambiguity_matches_sequence_base(self, tmpdir_path):
        p1n = P1[:5] + "N" + P1[6:]
        fasta = ">s\n" + planted_sequence() + "\n"
        assert run_search(f"S1\t{p1n}\t{P2}\t200\n", fasta, tmpdir_path) == []
        assert run_search(f"S1\t{p1n}\t{P2}\t200\n", fasta, tmpdir_path, iupac_mode=1) == [
            "s\t101..300\tS1\t(+)"
        ]

    def test_u_in_primer_is_not_t(self, tmpdir_path):
        i = P1.index("T")
        p1u = P1[:i] + "U" + P1[i + 1 :]
        fasta = ">s\n" + planted_sequence() + "\n"
        assert run_search(f"S1\t{p1u}\t{P2}\t200\n", fasta, tmpdir_path) == []
        assert run_search(f"S1\t{p1u}\t{P2}\t200\n", fasta, tmpdir_path, iupac_mode=1) == []

    def test_ambig_flags_recorded_only_in_iupac_mode(self, tmpdir_path):
        sts = write(os.path.join(tmpdir_path, "a.sts"), f"S1\t{P1[:-1]}R\t{P2}\t200\n")
        plain = MerPCR()
        plain.load_sts_file(sts)
        assert all(r.ambig_primer == 0 for r in plain.sts_records)
        iupac = MerPCR(iupac_mode=1)
        iupac.load_sts_file(sts)
        by_dir = {r.direct: r for r in iupac.sts_records}
        assert by_dir["+"].ambig_primer == 1
        assert by_dir["-"].ambig_primer == 2


class TestHitOrder:
    def test_hits_are_reported_in_scan_order(self, tmpdir_path):
        # The 15-mer primer of B reaches its hash word before A's 20-mer does,
        # so me-PCR reports B first although it starts later.
        rng = random.Random(5)
        seq = [rng.choice("ACGT") for _ in range(2000)]
        qa = P2
        seq[100:120] = P1
        seq[280:300] = rc(qa)
        pb = P1[3:18]
        qb = "TTGACCGGTTAACCGGTTAA"
        seq[400:420] = rc(qb)
        fasta = ">ord\n" + "".join(seq) + "\n"
        lines = run_search(f"A\t{P1}\t{qa}\t200\tA\nB\t{pb}\t{qb}\t317\tB\n", fasta, tmpdir_path)
        assert lines == ["ord\t104..420\tB\tB\t(+)", "ord\t101..300\tA\tA\t(+)"]

    def test_later_sts_with_same_hash_is_reported_first(self, tmpdir_path):
        rng = random.Random(5)
        seq = [rng.choice("ACGT") for _ in range(2000)]
        qc = "CCGGAATTCCGGAATTCCAA"
        seq[100:120] = P1
        seq[280:300] = rc(P2)
        seq[350:370] = rc(qc)
        fasta = ">same\n" + "".join(seq) + "\n"
        lines = run_search(f"A\t{P1}\t{P2}\t200\tA\nC\t{P1}\t{qc}\t270\tC\n", fasta, tmpdir_path)
        assert lines == ["same\t101..370\tC\tC\t(+)", "same\t101..300\tA\tA\t(+)"]


class TestThreading:
    def test_no_duplicate_hits_across_chunks(self, tmpdir_path):
        rng = random.Random(3)
        length = 400000
        seq = [rng.choice("ACGT") for _ in range(length)]
        sts_lines = []
        for n in range(30):
            p1 = "".join(rng.choice("ACGT") for _ in range(20))
            p2 = "".join(rng.choice("ACGT") for _ in range(20))
            size = rng.randint(100, 400)
            for _ in range(3):
                pos = rng.randint(0, length - size - 1)
                seq[pos : pos + 20] = p1
                seq[pos + size - 20 : pos + size] = rc(p2)
            sts_lines.append(f"STS{n}\t{p1}\t{p2}\t{size}\tA{n}")
        fasta = ">big\n" + "".join(seq) + "\n"
        sts_text = "\n".join(sts_lines) + "\n"
        single = run_search(sts_text, fasta, tmpdir_path, threads=1)
        multi = run_search(sts_text, fasta, tmpdir_path, threads=4)
        assert len(single) == len(set(single))
        assert multi == single


class TestStsFile:
    def test_line_length_limit_is_enforced(self, tmpdir_path):
        sts = write(os.path.join(tmpdir_path, "l.sts"), f"A\t{P1}\t{P2}\t200\t" + "x" * 1100 + "\n")
        assert MerPCR().load_sts_file(sts) is False
        assert MerPCR(max_sts_line_length=2000).load_sts_file(sts) is True

    def test_sts_without_hash_for_first_primer_is_skipped_entirely(self, tmpdir_path):
        sts = write(os.path.join(tmpdir_path, "n.sts"), f"A\t{'N' * 20}\t{P2}\t200\n")
        engine = MerPCR()
        assert engine.load_sts_file(sts)
        assert engine.sts_records == []

    def test_sts_loader_matches_engine(self, tmpdir_path):
        sts = write(
            os.path.join(tmpdir_path, "s.sts"), f"A\t{P1}\t{P2}\t200\tA\nB\t{P2}\t{P1}\t150\n"
        )
        records, table, max_size = STSLoader(
            wordsize=11, margin=50, default_pcr_size=240
        ).load_file(sts)
        engine = MerPCR()
        engine.load_sts_file(sts)
        assert [vars(r) for r in records] == [vars(r) for r in engine.sts_records]
        assert max_size == engine.max_pcr_size
        assert set(table) == set(engine.sts_table)


class TestUtils:
    def test_hash_value_uses_three_prime_end(self):
        engine = MerPCR(wordsize=11)
        for primer in ["ATCGATCGATCG", "NATCGTTACGGATCC", "ACGTNACGTACGTACGT", "ACGUACGTACGTAC"]:
            assert hash_value(primer, 11) == engine._hash_value(primer)

    def test_u_is_not_hashable(self):
        assert hash_value("ACGUACGUACG", 11) == (-1, 0)


class TestCli:
    def test_legacy_margin_flag(self):
        assert convert_mepcr_arguments(["a.sts", "b.fa", "-margin", "5"]) == [
            "a.sts",
            "b.fa",
            "-M",
            "5",
        ]


class TestFindHits:
    def test_find_hits_returns_structured_results(self, tmpdir_path):
        sts = write(os.path.join(tmpdir_path, "f.sts"), f"S1\t{P1}\t{P2}\t200\tAL\n")
        fa = write(os.path.join(tmpdir_path, "f.fa"), ">seq1 desc\n" + planted_sequence() + "\n")
        engine = MerPCR()
        assert engine.load_sts_file(sts)
        hits = engine.find_hits(engine.load_fasta_file(fa))
        assert len(hits) == 1
        hit = hits[0]
        assert (hit.seq_label, hit.start, hit.end) == ("seq1", 101, 300)
        assert hit.sts_id == "S1" and hit.alias == "AL" and hit.strand == "+"
        assert hit.size == 200
        assert hit.format() == "seq1\t101..300\tS1\tAL\t(+)"
        assert engine.hits == hits

    def test_search_output_matches_find_hits(self, tmpdir_path):
        sts = write(
            os.path.join(tmpdir_path, "f.sts"), f"S1\t{P1}\t{P2}\t200\nS2\t{P2}\t{P1}\t200\n"
        )
        fa = write(os.path.join(tmpdir_path, "f.fa"), ">seq1\n" + planted_sequence() + "\n")
        engine = MerPCR()
        engine.load_sts_file(sts)
        records = engine.load_fasta_file(fa)
        out = os.path.join(tmpdir_path, "out.txt")
        count = engine.search(records, output_file=out)
        with open(out) as fh:
            lines = fh.read().splitlines()
        assert count == len(lines) == len(engine.hits)
        assert lines == [h.format() for h in engine.find_hits(records)]


class TestCompressedInput:
    @pytest.mark.parametrize("suffix", [".gz", ".bz2", ".xz"])
    def test_compressed_sts_and_fasta_give_same_hits(self, tmpdir_path, suffix):
        import bz2
        import gzip
        import lzma

        openers = {".gz": gzip.open, ".bz2": bz2.open, ".xz": lzma.open}
        sts_text = f"S1\t{P1}\t{P2}\t200\tAL\n"
        fasta_text = ">seq1\n" + planted_sequence() + "\n"
        plain = run_search(sts_text, fasta_text, tmpdir_path)

        sts = os.path.join(tmpdir_path, "c.sts" + suffix)
        fa = os.path.join(tmpdir_path, "c.fa" + suffix)
        with openers[suffix](sts, "wt") as fh:
            fh.write(sts_text)
        with openers[suffix](fa, "wt") as fh:
            fh.write(fasta_text)

        engine = MerPCR()
        assert engine.load_sts_file(sts)
        records = engine.load_fasta_file(fa)
        assert engine._fasta_file_size is None
        assert (
            [h.format() for h in engine.find_hits(records)]
            == plain
            == ["seq1\t101..300\tS1\tAL\t(+)"]
        )

    def test_plain_file_still_uses_size_on_disk(self, tmpdir_path):
        fa = write(os.path.join(tmpdir_path, "p.fa"), ">seq1\n" + planted_sequence() + "\n")
        engine = MerPCR()
        engine.load_fasta_file(fa)
        assert engine._fasta_file_size == os.path.getsize(fa)
