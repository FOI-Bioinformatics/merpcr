"""
Tests for optional hydrolysis-probe support (--probe-column).

A marker may carry a probe in column 5 of the STS file. The probe is checked
on each candidate amplicon after the primer search. Lines of markers without
a probe must stay identical to the me-PCR output format.
"""

import random
import sys

import pytest

import merpcr.core.engine as engine_module
from merpcr import FASTARecord, MerPCR
from merpcr.cli import convert_mepcr_arguments, main
from merpcr.io.sts import STSLoader

PRIMER1 = "GATTACAGGCTAGCTTAC"
PRIMER2 = "CCATGGTTAGCAATCGGA"  # as written in the STS file
PROBE = "TTGACCGATAGGCATCCAGT"
FILL1 = "CAGTCTAAGC"
FILL2 = "TGCAATGTCA"

_COMPLEMENT = str.maketrans("ACGTRYKMBVDHN", "TGCAYRMKVBHDN")


def rc(seq: str) -> str:
    return seq.translate(_COMPLEMENT)[::-1]


def random_dna(length: int, seed: int) -> str:
    rng = random.Random(seed)
    return "".join(rng.choice("ACGT") for _ in range(length))


def make_amplicon(inner: str = PROBE) -> str:
    """Forward-strand amplicon: primer1, filler, inner part, filler, rc(primer2)."""
    return PRIMER1 + FILL1 + inner + FILL2 + rc(PRIMER2)


AMPLICON = make_amplicon()
SIZE = len(AMPLICON)  # 76
FLANK_L = random_dna(100, 1)
FLANK_R = random_dna(100, 2)


def make_sequence(amplicon: str = AMPLICON) -> str:
    return FLANK_L + amplicon + FLANK_R


def sts_line(probe: str = PROBE, alias: str = "my alias", sts_id: str = "ASSAY") -> str:
    fields = [sts_id, PRIMER1, PRIMER2, str(SIZE)]
    if probe is not None:
        fields.append(probe)
        if alias is not None:
            fields.append(alias)
    return "\t".join(fields) + "\n"


def write_sts(tmp_path, *lines: str, name: str = "assays.sts") -> str:
    path = tmp_path / name
    path.write_text("".join(lines))
    return str(path)


def run(tmp_path, sequence: str, *lines: str, **options):
    """Load the STS lines with probe_column set and return the hits."""
    options.setdefault("probe_column", True)
    engine = MerPCR(**options)
    assert engine.load_sts_file(write_sts(tmp_path, *lines))
    return engine.find_hits([FASTARecord(defline=">seq", sequence=sequence)])


class TestProbeParsing:
    def test_flag_off_keeps_column_five_in_alias(self, tmp_path):
        engine = MerPCR()
        assert engine.load_sts_file(write_sts(tmp_path, sts_line()))
        for record in engine.sts_records:
            assert record.probe is None
            assert record.alias == f"{PROBE}\tmy alias"

    def test_flag_on_reads_probe_and_alias(self, tmp_path):
        engine = MerPCR(probe_column=True)
        assert engine.load_sts_file(write_sts(tmp_path, sts_line()))
        assert [r.direct for r in engine.sts_records] == ["+", "-"]
        for record in engine.sts_records:
            assert record.probe == PROBE
            assert record.alias == "my alias"

    def test_alias_keeps_everything_after_fifth_tab(self, tmp_path):
        engine = MerPCR(probe_column=True)
        assert engine.load_sts_file(write_sts(tmp_path, sts_line(alias="a\tb c\t")))
        assert engine.sts_records[0].alias == "a\tb c\t"

    @pytest.mark.parametrize("placeholder", ["-", ""])
    def test_placeholder_means_no_probe(self, tmp_path, placeholder):
        engine = MerPCR(probe_column=True)
        assert engine.load_sts_file(write_sts(tmp_path, sts_line(probe=placeholder)))
        assert engine.sts_records[0].probe is None
        assert engine.sts_records[0].alias == "my alias"

    def test_four_column_line_has_no_probe_and_no_alias(self, tmp_path):
        engine = MerPCR(probe_column=True)
        assert engine.load_sts_file(write_sts(tmp_path, sts_line(probe=None)))
        assert engine.sts_records[0].probe is None
        assert engine.sts_records[0].alias is None

    def test_probe_without_alias(self, tmp_path):
        engine = MerPCR(probe_column=True)
        assert engine.load_sts_file(write_sts(tmp_path, sts_line(alias=None)))
        assert engine.sts_records[0].probe == PROBE
        assert engine.sts_records[0].alias is None

    def test_probe_is_upper_cased(self, tmp_path):
        engine = MerPCR(probe_column=True)
        assert engine.load_sts_file(write_sts(tmp_path, sts_line(probe=PROBE.lower())))
        assert engine.sts_records[0].probe == PROBE

    def test_invalid_probe_characters_are_an_error(self, tmp_path):
        engine = MerPCR(probe_column=True)
        assert not engine.load_sts_file(write_sts(tmp_path, sts_line(probe="ACGT ACGTEE")))

    def test_sts_loader_matches_engine(self, tmp_path):
        path = write_sts(tmp_path, sts_line())
        engine = MerPCR(probe_column=True)
        assert engine.load_sts_file(path)
        loader = STSLoader(wordsize=11, margin=50, default_pcr_size=240, probe_column=True)
        records, _, _ = loader.load_file(path)
        assert [vars(r) for r in records] == [vars(r) for r in engine.sts_records]

    def test_probe_mismatches_is_validated(self):
        with pytest.raises(ValueError):
            MerPCR(probe_mismatches=-1)
        with pytest.raises(ValueError):
            MerPCR(probe_mismatches=11)


class TestProbeMatching:
    def test_probe_on_plus_strand(self, tmp_path):
        hits = run(tmp_path, make_sequence(), sts_line())
        assert len(hits) == 1
        hit = hits[0]
        assert (hit.start, hit.end, hit.strand) == (101, 176, "+")
        assert hit.probe_found
        assert (hit.probe_start, hit.probe_end, hit.probe_strand) == (129, 148, "+")
        assert hit.probe_mismatches == 0
        assert hit.probe_seq == PROBE

    def test_probe_on_minus_strand(self, tmp_path):
        hits = run(tmp_path, make_sequence(make_amplicon(rc(PROBE))), sts_line())
        assert len(hits) == 1
        assert (hits[0].probe_start, hits[0].probe_end, hits[0].probe_strand) == (129, 148, "-")
        assert hits[0].probe_seq == rc(PROBE)

    def test_marker_on_minus_strand(self, tmp_path):
        hits = run(tmp_path, make_sequence(rc(AMPLICON)), sts_line())
        assert len(hits) == 1
        hit = hits[0]
        assert (hit.start, hit.end, hit.strand) == (101, 176, "-")
        # Probe strand is given relative to the sequence, not to the marker
        assert (hit.probe_start, hit.probe_end, hit.probe_strand) == (129, 148, "-")

    def test_amplicon_without_probe_is_dropped(self, tmp_path):
        sequence = make_sequence(make_amplicon(random_dna(len(PROBE), 3)))
        assert run(tmp_path, sequence, sts_line()) == []
        # The primer pair alone does give a hit
        assert len(run(tmp_path, sequence, sts_line(probe="-"))) == 1

    def test_probe_mismatch_limit(self, tmp_path):
        mutated = PROBE[:9] + ("A" if PROBE[9] != "A" else "C") + PROBE[10:]
        sequence = make_sequence(make_amplicon(mutated))
        assert run(tmp_path, sequence, sts_line()) == []
        hits = run(tmp_path, sequence, sts_line(), probe_mismatches=1)
        assert len(hits) == 1
        assert hits[0].probe_mismatches == 1
        assert hits[0].probe_seq == mutated

    def test_probe_limit_is_independent_of_primer_limit(self, tmp_path):
        mutated = PROBE[:9] + ("A" if PROBE[9] != "A" else "C") + PROBE[10:]
        sequence = make_sequence(make_amplicon(mutated))
        assert run(tmp_path, sequence, sts_line(), mismatches=2) == []

    def test_degenerate_probe_without_iupac_mode(self, tmp_path):
        # R covers the G at position 6 of the probe site
        assert PROBE[6] == "G"
        degenerate = PROBE[:6] + "R" + PROBE[7:]
        hits = run(tmp_path, make_sequence(), sts_line(probe=degenerate), iupac_mode=0)
        assert len(hits) == 1
        assert hits[0].probe_mismatches == 0

    def test_degenerate_probe_on_minus_strand(self, tmp_path):
        degenerate = PROBE[:6] + "R" + PROBE[7:]
        sequence = make_sequence(make_amplicon(rc(PROBE)))
        hits = run(tmp_path, sequence, sts_line(probe=degenerate))
        assert len(hits) == 1
        assert (hits[0].probe_strand, hits[0].probe_mismatches) == ("-", 0)

    def test_ambiguity_code_in_target_is_a_mismatch(self, tmp_path):
        with_n = PROBE[:9] + "N" + PROBE[10:]
        sequence = make_sequence(make_amplicon(with_n))
        assert run(tmp_path, sequence, sts_line()) == []
        hits = run(tmp_path, sequence, sts_line(), probe_mismatches=1)
        assert [h.probe_mismatches for h in hits] == [1]

    def test_plus_probe_may_not_overlap_left_primer(self, tmp_path):
        # Written as the forward strand, the probe anneals to the same strand
        # as the left primer, so a site reaching into that primer is rejected.
        assert run(tmp_path, make_sequence(), sts_line(probe=AMPLICON[10:30])) == []

    def test_minus_probe_may_overlap_left_primer(self, tmp_path):
        hits = run(tmp_path, make_sequence(), sts_line(probe=rc(AMPLICON[10:30])))
        assert len(hits) == 1
        assert (hits[0].probe_start, hits[0].probe_end, hits[0].probe_strand) == (111, 130, "-")

    def test_minus_probe_may_not_overlap_right_primer(self, tmp_path):
        assert run(tmp_path, make_sequence(), sts_line(probe=rc(AMPLICON[50:70]))) == []

    def test_plus_probe_may_overlap_right_primer(self, tmp_path):
        hits = run(tmp_path, make_sequence(), sts_line(probe=AMPLICON[50:70]))
        assert len(hits) == 1
        assert (hits[0].probe_start, hits[0].probe_end, hits[0].probe_strand) == (151, 170, "+")

    def test_leftmost_site_is_reported_for_equal_mismatches(self, tmp_path):
        amplicon = make_amplicon(PROBE + "ACCA" + PROBE)
        line = "\t".join(["ASSAY", PRIMER1, PRIMER2, str(len(amplicon)), PROBE]) + "\n"
        hits = run(tmp_path, make_sequence(amplicon), line)
        assert len(hits) == 1
        assert (hits[0].probe_start, hits[0].probe_end) == (129, 148)

    def test_site_with_fewest_mismatches_is_reported(self, tmp_path):
        mutated = "A" + PROBE[1:] if PROBE[0] != "A" else "C" + PROBE[1:]
        amplicon = make_amplicon(mutated + "ACCA" + rc(PROBE))
        line = "\t".join(["ASSAY", PRIMER1, PRIMER2, str(len(amplicon)), PROBE]) + "\n"
        hits = run(tmp_path, make_sequence(amplicon), line, probe_mismatches=2)
        assert len(hits) == 1
        assert (hits[0].probe_strand, hits[0].probe_mismatches) == ("-", 0)

    def test_probe_longer_than_amplicon(self, tmp_path):
        assert run(tmp_path, make_sequence(), sts_line(probe=random_dna(90, 4))) == []

    def test_lower_case_sequence(self, tmp_path):
        hits = run(tmp_path, make_sequence().lower(), sts_line())
        assert len(hits) == 1
        assert hits[0].probe_seq == PROBE


class TestProbeReporting:
    def test_output_line_with_probe(self, tmp_path):
        hits = run(tmp_path, make_sequence(), sts_line())
        assert [h.format() for h in hits] == [
            "seq\t101..176\tASSAY\tmy alias\t(+)\tprobe=129..148(+)\tmm=0"
        ]

    def test_output_line_with_probe_and_no_alias(self, tmp_path):
        hits = run(tmp_path, make_sequence(), sts_line(alias=None))
        assert [h.format() for h in hits] == ["seq\t101..176\tASSAY\t(+)\tprobe=129..148(+)\tmm=0"]

    def test_keep_probeless(self, tmp_path):
        sequence = make_sequence(make_amplicon(random_dna(len(PROBE), 3)))
        hits = run(tmp_path, sequence, sts_line(), require_probe=False)
        assert len(hits) == 1
        hit = hits[0]
        assert not hit.probe_found
        assert hit.probe_start is None and hit.probe_end is None
        assert hit.probe_strand is None and hit.probe_mismatches is None
        assert hit.probe_seq == ""
        assert hit.format() == "seq\t101..176\tASSAY\tmy alias\t(+)\tprobe=none"

    def test_marker_without_probe_keeps_mepcr_line(self, tmp_path):
        hits = run(tmp_path, make_sequence(), sts_line(probe="-"))
        assert [h.format() for h in hits] == ["seq\t101..176\tASSAY\tmy alias\t(+)"]
        assert not hits[0].probe_found

    def test_mixed_file(self, tmp_path):
        lines = [sts_line(sts_id="WITH"), sts_line(probe="-", sts_id="WITHOUT", alias=None)]
        with_flag = run(tmp_path, make_sequence(), *lines)
        plain_line = "\t".join(["WITHOUT", PRIMER1, PRIMER2, str(SIZE)]) + "\n"
        plain = run(tmp_path, make_sequence(), plain_line, probe_column=False)
        formatted = [h.format() for h in with_flag]
        assert sorted(formatted) == sorted(
            ["seq\t101..176\tWITH\tmy alias\t(+)\tprobe=129..148(+)\tmm=0"]
            + [h.format() for h in plain]
        )
        assert [h.format() for h in plain] == ["seq\t101..176\tWITHOUT\t(+)"]

    def test_search_counts_reported_lines(self, tmp_path, capsys):
        with_probe = make_sequence()
        without_probe = make_sequence(make_amplicon(random_dna(len(PROBE), 3)))
        engine = MerPCR(probe_column=True)
        assert engine.load_sts_file(write_sts(tmp_path, sts_line()))
        records = [
            FASTARecord(defline=">a", sequence=with_probe),
            FASTARecord(defline=">b", sequence=without_probe),
        ]
        assert engine.search(records) == 1
        assert engine.total_hits == 1
        assert capsys.readouterr().out == (
            "a\t101..176\tASSAY\tmy alias\t(+)\tprobe=129..148(+)\tmm=0\n"
        )

    def test_hits_are_stored_on_engine(self, tmp_path):
        engine = MerPCR(probe_column=True)
        assert engine.load_sts_file(write_sts(tmp_path, sts_line()))
        records = [FASTARecord(defline=">seq", sequence=make_sequence())]
        hits = engine.find_hits(records)
        assert engine.hits == hits
        assert engine.find_hits(records) == hits


def _large_sequence() -> str:
    """About 150 kb with amplicons on both strands, one of them without a probe site."""
    parts = [random_dna(20000, 10), AMPLICON]
    parts += [random_dna(30000, 11), rc(AMPLICON)]
    parts += [random_dna(30000, 12), make_amplicon(random_dna(len(PROBE), 3))]
    parts += [random_dna(30000, 13), make_amplicon(rc(PROBE))]
    parts += [random_dna(40000, 14)]
    return "".join(parts)


class TestProbeInvariance:
    def test_all_placeholders_equal_run_without_column(self, tmp_path):
        sequence = make_sequence()
        with_column = run(tmp_path, sequence, sts_line(probe="-"), sts_line(probe="", alias="x"))
        plain_lines = [
            "\t".join(["ASSAY", PRIMER1, PRIMER2, str(SIZE), "my alias"]) + "\n",
            "\t".join(["ASSAY", PRIMER1, PRIMER2, str(SIZE), "x"]) + "\n",
        ]
        plain = run(tmp_path, sequence, *plain_lines, probe_column=False)
        assert [h.format() for h in with_column] == [h.format() for h in plain]
        assert len(plain) == 2

    @pytest.mark.parametrize("require_probe", [True, False])
    def test_threads_give_same_hits(self, tmp_path, require_probe):
        sequence = _large_sequence()
        single = run(tmp_path, sequence, sts_line(), threads=1, require_probe=require_probe)
        multi = run(tmp_path, sequence, sts_line(), threads=3, require_probe=require_probe)
        assert [h.format() for h in multi] == [h.format() for h in single]
        assert len(single) == (3 if require_probe else 4)

    def test_pure_python_gives_same_hits(self, tmp_path, monkeypatch):
        sequence = _large_sequence()
        fast = run(tmp_path, sequence, sts_line(), require_probe=False)
        monkeypatch.setattr(engine_module, "HAS_CYTHON", False)
        slow = run(tmp_path, sequence, sts_line(), require_probe=False)
        assert [h.format() for h in slow] == [h.format() for h in fast]
        assert len(fast) == 4


class TestProbeCli:
    def _run_cli(self, monkeypatch, capsys, tmp_path, sequence, *extra):
        sts = write_sts(tmp_path, sts_line())
        fasta = tmp_path / "seq.fa"
        fasta.write_text(f">seq\n{sequence}\n")
        monkeypatch.setattr(sys, "argv", ["merpcr", sts, str(fasta), *extra])
        code = main()
        return code, capsys.readouterr().out

    def test_default_keeps_column_five_as_alias(self, monkeypatch, capsys, tmp_path):
        code, out = self._run_cli(monkeypatch, capsys, tmp_path, make_sequence())
        assert code == 0
        assert out == f"seq\t101..176\tASSAY\t{PROBE}\tmy alias\t(+)\n"

    def test_probe_column(self, monkeypatch, capsys, tmp_path):
        code, out = self._run_cli(monkeypatch, capsys, tmp_path, make_sequence(), "--probe-column")
        assert code == 0
        assert out == "seq\t101..176\tASSAY\tmy alias\t(+)\tprobe=129..148(+)\tmm=0\n"

    def test_probe_mismatches_and_keep_probeless(self, monkeypatch, capsys, tmp_path):
        mutated = PROBE[:9] + ("A" if PROBE[9] != "A" else "C") + PROBE[10:]
        sequence = make_sequence(make_amplicon(mutated))
        code, out = self._run_cli(monkeypatch, capsys, tmp_path, sequence, "--probe-column")
        assert (code, out) == (0, "")
        code, out = self._run_cli(
            monkeypatch, capsys, tmp_path, sequence, "--probe-column", "--probe-mismatches", "1"
        )
        assert out == "seq\t101..176\tASSAY\tmy alias\t(+)\tprobe=129..148(+)\tmm=1\n"
        code, out = self._run_cli(
            monkeypatch, capsys, tmp_path, sequence, "--probe-column", "--keep-probeless"
        )
        assert out == "seq\t101..176\tASSAY\tmy alias\t(+)\tprobe=none\n"

    def test_legacy_one_sided_argument_is_converted(self):
        assert convert_mepcr_arguments(["U=1", "a.sts", "b.fa"]) == ["-U", "1", "a.sts", "b.fa"]
