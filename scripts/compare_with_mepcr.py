#!/usr/bin/env python3
"""
Differential test of merPCR against the original me-PCR binary.

Random genomes and STS files are generated, planted with primer pairs, and
searched with both programs under randomly chosen parameters. Any difference
in standard output is reported and the offending input files are kept.

Usage:
    python scripts/compare_with_mepcr.py --mepcr /path/to/me-PCR -n 50
    python scripts/compare_with_mepcr.py --mepcr /path/to/me-PCR --pure-python

Size ranges and "-" sizes are not generated because me-PCR 1.0.6 rejects them
(see docs/VERIFICATION.md). Comment lines are not generated either, since
me-PCR 1.0.6 misreports the STS that follows a comment.
"""

import argparse
import os
import random
import subprocess
import sys
import tempfile

IUPAC = {
    "R": "AG",
    "Y": "CT",
    "M": "AC",
    "K": "GT",
    "S": "CG",
    "W": "AT",
    "B": "CGT",
    "D": "AGT",
    "H": "ACT",
    "V": "ACG",
    "N": "ACGT",
}
COMPL = str.maketrans("ACGTRYMKSWBDHVN", "TGCAYRKMSWVHDBN")


def rc(seq):
    return seq[::-1].translate(COMPL)


def random_primer(rng, iupac):
    length = rng.randint(12, 26)
    bases = [rng.choice("ACGT") for _ in range(length)]
    if iupac and rng.random() < 0.5:
        # Keep the 3' end unambiguous so that a hash can still be formed
        for _ in range(rng.randint(1, 2)):
            pos = rng.randrange(0, max(1, length - 11))
            bases[pos] = rng.choice(list(IUPAC))
    return "".join(bases)


def instantiate(rng, primer):
    """Replace IUPAC codes with one of their bases."""
    return "".join(rng.choice(IUPAC.get(b, b)) for b in primer)


def mutate(rng, seq, count):
    seq = list(seq)
    for _ in range(count):
        pos = rng.randrange(len(seq))
        seq[pos] = rng.choice([b for b in "ACGT" if b != seq[pos]])
    return "".join(seq)


def make_case(rng, workdir, iupac):
    nseq = rng.choice([1, 1, 1, 2, 3])
    records = []
    sts_lines = []
    nsts = rng.randint(3, 25)
    for n in range(nsts):
        p1 = random_primer(rng, iupac)
        p2 = random_primer(rng, iupac)
        size = rng.randint(len(p1) + len(p2), 500)
        if rng.random() < 0.15:
            size = rng.choice([0, len(p1) + len(p2) - 3])  # default size / too small
        alias = rng.choice(["", "", f"\tA{n}", f"\tA{n}\textra field", "\t"])
        sts_lines.append(f"STS{n}\t{p1}\t{p2}\t{size}{alias}")
    for s in range(nseq):
        length = rng.choice([3000, 20000, 60000, 150000, 400000])
        genome = [rng.choice("ACGT") for _ in range(length)]
        # a few N blocks and odd characters
        for _ in range(rng.randint(0, 3)):
            start = rng.randrange(length)
            for i in range(start, min(length, start + rng.randint(1, 40))):
                genome[i] = "N"
        for line in sts_lines:
            sid, p1, p2, size = line.split("\t")[:4]
            size = int(size)
            if size < len(p1) + len(p2):
                size = len(p1) + len(p2) if size > 0 else 240
            for _ in range(rng.randint(0, 3)):
                actual = max(len(p1) + len(p2), size + rng.randint(-70, 70))
                if actual + 1 >= length:
                    continue
                pos = rng.randrange(0, length - actual)
                left = instantiate(rng, p1)
                right = rc(instantiate(rng, p2))
                left = mutate(rng, left, rng.choice([0, 0, 0, 1, 2]))
                right = mutate(rng, right, rng.choice([0, 0, 0, 1, 2]))
                if rng.random() < 0.5:
                    # reverse orientation: primer2 forward, primer1 reverse complement
                    left, right = instantiate(rng, p2), rc(instantiate(rng, p1))
                genome[pos : pos + len(left)] = left
                genome[pos + actual - len(right) : pos + actual] = right
        text = "".join(genome)
        if rng.random() < 0.3:
            text = text.lower()
        width = rng.choice([60, 70, 80, 0])
        if width:
            text = "\n".join(text[i : i + width] for i in range(0, len(text), width))
        records.append(f">seq{s} description {s}\n{text}\n")
    sts_path = os.path.join(workdir, "case.sts")
    fa_path = os.path.join(workdir, "case.fa")
    with open(sts_path, "w") as fh:
        fh.write("\n".join(sts_lines) + "\n")
    with open(fa_path, "w") as fh:
        fh.write("".join(records))
    return sts_path, fa_path


def random_params(rng, iupac):
    params = []
    if rng.random() < 0.7:
        params.append(f"M={rng.choice([0, 5, 20, 50, 100, 300])}")
    if rng.random() < 0.6:
        params.append(f"N={rng.randint(0, 3)}")
    if rng.random() < 0.5:
        params.append(f"W={rng.randint(7, 12)}")
    if rng.random() < 0.4:
        params.append(f"X={rng.randint(0, 3)}")
    if rng.random() < 0.6:
        params.append(f"T={rng.randint(1, 6)}")
    if rng.random() < 0.3:
        params.append(f"Z={rng.choice([100, 240, 400])}")
    params.append(f"I={1 if iupac else 0}")
    return params


def run_mepcr(binary, sts, fa, params):
    res = subprocess.run([binary, sts, fa] + params, capture_output=True, text=True)
    return res.stdout, res.returncode, res.stderr


def run_merpcr(sts, fa, params, pure_python):
    code = (
        "import sys\n"
        "import merpcr.core.engine as e\n"
        f"e.HAS_CYTHON = {not pure_python}\n"
        "from merpcr.cli import main\n"
        "sys.argv = ['merpcr'] + sys.argv[1:]\n"
        "sys.exit(main())\n"
    )
    res = subprocess.run(
        [sys.executable, "-c", code, sts, fa] + params, capture_output=True, text=True
    )
    return res.stdout, res.returncode, res.stderr


def main():
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    ap.add_argument("--mepcr", required=True, help="path to the me-PCR executable")
    ap.add_argument("-n", "--cases", type=int, default=30)
    ap.add_argument("--seed", type=int, default=1)
    ap.add_argument(
        "--pure-python", action="store_true", help="disable the Cython extension in merPCR"
    )
    ap.add_argument("--keep", default=None, help="directory in which to keep failing cases")
    args = ap.parse_args()

    rng = random.Random(args.seed)
    failures = 0
    total_hits = 0
    for case in range(args.cases):
        iupac = rng.random() < 0.3
        workdir = tempfile.mkdtemp(prefix="merpcr_cmp_")
        sts, fa = make_case(rng, workdir, iupac)
        params = random_params(rng, iupac)
        me_out, me_rc, me_err = run_mepcr(args.mepcr, sts, fa, params)
        mer_out, mer_rc, mer_err = run_merpcr(sts, fa, params, args.pure_python)
        nhits = me_out.count("\n")
        total_hits += nhits
        if me_out == mer_out and me_rc == mer_rc:
            print(f"case {case:3d}: OK   hits={nhits:4d} params={' '.join(params)}")
            for name in os.listdir(workdir):
                os.unlink(os.path.join(workdir, name))
            os.rmdir(workdir)
        else:
            failures += 1
            print(
                f"case {case:3d}: DIFF hits={nhits:4d} params={' '.join(params)} files in {workdir}"
            )
            me_lines = me_out.splitlines()
            mer_lines = mer_out.splitlines()
            for line in sorted(set(me_lines) - set(mer_lines))[:5]:
                print(f"    only me-PCR : {line}")
            for line in sorted(set(mer_lines) - set(me_lines))[:5]:
                print(f"    only merPCR : {line}")
            if me_lines and mer_lines and set(me_lines) == set(mer_lines):
                print("    (same hits, different order)")
            if me_rc != mer_rc:
                print(f"    exit codes: me-PCR {me_rc}, merPCR {mer_rc}")
                print(f"    me-PCR stderr: {me_err.strip()[:200]}")
                print(f"    merPCR stderr: {mer_err.strip()[-200:]}")
    print(f"\n{args.cases - failures}/{args.cases} cases identical, {total_hits} hits compared")
    return 1 if failures else 0


if __name__ == "__main__":
    sys.exit(main())
