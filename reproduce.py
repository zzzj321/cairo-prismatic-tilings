"""Reproduce the eight verification examples in a fresh output directory."""

import argparse
import json
from pathlib import Path
import shutil
import subprocess
import sys
from tempfile import TemporaryDirectory


STEPS = (
    "check_groups_round32.py",
    "develop_matrix_round36.py",
    "develop_words_round36.py",
    "dualize_all_round36.py",
    "check_physical_groups_round36.py",
)


def read_json(path):
    return json.loads(path.read_text(encoding="utf-8"))


def check_results(source, work):
    witnesses = read_json(source / "witness-inputs.json")["witnesses"]
    discrete = read_json(work / "group-checks-round32.json")
    report = read_json(work / "physical-groups-round36.json")
    geometric = {case["physical_group"]: case for case in report["results"]}
    if len(report["results"]) != len(witnesses) or set(geometric) != set(witnesses):
        raise RuntimeError("The geometric check did not recover the eight expected groups.")

    for group, witness in witnesses.items():
        if "word" in witness:
            case = discrete["words"][group]
            word = "".join(map(str, witness["word"]))
            coordinates = read_json(work / f"word-{word}-round36.json")
            values_match = case["word"] == coordinates["word"] == witness["word"]
        else:
            case = discrete["matrices"][group]
            coordinates = read_json(work / f"matrix-{group}-round36.json")
            values_match = case["input"] == coordinates["matrix"] == witness["matrix"]
        physical = geometric[group]
        index = physical["translation_coset_count"]
        counts = witness["primitive_pentagons"]
        if not values_match or not physical["matches_expected"] or any(
            case[f"primitive_{kind}"] != counts[kind]
            or physical["vertex_types"][kind] != index * counts[kind]
            for kind in ("C", "P")
        ):
            raise RuntimeError(f"The input or primitive-cell counts differ for {group}.")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=Path("reproduced"),
                        help="fresh directory for generated JSON files (default: reproduced)")
    target = parser.parse_args().output.resolve()
    if target.exists():
        parser.error(f"the output directory already exists: {target}")

    source = Path(__file__).resolve().parent
    target.parent.mkdir(parents=True, exist_ok=True)
    with TemporaryDirectory(prefix="cairo-prismatic-", dir=target.parent) as temporary:
        work = Path(temporary)
        for name in STEPS:
            shutil.copyfile(source / name, work / name)
        for name in STEPS:
            subprocess.run([sys.executable, "-E", "-B", str(work / name)],
                           cwd=work, check=True)
        check_results(source, work)
        results = work / "results"
        results.mkdir()
        for path in work.glob("*.json"):
            path.replace(results / path.name)
        results.rename(target)

    print(f"Verified eight groups and primitive-cell counts. Outputs: {target}")


if __name__ == "__main__":
    main()
