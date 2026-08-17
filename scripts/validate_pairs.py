#!/usr/bin/env python3
"""Independent pair-level validation of zcontact's geometric engine with Gemmi.

The script asks zcontact for its selected, altloc-resolved atom inventory, writes
a surrogate mmCIF with unique atom names and blank altlocs, and compares
unordered atom-index pair sets. This deliberately validates contact geometry
separately from each program's conformer policy and identifier formatting.
"""

from __future__ import annotations

import argparse
import csv
import pathlib
import shutil
import subprocess
import sys
import tempfile


def run(command: list[str], *, stdout=None) -> None:
    subprocess.run(command, check=True, stdout=stdout)


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("structure")
    parser.add_argument("--zcontact", default=None)
    parser.add_argument("--gemmi", default="gemmi")
    parser.add_argument("--cutoff", type=float, default=4.0)
    parser.add_argument("--select", default="polymer,heavy")
    parser.add_argument("--model", type=int, default=1)
    args = parser.parse_args()

    root = pathlib.Path(__file__).resolve().parent.parent
    zcontact = pathlib.Path(args.zcontact) if args.zcontact else root / "zig-out/bin/zcontact"
    if not zcontact.is_file():
        parser.error(f"{zcontact} not found; run zig build -Doptimize=ReleaseFast")
    if shutil.which(args.gemmi) is None:
        parser.error(f"Gemmi executable not found: {args.gemmi}")
    if not (args.cutoff > 0):
        parser.error("--cutoff must be > 0")

    with tempfile.TemporaryDirectory(prefix="zcontact-validate-") as directory:
        tmp = pathlib.Path(directory)
        atoms_path = tmp / "atoms.tsv"
        contacts_path = tmp / "zcontact.tsv"
        surrogate_path = tmp / "selected.cif"
        gemmi_path = tmp / "gemmi.txt"
        with contacts_path.open("w") as output:
            run(
                [
                    str(zcontact), "--mode", "atom", "--cutoff", str(args.cutoff),
                    "--select", args.select, "--model", str(args.model),
                    "--atoms-output", str(atoms_path), args.structure,
                ],
                stdout=output,
            )

        with atoms_path.open(newline="") as stream:
            atoms = list(csv.DictReader(stream, delimiter="\t"))
        if not atoms:
            raise RuntimeError("zcontact selected no atoms")
        names = {f"Z{int(atom['index']):X}": int(atom["index"]) for atom in atoms}

        with surrogate_path.open("w") as out:
            out.write(
                "data_zcontact_selected\nloop_\n"
                "_atom_site.group_PDB\n_atom_site.id\n_atom_site.type_symbol\n"
                "_atom_site.label_atom_id\n_atom_site.label_alt_id\n"
                "_atom_site.label_comp_id\n_atom_site.label_asym_id\n"
                "_atom_site.label_entity_id\n_atom_site.label_seq_id\n"
                "_atom_site.Cartn_x\n_atom_site.Cartn_y\n_atom_site.Cartn_z\n"
                "_atom_site.occupancy\n_atom_site.B_iso_or_equiv\n"
                "_atom_site.auth_seq_id\n_atom_site.auth_comp_id\n"
                "_atom_site.auth_asym_id\n_atom_site.auth_atom_id\n"
                "_atom_site.pdbx_PDB_model_num\n"
            )
            for serial, atom in enumerate(atoms, 1):
                name = f"Z{int(atom['index']):X}"
                residue = int(atom["residue_index"]) + 1
                element = atom["element"] or "C"
                out.write(
                    f"ATOM {serial} {element} {name} . GLY A 1 {residue} "
                    f"{atom['x']} {atom['y']} {atom['z']} 1.0 0.0 "
                    f"{residue} GLY A {name} 1\n"
                )
            out.write("#\n")

        with gemmi_path.open("w") as output:
            run(
                [args.gemmi, "contact", "-d", str(args.cutoff), "--ignore=1", "--nosym", str(surrogate_path)],
                stdout=output,
            )

        z_pairs: dict[tuple[int, int], float] = {}
        with contacts_path.open(newline="") as stream:
            for row in csv.DictReader(stream, delimiter="\t"):
                pair = tuple(sorted((int(row["index1"]), int(row["index2"]))))
                z_pairs[pair] = float(row["distance_A"])

        gemmi_pairs: dict[tuple[int, int], float] = {}
        with gemmi_path.open() as stream:
            for line_number, line in enumerate(stream, 1):
                fields = line.split()
                atom_names = [field for field in fields[:-1] if field in names]
                if len(atom_names) != 2:
                    raise RuntimeError(f"cannot parse Gemmi line {line_number}: {line.rstrip()}")
                pair = tuple(sorted((names[atom_names[0]], names[atom_names[1]])))
                gemmi_pairs[pair] = float(fields[-1])

        z_set = set(z_pairs)
        gemmi_set = set(gemmi_pairs)
        missing = sorted(gemmi_set - z_set)
        extra = sorted(z_set - gemmi_set)
        shared = z_set & gemmi_set
        max_distance_delta = max((abs(z_pairs[pair] - gemmi_pairs[pair]) for pair in shared), default=0.0)
        print(f"atoms\t{len(atoms)}")
        print(f"zcontact_pairs\t{len(z_set)}")
        print(f"gemmi_pairs\t{len(gemmi_set)}")
        print(f"missing_pairs\t{len(missing)}")
        print(f"extra_pairs\t{len(extra)}")
        print(f"max_display_distance_delta_A\t{max_distance_delta:.6f}")
        if missing or extra or max_distance_delta > 0.0055:
            if missing:
                print(f"first_missing\t{missing[:5]}", file=sys.stderr)
            if extra:
                print(f"first_extra\t{extra[:5]}", file=sys.stderr)
            return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
