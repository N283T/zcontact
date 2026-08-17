# zcontact

`zcontact` is a fast, standalone Zig CLI for finding atomic and residue
contacts in biomolecular structures. It reads PDB and mmCIF, writes
deterministic TSV, processes whole structure directories in parallel, and has
no runtime dependencies.

## Build

The tested toolchain is Zig 0.16.0. Zig's pre-1.0 standard-library APIs are
not guaranteed to remain source-compatible with later releases.

```sh
zig build
zig build test --summary all
zig build bench
./zig-out/bin/zcontact --help
```

## Quick start

```sh
# Residue contacts among polymer heavy atoms (the defaults)
zcontact structure.cif > contacts.tsv

# All atom pairs no farther apart than 3.5 Å
zcontact --mode atom --cutoff 3.5 structure.pdb

# Unordered contacts across two chains
zcontact --select1 chain:A --select2 chain:B structure.cif

# Protein C-alpha residue contact map
zcontact --select protein,ca --cutoff 8.0 structure.pdb

# Deterministic, atomic per-structure outputs for a directory
zcontact batch structures/ -o contacts/ -j 8
```

The complete CLI is shown by `zcontact --help`. `-o PATH` writes the TSV to a
file instead of standard output.

`zig build bench` runs the ReleaseFast deterministic cell-list microbenchmark
at 1,000, 10,000, and 100,000 synthetic atoms. Treat it as a regression aid,
not as a published cross-tool benchmark.

For development profiling without expanding the stable CLI, use the separate
ReleaseFast real-file profiler:

```sh
zig build profile -- --iterations 9 structure.cif.gz > profile.tsv
zig build profile -- --mode atom --iterations 9 structure.pdb > atom-profile.tsv
```

It performs one unmeasured warmup per file, then reports per-iteration TSV rows
for read/decompress, raw parse, altloc resolution, residue assignment,
selection, grid construction, search/aggregation, result sorting, and TSV
formatting. It also records atom, cell, candidate-pair, accepted-pair, result,
and byte counts. This is a development harness, not a supported `zcontact`
subcommand or a cross-machine benchmark.

## Scientific definition

- An **atom contact** is an unordered pair of selected atoms whose Cartesian
  Euclidean distance is `<= cutoff`. Coordinates and cutoffs are in Å.
- A **residue contact** exists when at least one selected atom pair between the
  two residues is an atom contact. It is reported once, with the minimum atom
  distance and the atom names that realize that minimum.
- The default cutoff is 4.0 Å, the default selection is `polymer,heavy`, and
  the default scope is `inter-residue`. Thus hydrogens, `HETATM` records, and
  pairs within one residue are excluded by default.
- `--scope all` includes intra-residue pairs. An atom is never paired with
  itself.
- `--select1`/`--select2` define an unordered bipartite contact: a pair is kept
  if either orientation matches. Output side 1 matches `--select1` and side 2
  matches `--select2`; when both orientations match, input-file order wins.

These are geometric contacts, not covalent bonds, hydrogen bonds, van der
Waals overlaps, crystallographic symmetry contacts, or biological-assembly
contacts.

## Selections

A selection is a comma-separated conjunction (logical AND):

| Clause | Meaning |
| --- | --- |
| `all` | all `ATOM` and `HETATM` records |
| `polymer` | PDB/mmCIF `ATOM` records |
| `protein` | `ATOM` records with a standard/ambiguous amino-acid name |
| `hetero` | `HETATM` records |
| `ligand` | `HETATM` records except common water names |
| `water` | HOH, WAT, H2O, or DOD residue name |
| `heavy` | elements other than H and D |
| `backbone` | atom name N, CA, C, or O |
| `sidechain` | atom names other than N, CA, C, and O |
| `ca` | atom name CA |
| `chain:A` | exact author chain/asym ID |
| `name:CG` | atom name (ASCII case-insensitive) |
| `element:N` | element symbol (ASCII case-insensitive) |
| `resname:ATP` | residue/component name (ASCII case-insensitive) |
| `resseq:42` | exact signed author residue number |
| `icode:A` | exact insertion code (`_` means blank) |

`chain:_` selects a blank author chain identifier.

For example, `protein,heavy,chain:A` means standard amino-acid heavy atoms in
chain A. `protein` recognizes the 20 standard names plus ASX, GLX, SEC, and PYL;
use `polymer` for nucleic acids and other polymer `ATOM` records, or an explicit
`all,resname:...` selection for a modified residue stored as `HETATM`. Without a
CCD, record type and this explicit name list are more reproducible than guessing
chemistry.
The MVP deliberately has no boolean expression language; use the two selection
options for the common between-group case.

## Structure interpretation

- PDB identifiers use chain ID, signed residue sequence number, insertion code,
  residue name, and atom name from their standard fixed columns.
- mmCIF output identifiers prefer `auth_*` columns and fall back to `label_*`.
  Coordinates and element symbols come from `_atom_site`.
- One model is processed per invocation (`--model 1` by default). A missing
  requested model is an error.
- Alternate records for the same `(model, chain, residue number, insertion
  code)` position are resolved to one residue-coherent non-blank conformer.
  The label with the highest summed atom occupancy wins; ties prefer `A`, then
  lexical order. Blank-altloc atoms are shared, and blank wins if the same atom
  site is also present in the selected conformer. This also resolves
  microheterogeneous residue names without creating A/B chimeras.
- Residues with insertion codes remain distinct. `ATOM` and `HETATM` are both
  parsed, even when the active selection later excludes one class.
- Internal residue identity uses PDB `TER` segments or mmCIF `label_asym_id`,
  preventing repeated/blank author chain identifiers from merging molecules.
  Display and `chain:` selection continue to use author identifiers.

## Output and exit behavior

TSV always has a header. Atom mode includes stable input-order atom indices,
both complete atom identifiers (including insertion code and selected altloc),
PDB serials or mmCIF atom IDs (synthesized only when a mmCIF ID is not a valid
`u32`), and `distance_A`; residue mode includes unique internal
`residue_index1/2`, both displayed residue identifiers, `min_distance_A`, and
the closest atom names/altlocs. Rows follow deterministic input order and
distances are printed to 0.001 Å.

Exit status is 0 on success, 2 for CLI/selection errors, and 1 for input,
parsing, or output failures. Diagnostics go to standard error.

## Batch processing

```sh
zcontact batch INPUT_DIR --output-dir OUTPUT_DIR \
  --mode residue --cutoff 4.0 --select polymer,heavy --threads 10
```

Batch discovery is non-recursive and accepts recognized PDB/mmCIF suffixes.
Inputs are sorted bytewise before scheduling. Each input produces
`FILENAME.contacts.tsv`; workers may finish in any order, but contact contents
and manifest row order are deterministic. Manifest latency values naturally
vary between runs. `manifest.tsv` records status, atom/residue/contact counts,
input/output bytes, latency, and errors. `run.tsv` fixes the scientific
configuration and schema for provenance. The default worker count is the
smaller of four and the detected logical CPU count; `-j` explicitly overrides
it up to the detected CPU count.

Outputs are written to a same-directory temporary file and atomically renamed.
Existing outputs are refused by default. Use `--overwrite` to replace them or
`--resume` to reuse them. Resume requires an exact `run.tsv` configuration
match and validates each output against its deterministic
`FILENAME.contacts.tsv.meta.tsv` sidecar (input/output SHA-256, byte counts,
scientific options, and result counts; the input hash covers decompressed
structure bytes). Missing, stale, truncated, or modified
outputs are recomputed rather than silently skipped. Per-file failures are
recorded and do not stop other files, but make the final exit status 1.

## Comparison with Gemmi

[Gemmi](https://gemmi.readthedocs.io/) is a useful independent oracle for the
MVP's default atom-contact semantics:

```sh
zcontact --mode atom --cutoff 4.0 structure.cif > zcontact.tsv
gemmi contact -d 4 --ignore=1 --nosym --noh --noligand structure.cif > gemmi.txt
```

The flags align the key policies: 4 Å, ignore same-residue pairs, asymmetric
unit only, no H/D, and no ligands/water. Formatting and altloc reporting differ,
so compare normalized identifiers and distances rather than raw lines.

`scripts/validate_pairs.py STRUCTURE` performs a stronger pair-level check. It
exports zcontact's selected, conformer-resolved atoms to a surrogate mmCIF,
runs Gemmi on exactly those atoms, and compares unordered atom-index pair sets.
This cleanly separates contact-engine validation from conformer-policy
differences. Tests on 1CRN, insertion-code 1IGT, multi-model 1D3Z, altloc 1AKE,
and microheterogeneous 3NIR all have zero missing or extra pairs. The simpler
raw count check remains available as `scripts/compare-gemmi-counts.sh`.

The full AlphaFold DB v6 *E. coli* proteome (4,370 PDB plus 4,370 mmCIF files,
10,520,167 atoms) completes with zero failures. Corresponding contact TSV
contents are byte-identical between PDB and mmCIF, and between 1-thread and
10-thread runs; filenames and provenance sidecars intentionally differ. See
[`docs/VALIDATION.md`](docs/VALIDATION.md) for commands and measurements.

## MVP limitations

- Plain or gzip-compressed PDB/mmCIF (`.pdb[.gz]`, `.ent[.gz]`,
  `.cif[.gz]`, `.mmcif[.gz]`). Compression is detected from magic bytes.
- Asymmetric-unit Cartesian coordinates only; no unit-cell symmetry or assembly
  expansion.
- One model per run; models are not compared with each other.
- A deterministic uniform cell list limits distance checks to neighboring
  cutoff-sized cells. Memory use is O(n + contacts); input and gzip
  decompression are whole-file and each worker materializes its contact result,
  so explicitly high `-j` values multiply peak memory. Decompressed input is
  bounded at 1 GiB per file. Periodic boxes are not supported.
- No MD trajectory input and no dependency on `ztraj` in this static-structure
  MVP.

See [DESIGN.md](DESIGN.md) for investigation notes and deferred decisions.
