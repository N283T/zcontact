# MVP design notes

## Investigation

The repository initially contained only `README.md` and `AGENTS.md`. The first
session inspected the local Zig 0.16 projects `zsasa`, `zdssp`, `zreduce`, and
`ztraj` under `/Users/nagaet/ghq/github.com/N283T`:

- `zsasa`, `zdssp`, and `zreduce` demonstrate self-contained PDB/mmCIF parsing,
  explicit model/altloc policies, a library module plus thin CLI, Zig 0.16 build
  layout, and unit plus real-file validation.
- `ztraj` already includes richer trajectory selections, neighbor lists, and
  trajectory contacts. Importing it would couple a static contact CLI to a much
  larger trajectory/SASA dependency graph, so this MVP remains independent.
- All relevant sibling code is MIT-licensed. The implementation here was kept
  small and written locally rather than importing sibling source as a library.

The local data tree contains AlphaFold DB proteomes in PDB/mmCIF form and a
divided experimental PDB mmCIF mirror. Clean AlphaFold structures were used for
smoke tests; experimental `1crn` was checked in both mmCIF and converted-PDB
forms. Local datasets were not changed.

## Decisions

1. **Minimum-distance geometric contacts first.** This is unambiguous, easy to
   compare with Gemmi, and supports both atom lists and residue contact maps.
   Chemistry-aware definitions can be separate future modes rather than hidden
   changes to `--cutoff`.
2. **One inclusive cutoff.** `distance <= cutoff`, in Å. Squared distances are
   compared internally and square roots are taken only for reported contacts.
3. **Small conjunction selections.** The MVP covers protein/heavy/backbone/CA
   and exact identifier filters. A full parser with NOT/OR/parentheses is
   intentionally deferred until concrete use cases require it.
   The default is `polymer,heavy`: `polymer` has the format-defined `ATOM`
   meaning, while `protein` uses an explicit standard amino-acid name list.
4. **Author identifiers for mmCIF display.** These compare naturally with
   conventional PDB output. Missing/unknown author values fall back per row to
   label identifiers; internal residue identity retains `label_asym_id` so
   repeated or blank author IDs do not merge distinct molecules.
5. **Residue-coherent altloc collapse.** Per-atom highest occupancy creates
   physically impossible A/B chimeras and mishandles microheterogeneity. One
   non-blank label is therefore selected per internal residue position by
   summed occupancy; ties prefer A then lexical order. Blank atoms are shared.
6. **Single-model invocation.** Pooling multiple models would create false
   cross-model contacts. Per-model batch/ensemble aggregation needs a separate
   output definition.
7. **TSV only.** It is stable, inspectable, and script-friendly. JSON and sparse
   matrix formats are deferred until consumers establish the schemas they need.

## Validation policy

- Unit tests cover selection conjunctions, inclusive geometry, cell-list versus
  brute force, residue minimum aggregation, PDB model filtering, coherent
  altloc choice, TER segmentation, auth/label fallback, gzip bounds, and batch
  atomic/resume behavior.
- `tests/mini.pdb` and `tests/mini.cif` are redistributable synthetic fixtures
  containing altloc, model, `HETATM`, and author-vs-label identifier cases.
- Gemmi `contact` is the current external comparison oracle. Comparisons must
  explicitly align hydrogen, ligand, same-residue, symmetry, model, and cutoff
  policies before interpreting differences.
- `scripts/validate_pairs.py` compares exact unordered pairs against Gemmi after
  normalizing the selected atom inventory. It covers clean, insertion-code,
  multi-model, altloc, and microheterogeneous experimental structures.
- Full AFDB *E. coli* PDB/mmCIF batch and cross-thread SHA-256 equality are
  recorded in `docs/VALIDATION.md`.

## Deferred scope

The first-session follow-up replaced the quadratic scan with a deterministic
cutoff-sized cell list, added bounded gzip input, and added `zig build bench`.
The synthetic benchmark scales from 0.72 ms at 1,000 atoms through 59 ms at
100,000 atoms on the initial development machine; these are local regression
numbers, not cross-tool performance claims. ReleaseFast smoke tests on 4,766-
and 10,882-atom AlphaFold models also complete below the 0.01 s resolution of
`/usr/bin/time` on that machine.

Assembly/symmetry expansion, periodic boxes, chemistry-aware contact types,
JSON/matrix schemas, recursive/multi-root batch discovery, and trajectory I/O
remain separate features. They should be added only with concrete scientific
definitions and fixtures; trajectory support may narrowly reuse `ztraj`.
