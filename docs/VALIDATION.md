# Validation and performance record

This record captures the initial publication-oriented validation on 2026-08-17.
Source datasets under `/Users/nagaet/pdb` were read-only.

## Environment

- Apple M4, 10 logical CPUs, 32 GiB RAM
- macOS/Darwin 25.5.0 arm64
- Zig 0.16.0, `ReleaseFast`
- Gemmi 0.7.4

## Unit and integration tests

```sh
zig fmt --check build.zig build.zig.zon src/*.zig tests/*.zig
zig build test --summary all
zig build -Doptimize=ReleaseFast
```

The 20-test suite covers PDB/mmCIF parsing, gzip bounds, model selection, insertion
codes, row-level auth/label fallback, coherent altloc selection, TER-separated
blank chains, selections, inclusive cutoffs, cell-list/brute-force equality,
residue minimum distances, deterministic parallel batch, atomic outputs, and
resume configuration rejection. Batch tests also corrupt an output and change
an input, proving that hash-validated resume recomputes stale files while
preserving counts for valid skipped files.
The same `zig build test` gate also runs a portable black-box CLI smoke script
covering stdout/stderr conventions, residue headers, asymmetric side
orientation, empty-selection rejection, preservation of an existing output on
parse failure, and batch selection-default parity.

## Independent pair-level Gemmi oracle

Raw Gemmi contact output uses different altloc semantics. The validator first
asks zcontact for the complete selected atom inventory, writes those coordinates
with blank altlocs and unique atom names to a surrogate mmCIF, then compares
unordered pair sets. This independently tests the neighbor search and distance
predicate without confusing conformer policies.

```sh
scripts/validate_pairs.py STRUCTURE.cif.gz
```

All used the default heavy-polymer selection, 4 Å inclusive cutoff, model 1,
and inter-residue scope:

| Structure | Edge case | Atoms | Pairs | Missing | Extra | Max displayed distance delta |
| --- | --- | ---: | ---: | ---: | ---: | ---: |
| 1CRN | clean experimental | 327 | 1,070 | 0 | 0 | 0.005 Å |
| 1IGT | insertion codes | 10,214 | 29,826 | 0 | 0 | 0.005 Å |
| 1D3Z | multi-model, model 1 | 602 | 1,632 | 0 | 0 | 0.005 Å |
| 1AKE | alternate locations | 3,312 | 10,276 | 0 | 0 | 0.005 Å |
| 3NIR | altloc + microheterogeneity | 327 | 1,098 | 0 | 0 | 0.005 Å |

Gemmi prints distances to 0.01 Å and zcontact to 0.001 Å, so 0.0055 Å is the
comparison tolerance. Gemmi 0.7.4 `contact --count` is not used because its
reported count was half the normal output row count in these tests.

Conformer policy is tested separately from that geometry oracle. In experimental
6Y2T, chain `AAA` residue 40 has equal-occupancy A/B sites that make a
per-atom winner rule choose a physically impossible mixture. The selected atom
inventory contained all 14 THR atoms from conformer B (and none from A),
including OG1/CG2 where a site-wise tie would otherwise prefer A. The synthetic
unit fixture additionally forces opposing per-site occupancy winners and checks
that the residue-level summed-occupancy choice remains coherent.

## Full *E. coli* AFDB v6 batch

Corpus:

- 4,370 PDB files, 878 MB
- 4,370 mmCIF files, 1.266 GB
- 10,520,167 selected atoms
- 1,349,634 residues
- 5,505,835 residue contacts

Representative command:

```sh
/usr/bin/time -l ./zig-out/bin/zcontact batch \
  /Users/nagaet/pdb/afdb/UP000000625_83333_ECOLI_v6/pdb \
  --output-dir /tmp/zcontact-ecoli-pdb --threads 10 --quiet
```

| Input | Workers | Wall | Peak RSS | Success | Failure | Output |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| PDB | 10 | 9.27 s | 151 MB | 4,370 | 0 | 251,623,549 B |
| PDB | 1 | 30.65 s | 48 MB | 4,370 | 0 | 251,623,549 B |
| mmCIF | 10 | 8.43 s | 169 MB | 4,370 | 0 | 251,623,549 B |

Two strong determinism checks passed:

1. SHA-256 for every one of 4,370 PDB result files was identical with one and
   ten workers.
2. SHA-256 matched accession-by-accession for all 4,370 corresponding PDB and
   mmCIF contact files.

Times are warm-cache, single-run development measurements after adding unique
residue indices and hash-validated resume sidecars; they should not be
presented as cross-tool benchmarks. The output byte column counts contact TSVs,
not metadata sidecars. The manifest provides per-file latency and counts for
repeated statistical benchmarking; latency itself is not deterministic.

## Synthetic scaling benchmark

`zig build bench` performs warmup plus repeated ReleaseFast measurements on a
bounded-density synthetic chain:

| Atoms | Contacts | Time/iteration |
| ---: | ---: | ---: |
| 1,000 | 1,997 | 0.78 ms |
| 10,000 | 19,997 | 7.46 ms |
| 100,000 | 199,997 | 59.7 ms |

This is a regression harness, not a biological or cross-tool benchmark.
