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

## Profiling-first optimization record

The `zig build profile -- ...` development harness was added after the v0.1.0
baseline was tagged. It is compiled with `ReleaseFast`, performs one unmeasured
warmup per input, writes measured rows as TSV, and keeps profiling out of the
stable public CLI. Instrumentation is opt-in: ordinary parser/contact calls do
not read clocks or increment profiling counters.

The representative corpus deliberately combines small and large structures,
gzip and plain input, experimental and AFDB data, and equivalent PDB/mmCIF
representations:

- experimental `1crn.cif.gz` (327 selected atoms),
- experimental `1igt.cif.gz` (12,956 parsed atoms; insertion-code coverage),
- AFDB `AF-P76347-F1-model_v6.cif` and its corresponding PDB (17,331 atoms).

Initial stage measurements showed that parse consumed 36--52% of end-to-end
time on the larger files. Splitting parse revealed that altloc resolution alone
used 23--34% despite the AFDB structures having no alternate conformers.
Contact search/aggregation used another 28--43%; atom-mode TSV formatting used
15--23%. Selection, grid construction, and final sorting were individually
small. These results did not support distance-kernel SIMD as the first change.

Two measured structural optimizations followed:

1. Blank-altloc structures now use one stable in-place site-deduplication pass.
   The conformer-aware path preserves first-site order with placeholders, which
   removes a second site map and the final sort. Altloc-stage median time fell
   56--73% on the four representative inputs.
2. Contact traversal now visits pairs of occupied cells rather than performing
   27 neighbor hash lookups for every atom. On the 17,331-atom AFDB structure,
   neighbor lookups consequently scale with 5,626 occupied cells rather than
   all atoms. Candidate and accepted-pair counts are unchanged. Search-stage
   median time fell 38--47% across the corpus.

Nine-iteration warm-cache medians after both changes, compared with the
instrumented pre-optimization baseline, were:

| Mode | Structure | Before | After | Change |
| --- | --- | ---: | ---: | ---: |
| residue | 1CRN gzip mmCIF | 1.363 ms | 1.013 ms | -25.7% |
| residue | 1IGT gzip mmCIF | 19.065 ms | 13.022 ms | -31.7% |
| residue | AF-P76347 mmCIF | 23.844 ms | 15.878 ms | -33.4% |
| residue | AF-P76347 PDB | 20.576 ms | 12.959 ms | -37.0% |
| atom | 1IGT gzip mmCIF | 21.098 ms | 16.227 ms | -23.1% |
| atom | AF-P76347 mmCIF | 28.219 ms | 21.093 ms | -25.3% |
| atom | AF-P76347 PDB | 24.887 ms | 17.984 ms | -27.7% |

The sub-millisecond 1CRN atom-mode result was noise-sensitive and is omitted
from the improvement table. These are local regression measurements, not
cross-tool claims.

Correctness gates after optimization were stronger than count comparison:

- 10 default atom/residue outputs over 1CRN, 1IGT, 1D3Z, and both AF-P76347
  formats were byte-identical to the tagged v0.1.0 binary;
- an asymmetric chain-selection output and selected-atom inventory were also
  byte-identical to v0.1.0;
- the five Gemmi oracle structures retained zero missing and zero extra pairs;
- the complete 20-test unit/integration/CLI suite passed.

The final profiles still show mixed costs: parse is 29--53% on the larger
cases, search/aggregation is 17--37%, and atom TSV formatting is 21--33%.
There is therefore still no evidence that explicit SIMD of the scalar distance
expression should take priority over parser or output work.

### Second profiling pass

A second pass targeted those remaining mixed costs without changing the public
CLI or output schema:

1. The mmCIF parser now starts tokenization from a validated lexical hint near
   the atom-site loop, while retaining position zero as the complete grammar
   fallback. Atom rows use a fixed-size slice array after the tag count is
   known rather than repeatedly appending to a dynamic row list. This reduced
   raw mmCIF parsing by 5--8% on the larger files.
2. Residue indices are assigned before conformer resolution. Altloc scoring and
   atom-site deduplication can consequently hash compact integer residue IDs
   instead of repeated model/chain/sequence/insertion structures. This reduced
   the altloc stage by 54--57% on the large AFDB PDB/mmCIF pair and 27% on the
   alternate-containing 1IGT structure.
3. Atom TSV rows use one formatter invocation instead of three. Output-stage
   time fell about 6% on the larger atom-mode cases. A subsequent manual stack
   line-builder experiment made output 4--7% slower and was rejected.

Nine-to-eleven-iteration warm-cache medians relative to commit `40e8991` were:

| Mode | Structure | Before | After | Change |
| --- | --- | ---: | ---: | ---: |
| residue | 1CRN gzip mmCIF | 1.049 ms | 0.878 ms | -16.3% |
| residue | 1IGT gzip mmCIF | 13.306 ms | 12.857 ms | -3.4% |
| residue | AF-P76347 mmCIF | 15.832 ms | 13.736 ms | -13.2% |
| residue | AF-P76347 PDB | 12.410 ms | 10.778 ms | -13.2% |
| atom | 1IGT gzip mmCIF | 16.025 ms | 15.451 ms | -3.6% |
| atom | AF-P76347 mmCIF | 20.987 ms | 18.846 ms | -10.2% |
| atom | AF-P76347 PDB | 17.907 ms | 15.878 ms | -11.3% |

Fourteen default atom/residue outputs across the five Gemmi structures and the
AF-P76347 PDB/mmCIF pair, plus asymmetric selection and atom inventory output,
were byte-identical to `40e8991`. The 20-test suite and all five Gemmi pair-set
comparisons also remained exact.

### Full-corpus post-optimization batch

The complete AFDB *E. coli* corpus was rerun after commit `a285944`, using the
same ReleaseFast binary and known PDB/mmCIF directories. All three runs produced
4,370 successes, zero failures, 10,520,167 atoms, 1,349,634 residues, 5,505,835
contacts, and 251,623,549 bytes of contact TSV output.

| Input | Workers | Wall | Peak RSS | Success | Failure |
| --- | ---: | ---: | ---: | ---: | ---: |
| PDB | 1 | 25.60 s | 15.5 MB | 4,370 | 0 |
| PDB | 10 | 4.20 s | 50.3 MB | 4,370 | 0 |
| mmCIF | 10 | 6.07 s | 48.9 MB | 4,370 | 0 |

These are warm-cache local regression runs rather than simultaneous controlled
A/B measurements, so wall-time differences from the earlier record are only
indicative. Peak RSS nevertheless fell from the previously recorded 151--169
MB at ten workers to about 49--50 MB; compact parser hash keys materially reduce
per-worker transient memory.

SHA-256 was recomputed for every contact TSV. All 4,370 hashes matched between
the one- and ten-worker PDB runs, and all 4,370 matched accession-by-accession
between PDB and mmCIF. This confirms that the profiling-guided parser and search
changes retain both scheduling and format determinism at full-corpus scale.
