# Changelog

Notable user-facing and project-level changes are recorded here. Detailed
scientific validation and local performance measurements remain in
[`docs/VALIDATION.md`](docs/VALIDATION.md).

## [Unreleased]

### Added

- Public Ubuntu and macOS CI for formatting, tests, ReleaseFast builds,
  version checks, benchmark smoke tests, and profiler smoke tests.
- SHA-256 verification of official Zig compiler archives used by CI.

## [0.1.1] - 2026-08-17

### Added

- A development-only real-file stage profiler available through
  `zig build profile -- ...`.
- Deterministic randomized differential tests against a scalar all-pairs
  reference across atom/residue modes, selections, scopes, and cutoffs.
- Parser error-path coverage for malformed mmCIF input.

### Changed

- Reduced conformer-resolution overhead with stable in-place blank-altloc
  handling and compact residue-index hash keys.
- Reduced contact-search overhead by traversing occupied cell pairs instead of
  repeating neighbor-cell lookups for every atom.
- Reduced mmCIF atom-loop parsing and atom-mode TSV formatting overhead.
- Reduced peak memory in parallel batch processing while preserving output.

## [0.1.0] - 2026-08-17

### Added

- Standalone Zig CLI for atom and residue contact analysis.
- PDB, mmCIF, and bounded gzip input with explicit model selection and
  residue-coherent alternate-conformer handling.
- Inclusive distance cutoffs, symmetric and asymmetric selections, and
  deterministic TSV output.
- Cutoff-sized cell-list search and residue minimum-distance aggregation.
- Parallel non-recursive batch processing with atomic outputs, manifests,
  SHA-256 sidecars, resume validation, and overwrite protection.
- Unit, integration, CLI, Gemmi pair-set, and AFDB corpus validation.

[Unreleased]: https://github.com/N283T/zcontact/compare/v0.1.1...HEAD
[0.1.1]: https://github.com/N283T/zcontact/compare/v0.1.0...v0.1.1
[0.1.0]: https://github.com/N283T/zcontact/tree/v0.1.0
