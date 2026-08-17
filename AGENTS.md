# Project Instructions

## Project scope

- This repository is an independent member of the ZBioKit family.
- Its initial scope is **atomic and residue contact analysis for biomolecular structures**.
- Keep the command-line tool useful on its own. The unified ZBioKit distribution is out of scope for now.
- Prefer a small, testable MVP and evolve the interface from concrete use cases rather than building a large framework first.

## Data and validation

- Local structural-biology test data is available under `/Users/nagaet/pdb`.
- AlphaFold DB structures under `/Users/nagaet/pdb/afdb` are generally clean and are good fixtures for initial correctness and smoke tests.
- Also test representative experimental PDB/mmCIF structures when behavior depends on alternate locations, insertion codes, missing atoms, multiple models, assemblies, ligands, or other real-world edge cases.
- Do not modify source datasets under `/Users/nagaet/pdb`. Keep small redistributable fixtures in this repository when tests require stable inputs.
- Established third-party tools may be used as comparison oracles and for benchmarking. Use already-installed tools or project-local/declarative environments; do not install packages globally.

## Relationship to sibling projects

- Useful sibling repositories live under `/Users/nagaet/ghq/github.com/N283T`, especially `zsasa`, `zdssp`, `zreduce`, and `ztraj`.
- Inspect them for Zig project layout, parsing approaches, performance techniques, testing patterns, and optimization ideas where useful.
- Keep this repository self-contained: do not import sibling projects as implementation libraries merely to share convenience code.
- A deliberate exception is molecular-dynamics file I/O: reusing or depending on `ztraj` may be considered when it materially avoids duplicating format support. Keep that dependency narrow and explicit.
- Copying a small, well-understood implementation into this repository is preferable to introducing tight cross-repository coupling, but preserve applicable licenses and attribution.

## Engineering priorities

- Correctness and clearly documented scientific definitions come before micro-optimization.
- Make selection rules, units, cutoffs, atom correspondence, and edge-case behavior explicit in the CLI and documentation.
- Add tests before or alongside optimization. Benchmark against established tools where practical.
- Favor deterministic, script-friendly output and stable exit behavior.
