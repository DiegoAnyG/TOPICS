# Objectives

## Purpose

Generate, evaluate and compare POI–PROTAC–E3 ternary-complex models with reproducible,
traceable outputs. The initial scientific scope is heterobifunctional PROTACs.

## Initial implementation sequence

1. Build a small independent Python CLI with validated inputs and a run manifest.
2. Implement molecular preparation and linker assembly with explicit attachment atoms,
   preserved stereochemistry and an ensemble of conformers with recorded seeds.
3. Generate ternary candidates from binary complexes using attachment restraints and
   clash assessment. Select and validate the modeling backend before promising capabilities.
4. Distinguish POI–PROTAC, E3–PROTAC and POI–E3 interactions with full chain/residue identities.
5. Cluster and compare candidates with explicit metrics; export CSV, JSON and PyMOL inputs.
6. Evaluate geometry relative to ubiquitination machinery when suitable models are available.

## First milestone acceptance

- Invalid structures, molecular identities and attachment mappings produce actionable errors.
- Run manifests record source provenance, hashes, configuration, seeds and software versions.
- Conformer ensembles retain molecule identity and stereochemistry; failures are explicit.
- Focused tests verify validation, repeatability and cache invalidation where implemented.
- No synthetic result is presented as a validated scientific prediction.

## Later stages

Graphical interface, larger libraries, GPU execution, automatic exit-vector selection,
quantitative cooperativity models and MD topology generation require separate implementation
and validation. A benchmark must reserve experimental structures for evaluation and report
sampling variability. Performance claims require measurements on stated hardware.

## Boundaries

Geometric accessibility does not prove ubiquitination or cellular degradation. Cooperativity
requires a defined thermodynamic model or compatible measured affinities. Do not assume a
universal lysine-distance cutoff or derive entropy from conformer count alone.

Windows and Linux are intended deployment candidates inherited from the user's context;
packaging design and supported platform guarantees remain to be established for TOPICS.

## Implementation status — 2026-09-11

The CLI, provenance, full-PROTAC conformer generation, rigid assembly baseline, proximity
contact counts and reports are implemented. Fragment/linker construction, typed interactions,
clustering and full ubiquitination geometry remain future work. Measured control errors are
large; the implemented baseline must not be described as an accurate ternary predictor.
