# TOPICS Accuracy Implementation Checklist

Proposed on 2026-09-12. P0 is implemented with the scope and limits in [ASSESSMENT.md](ASSESSMENT.md);
scientific sampling/ranking changes in P1–P4 remain pending. The rationale,
method comparison, estimates and citations are in [ACCURACY_RESEARCH.md](ACCURACY_RESEARCH.md).
Keep baseline commit `b2ae904` and its measured outputs available for comparison.

## P0 — Evaluation and protected benchmark

- [x] Add optional DockQ evaluation with explicit chain mapping, role-defined RMSDs and
  symmetry-aware ligand evaluation. Retain existing metrics under their original definitions.
- [x] Define chemical-validity checks, complete-head RMSD and component-specific clashes.
  Cases with absent native PPI contacts need explicit metric applicability and fallback.
- [x] Curate an initial inventory aiming at 20–40 eligible ternary cases. Record all
  exclusions; separate binary structures, molecular glues and incomplete ligands.
- [x] Freeze group-based development/validation/test splits, training-overlap annotations,
  input difficulty, success criteria, seed policy and equal compute budgets.
- [x] Verify native identity, rigid-transform invariance, residue/atom identity and failure
  reporting. Confirm evaluation references cannot enter generation or ranking.

Gate: the evaluator and inventory are reproducible before changing scientific scoring.
The existing six-pose DockQ audit is a pilot, not this complete adapter.

Delivered: all 192 existing candidates assessed, 43 entries reviewed, 32 eligible cases with
protected splits. Complete-head metrics accept curated definitions but remain unavailable
for the existing ring-only controls. P1 must curate those definitions and runnable benchmark
inputs; full benchmark execution also requires budget enforcement. Training overlap is unknown,
not assumed absent. Missing protein residues require curated complete mappings, not automatic alignment.

## P1 — Complete heads and linker sampling

- [ ] Extend input mappings to complete chemical heads, linker and attachment bonds;
  preserve schema-1 input compatibility as a labeled legacy path.
- [ ] Add SDF/SMILES inputs with stereo-preserving atom maps and ambiguous-match errors.
- [ ] Preserve each head relative to its own protein, allowing the two units to move
  relative to each other. Never constrain the native ternary placement.
- [ ] Compare constrained linker-driven assembly and protein-pose-driven linker closure
  with the unchanged free-conformer baseline under matched budgets.
- [ ] Record head distortion, closure failures, conformer diversity and feasible yield.

Gate: chemical and head-preservation properties pass; improvements in feasible yield and
oracle structural agreement are measured on development/validation, not asserted in tests.

## P2 — Interface search and selective refinement

- [ ] Add coarse protein-interface proposals and linker-compatible local moves.
- [ ] Use typed steric checks with cofactor/protonation handling and explicit failures.
- [ ] Refine a small diverse subset; begin with ligand constraints, then add restricted
  side-chain relaxation with a validated optional OpenMM environment.
- [ ] Cluster candidates and retain representatives from distinct feasible orientations.

Gate: paired gains persist across seeds within the budget; refinement does not destroy
binding geometry or chemistry. If feasible poses remain absent, diagnose before scaling.

## P3 — Ranking and external generators

- [ ] Replace lexicographic feature dominance with a documented feasibility gate and
  interpretable physical consensus; tune only on grouped validation data.
- [ ] Measure success@1/5/10 separately from oracle success and chemical validity.
- [ ] Pilot isolated DeepTernary with checkpoint/input provenance and native ranking.
- [ ] Benchmark Boltz-2 if memory allows; keep sequence-input and binary-input tasks separate.
- [ ] Record end-to-end timing, peak memory, failed runs and calibration. Add more engines
  only if the comparison resolves an outstanding scientific question.

Gate: ranking gains reproduce beyond the two development controls; local hardware limits
are measured. No automatic probability-of-correctness claim from confidence scores.

## P4 — Frozen evaluation and scientific outputs

- [ ] Freeze the selected pipeline and run the protected test set with every eligible case
  in the denominator, including timeouts and zero-feasible-pose outcomes.
- [ ] Report paired system/group-level uncertainty, per-family results and all ablations.
- [ ] Extend the offline report with clusters, validity, interface details and uncertainty.
- [ ] Export publication figures in SVG/PDF, case-level source data and reproducible Methods.
- [ ] Update notes and brief English changelog items; commit verified stages and push them.

Gate: communicate measured structural improvement and its limits. Degradation,
cooperativity and ubiquitination remain separately validated endpoints.

## First implementation unit

Complete P0 first: evaluate every existing candidate through the independent adapter,
define joint structural/chemical metrics, and commit the curated benchmark manifest with
protected groups. P1 follows without selecting parameters against final test structures.
Use a separate short commit for each completed unit; do not replace the baseline outputs.
