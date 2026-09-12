# Next session — start here

Updated 2026-09-12. This is a handover, not a historical log.

## State

- TOPICS 0.1 retains its conformer-driven baseline, interactive/batch CLI, provenance,
  optional CUDA distances and offline scientific reports. Prediction/ranking have not changed.
- P0 adds `topics assess`: all-candidate independent DockQ, symmetry ligand RMSD, complete
  identity mappings, optional whole-head metrics, stereo/ligand geometry and component sterics.
- The optional `.venv-evaluation` is installed separately, using requirements-evaluation.lock;
  production uv.lock remains unchanged. See docs/ASSESSMENT.md for CLI/environment selection.
- 43 PDB entries were reviewed; benchmarks/inventory-v1.json freezes 32 eligible cases,
  11 exclusions and 16/10/6 development/validation/test assignments with protected groups.
- The inventory is structurally curated, not a ready-to-execute benchmark. Complete chemical
  head/linker definitions, independent input sources and process budget enforcement remain pending.
- All 192 original control candidates have independent assessments. See docs/P0_RESULTS.md
  and the separate local runs/accuracy-p0 reports. All original control files are preserved.
- Ligand geometry/stereo passes for both ensembles, but every candidate has intercomponent
  steric violations. Only one 5T35 candidate reaches acceptable PPI DockQ, at rank 24.
- Whole-head and joint endpoints remain unassessed for these ring-only baseline inputs.
  Geometry controls prove implementation consistency, not predictive accuracy.
- The user requested concise English commits and changelog items pushed to the configured origin.

## Next work

1. Begin P1: curate complete chemical heads, linker atoms and attachment bonds from chemistry;
   preserve schema-1 input compatibility and independently framed protein/head binding poses.
2. Add stereo-preserving SDF/SMILES mappings and constrained linker/protein-pose sampling.
   Compare against the unchanged baseline on development/validation under matched budgets.
3. Revisit the five CCD-preparation exclusions (QIY, TOO and YF8) without changing chemical
   identity or protected groups. Amend inventory versions explicitly if eligibility changes.
4. Prepare runnable benchmark inputs, independent binary/apo sources and a case-level runner
   that enforces wall-time/memory limits and counts preparation failures/timeouts in the denominator.
5. Follow P2/P3 interface search, selective refinement, clustering and ranking only after P1
   chemistry/head-preservation checks pass. Optional learned engines need checkpoint overlap
   and memory audits. No protected test prediction has been run; keep it that way until freeze.
6. Validate Windows installation/runtime. Current optional adapter checks ran on Linux/WSL;
   CUDA is unnecessary for assessment and its existing runtime test may skip in a sandbox.

Read docs/ACCURACY_PLAN.md, docs/ASSESSMENT.md and benchmarks/README.md before implementation.
Never start Docker Desktop automatically. Do not reuse output directories for new work.
