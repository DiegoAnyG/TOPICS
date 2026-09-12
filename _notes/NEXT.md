# Next session — start here

Updated 2026-09-12. This is a handover, not a historical log.

## State

- TOPICS 0.1 has interactive/batch CLI, conformer-driven assembly, explicit atom mappings,
  provenance, optional CUDA distances, reference evaluation and offline reports.
- uv.lock fixes the environment. The local .venv is installed with dev and GPU extras.
  Start with `uv run topics`; recreate with `uv sync --extra dev --extra gpu` when needed.
- Two retrospective controls (5T35 and 6BOY), 32 conformers per seed and seeds 42/43/44,
  expose large prediction errors. Read docs/VALIDATION.md before making accuracy claims.
- Geometry/mapping positive controls pass; they do not establish predictive accuracy.
- The core methods and input contract are in docs/METHODS.md.
- The user requested short English commits and changelog items, pushed to the configured origin.

## Accuracy research — 2026-09-12

- Read docs/ACCURACY_RESEARCH.md and docs/ACCURACY_PLAN.md for the current priorities.
- The six-pose independent audit is in docs/accuracy-audit.json, reproduced by docs/audit_accuracy.py.
  Both rank-1 PPI DockQ scores are poor; all 96 candidates per control have current clash pairs.
- Complete-head binding geometry is lost despite small ring-fit RMSD. The lexicographic
  ranking is dominated by clash count. No scientific engine changes were made during research.
- The research report includes 24 references and Markdown, offline HTML and PDF versions.
  Optional DockQ was used in an isolated temporary environment, not added to uv.lock.

## Next work

1. Complete P0: general independent evaluation, chemical validity and a curated benchmark
   manifest with protected family/scaffold groups, input strata and fixed budgets.
2. Complete P1: explicit complete-head/linker mappings, preserved binary binding geometry,
   SDF/SMILES input and constrained linker/protein-pose sampling; retain schema-1 compatibility.
3. Evaluate coarse PPI search, selective refinement, clustering and physical ranking with
   development/validation ablations before a frozen test evaluation.
4. Pilot DeepTernary in isolation and use Boltz-2 as an optional comparator after memory
   checks. Published metrics differ by input protocol and best-of-N versus rank-1 selection.
5. Validate Windows installation/runtime; only Linux/WSL and local CUDA were tested here.

Never start Docker Desktop automatically. Do not reuse an output directory for a new run.
