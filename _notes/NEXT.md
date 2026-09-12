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

## Next work

1. Design a larger held-out benchmark before changing sampling or ranking using these controls.
2. Improve physically justified sampling and interface refinement; the current rigid baseline
   has no protein flexibility, PPI energy or calibrated scoring.
3. Consider a chemically explicit fragment/linker builder and general SDF input. Currently
   the program requires a complete PROTAC CCD and local head poses.
4. Add residue-typed interaction profiling, clustering and full-complex cofactor handling
   when their methods and validation are defined.
5. Validate Windows installation/runtime; only Linux/WSL and local CUDA were tested here.

Never start Docker Desktop automatically. Do not reuse an output directory for a new run.
