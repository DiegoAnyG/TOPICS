# Decisions

Newest first. Record each decision, its reason and its cost. Announce and document reversals.

## 2026-09-11 — Independent project

Decision: keep TOPICS outside PoliScreen, with its own AGENTS.md and working notes.
Reason: the user requested a new project and explicit project rules.
Cost: shared components need review and adaptation rather than implicit imports or copying.

## 2026-09-11 — Adopt the seven general PoliScreen rules

Decision: preserve scientific integrity, privacy, runtime path discovery, English project
files, workspace hygiene, compatibility and runtime verification. Adapt project names and
remove claims about configuration and tests that do not yet exist in TOPICS.
Reason: the user explicitly requested the documented rules and their maintenance workflow.
Cost: notes and evidence must be maintained alongside implementation.
Historical PoliScreen tasks and retracted findings are not copied as active TOPICS decisions.

## 2026-09-11 — Minimal validated core first

Decision: use the previously discussed CLI, input validation, provenance and conformer
milestone before ternary modeling, UI and large-scale screening.
Reason: each scientific stage needs a verifiable input/output contract.
Cost: advanced capabilities remain pending until their methods and tests exist.

## 2026-09-11 — Measured baseline and constrained resource defaults

Decision: implement rigid conformer-driven assembly from known local binding poses, with
one CPU thread by default and optional CUDA distance blocks. Retain all converged conformers
and rank with explicit geometry rules; do not use reference error for selection.
Reason: this is a concrete, testable baseline for the requested terminal application.
Cost: the method lacks protein flexibility, PPI energies and an intact ligase assembly.

Decision: report the large observed errors in both retrospective controls without tuning
parameters against their native structures. Treat both as development controls, not held-out
validation. Record top-ranked and best sampled results separately.
Reason: implementation correctness and scientific accuracy are separate claims.

Decision: standalone HTML and vector/600-dpi figures are generated from the same saved data.
Reason: interactive inspection and scientific export were explicitly requested.
Cost: embedding Plotly makes a report several megabytes but allows offline use.
