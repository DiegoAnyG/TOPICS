# Log

## 2026-09-11

- Established the project documentation and mandatory context-reading workflow.
- Preserved the seven general rules with TOPICS-specific references and explicit scientific
  interpretation requirements.

## Lessons

Adapted from PoliScreen's history; these are methodological lessons, not TOPICS findings.

- Read the tree before claiming that a capability is absent or already implemented.
- A plausible cause is a hypothesis until a discriminating measurement supports it.
- Check later corrections before repeating historical claims.
- Trace reproducibility across the entire pipeline, including conversion and caches.
- Check response schemas before treating an external service's output as empty.
- Distinguish numerical repeatability from scientific accuracy and sampling uncertainty.
- Record descriptor definitions consistently across filtering, ranking and reports.
- Review licences and dependency constraints before reusing code, tools or datasets.

## Skill installation — 2026-09-11

Installed the standalone `ponytail` skill from `DietrichGebert/ponytail`,
commit `356918eba965ee1eac64bd3a7f0dd02108350de5`, path `skills/ponytail`, into the user skill directory.
No plugin lifecycle hooks were installed. The skill is available on the next turn.

## 2026-09-11 — TOPICS 0.1

- Implemented an independent Python package, interactive terminal and reproducible uv environment.
- Downloaded 5T35 and 6BOY plus their CCDs through PDB tools; recorded sources and explicit selections.
- Verified native-geometry reconstruction separately from de novo conformer-based placement.
- Recorded large prediction errors; no scoring adjustment was made to favor these controls.
- Verified float64 CPU/CUDA distance agreement and a four-conformer ranking/coordinate comparison.
- Kept scientific runs regenerable and out of Git; small validation summaries are versioned.

Lesson: a near-zero native-geometry control proves mapping consistency, not predictive power.

## 2026-09-12 — Delivery

- Finalized documentation formatting and prepared the verified application and local reports for the project workspace.
- Scientific code and the recorded control results remain unchanged from the verified 2026-09-11 runs.

## 2026-09-12 — Accuracy investigation

- Reviewed the engine, chemistry, geometry, structure selection and evaluation against
  published ternary predictors, independent benchmarks and current upstream documentation.
- Added 24-source research and a staged plan; explicit ligand inputs, binary/ternary mixtures,
  best-of-N selection and training overlap explain important differences among reported results.
- Audited native, rank-1 and lowest-E3-RMSD candidates with DockQ 2.1.3 in an isolated environment.
  Rank-1 DockQ values were 0.01585 and 0.06362; native self-comparisons approached one.
- Measured ligand binding-region deviations after independent protein alignment; the
  5T35 rank-1 E3-associated atoms deviate by 7.50 A despite 0.277 A maximum ring-fit RMSD.
- Repeated the entire audit and obtained byte-identical JSON; checked exported native
  identity, saved clash counts, document links, offline navigation and a 12-page PDF.
- Scientific source, production dependencies and original prediction runs remain unchanged.

Lesson: a PPI DockQ of acceptable quality can coexist with severe intercomponent clashes;
structural similarity and chemical validity require separate, jointly reported endpoints.
Lesson: preserve complete chemical binding groups; a small-ring alignment does not validate
the pose of the rest of the ligand or the orientation of a distant protein domain.
