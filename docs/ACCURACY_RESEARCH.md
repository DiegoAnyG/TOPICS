# TOPICS Accuracy Research and Development Plan

## 1. Assessment

TOPICS needs a better model of the coupled protein–PROTAC–protein geometry before a larger conformer count can yield dependable predictions. Its current implementation samples a free PROTAC, fits small atom selections to two binding sites, and places otherwise rigid proteins around it. This provides a reproducible baseline, but does not preserve the complete bound ligand groups, search protein interfaces independently, or optimize the assembled complex in its molecular environment.

The recommended sequence is: establish independent evaluation; preserve complete binding groups while sampling linker and protein orientations; reject chemically invalid poses; refine a small, diverse ensemble; then improve ranking against a separate validation set. Test DeepTernary as an optional alternative generator and Boltz-2 as an external comparator. Keep their original confidence measures separate from TOPICS heuristics until validation supports combining them.

This assessment covers heterobifunctional PROTAC structure prediction and evidence available through 12 September 2026. Molecular glues and cellular degradation prediction are adjacent tasks with different input assumptions. The quantitative baseline is TOPICS commit `b2ae904bcd605b58c2e29e533e16305d8d2e5cf9`; the new audit evaluates saved coordinates without changing generation or ranking. Improvements below are proposals, not measured performance gains.

## 2. What the current controls establish

### Structural accuracy and physical validity

Each control contains 96 accepted conformers: 32 per seed, using seeds 42, 43 and 44. Both protein conformations and local binding poses originate from the same ternary crystal, making this an optimistic retrospective reconstruction task. The native relative arrangement is withheld from assembly and ranking, but these are not independent binary-input or blind predictions. The existing [validation record](VALIDATION.md) documents the experiment and limitations.

| Measurement | 5T35: BRD4 BD2 / MZ1 / VHL | 6BOY: BRD4 BD1 / dBET6 / CRBN |
|---|---:|---:|
| Crystal resolution, Å | 2.70 | 3.33 |
| Rank-1 E3 Cα RMSD after fitting POI, Å | 51.18 | 48.88 |
| Lowest sampled E3 Cα RMSD, Å | 11.73 | 14.74 |
| Rank of that lowest-RMSD candidate | 24 | 6 |
| Rank-1 PPI native-contact recall | 0 | 0 |
| Candidates with zero current clash pairs | 0 / 96 | 0 / 96 |
| Rank-1 current clash pairs | 22 | 68 |
| Lowest-E3-RMSD candidate clash pairs | 238 | 176 |
| Native complex current clash pairs | 0 | 0 |

The current clash definition counts intercomponent heavy-atom pairs below 2 Å. It is a coarse geometric diagnostic, not a complete steric energy. Nevertheless, the native complexes pass this same check and every prediction fails it. A ranking change cannot create a clash-free pose absent from the sampled ensemble. Conversely, proximity to the native orientation does not guarantee a usable geometry: the lowest-RMSD candidates still have extensive overlaps.

An independent audit with DockQ 2.1.3 gives the following PPI scores. Only rank 1, the existing lowest-E3-RMSD candidate, and native self-comparisons were evaluated with DockQ. These values therefore do not establish the maximum DockQ among all 96 predictions.

| Audited candidate | DockQ | Interface RMSD, Å | Current total clash pairs |
|---|---:|---:|---:|
| 5T35 rank 1 | 0.0159 | 11.81 | 22 |
| 5T35 rank 24, lowest E3 RMSD | 0.3163 | 4.75 | 238 |
| 6BOY rank 1 | 0.0636 | 7.71 | 68 |
| 6BOY rank 6, lowest E3 RMSD | 0.2088 | 4.90 | 176 |
| Native self-comparison, either system | approximately 1 | approximately 0 | 0 |

DockQ uses 0.23, 0.49 and 0.80 as boundaries for acceptable, medium and high PPI quality.[^1] The 5T35 rank-24 candidate crosses the first boundary while remaining physically unsuitable. TOPICS must report structural agreement and physical validity together. DockQ's own clash count differs from TOPICS' all-component atom-pair count; their values must retain distinct names. DockQ also treats the larger E3 chain as receptor here, whereas the existing E3 RMSD fits POI first. The RMSDs are consequently not interchangeable.

### Ranking mostly considers one feature

`rank_candidates()` in [engine.py](../src/topics/engine.py) sorts lexicographically by clash count, head-fit RMSD, negative PPI contact count, free-ligand MMFF energy, and candidate identifier. Later criteria only operate when all previous values tie. Clash counts are unique for 76/96 candidates in 5T35 and 88/96 in 6BOY; their placement is fixed before the other features matter. Even in a clash tie, a continuous head-fit value usually prevents contact or energy terms from influencing order.

This behavior is reproducible, but it is not a balanced physical score. Within these two ensembles, descriptive Spearman correlations between head-fit RMSD and E3 error are approximately −0.05 and +0.03. They do not establish a general statistical result, but provide no support for treating this local fit as a global accuracy proxy. Raw PPI contacts also increase when proteins interpenetrate, making a larger contact count unsafe as an unconditional reward.

### The bound ligand groups are insufficiently constrained

The examples fit a five-atom triazole ring on the POI side and five or six atoms on the E3 side. The complete PROTAC remains flexible during embedding and free-molecule minimization. A good fit for one ring cannot preserve all other atoms of the binding group or its exit-vector geometry.

The audit independently fits each predicted protein onto its native counterpart and measures the ligand atoms associated with that local binding site. For this diagnostic only, an atom is included when it lies within 4 Å of its native protein and more than 5 Å from the other native protein. These reference-derived sets are not chemical fragment definitions and must never enter blind prediction inputs.

| Rank-1 local diagnostic | 5T35 | 6BOY |
|---|---:|---:|
| Maximum fitted small-ring RMSD, Å | 0.277 | 0.352 |
| POI-associated ligand atom count | 17 | 15 |
| POI-associated ligand RMSD, Å | 4.11 | 1.03 |
| E3-associated ligand atom count | 20 | 16 |
| E3-associated ligand RMSD, Å | 7.50 | 4.61 |

The loss of local binding geometry is demonstrated. How much each proposed constraint improves final docking remains to be tested. Nearly planar rings are not, by themselves, an invalid Kabsch fit: noncollinear planar points determine a proper rigid alignment. The main concern is the limited chemical coverage and the amplification of small orientation changes at distant protein atoms, not an established singular-matrix bug.

### Other limitations visible in the code

`chemistry.py` minimizes isolated PROTAC conformers; its MMFF energy omits protein interactions and solvent. `geometry.py` counts distance contacts without donor/acceptor direction, charge, desolvation or atom-specific radii. `structures.py` removes non-amino-acid components from protein clouds, including waters and metals. No side-chain relaxation, protein-domain ensemble, clustering or bound-state strain assessment is implemented.

`benchmark.py` currently requires native PPI contacts and exact reference atom ordering. This is suitable for the two curated controls but will reject legitimate inputs requiring sequence mapping, missing-residue reconciliation or evaluation without a direct PPI interface. These are specific engineering limitations, not explanations that make the observed errors acceptable.

## 3. Relevant methods and the strength of their evidence

### Established docking pipelines

Rovers and Schapira's 2024 benchmark compares PRosettaC, MOE and ICM. Its central practical finding is that obtaining a near-native candidate and selecting it reliably remain distinct difficulties.[^2] The accompanying public repository is a useful starting point for annotated controls and comparable input preparation.[^3] Imported cases still need explicit review of chains, ligand completeness, assembly and experimental provenance.

PRosettaC combines linker-derived restraints, protein docking, local refinement, PROTAC placement and clustering.[^4] It supports the decision to add a protein-interface search rather than relying exclusively on free PROTAC conformations. Its public implementation depends on PatchDock, Rosetta, OpenBabel and scheduler-oriented scripts; the wrapper's MIT license does not settle the terms or installation requirements of those dependencies.[^5] Use it as a benchmark comparator where practical, rather than make the entire stack a laptop dependency.

| Method | Useful contribution to TOPICS | Qualification |
|---|---|---|
| PRosettaC | Constrained protein docking, local refinement and cluster representatives | External dependencies; published success does not imply reliable rank 1.[^4], [^5] |
| Ignatov et al., 2023 | Simultaneous interface and linker search using complementary half-linker conformational clouds | Twelve crystal structures; high-quality models among the top ten, with selection uncertainty. Ubiquitination filters are a separate biological assumption.[^6] |
| BOTCP, 2023 | Bayesian optimization to use evaluations more selectively | Original method and bibliographic record identified; detailed runtime and deployment claims are not independently established here. Treat as a sampling design reference.[^7] |
| MEGA PROTAC, 2025 | Fast docking proposals, inexpensive filters first, local search, contact clustering and rank aggregation | Twenty-two cases; comparisons include BOTCP before refinement. Filtering reduced the reported mean maximum DockQ from 0.646 to 0.554: early rejection can lose good candidates.[^8] |
| COMPASS, 2026 | Linker compatibility and optional geometry of the complete ubiquitination machinery | Twenty crystal controls and a separate 112-PROTAC activity study. Published workflow includes Glide and lengthy initial PPI docking; unsuitable as an assumed cheap default.[^9] |

The reusable principle is to search linker compatibility and protein orientation together, retaining diverse candidates until enough evidence exists to reject them. It does not require copying an entire pipeline. A first implementation can use the installed RDKit and SciPy dependencies, while an external docking backend is evaluated separately.

COMPASS reports 93% recall on degradation endpoints; recall is not precision or proof that a predicted structure is correct. Its initial protein docking and clustering took 18–24 hours per system; the faster assembly stage followed this preparation on a 32-core, 96-GB node. Its structural comparisons use selected cluster representatives and an E3-based alignment. These distinctions prevent comparing its headline numbers directly with TOPICS rank-1 E3 RMSD.[^9]

### Learned ternary generators and cofolding models

DeepTernary is the most directly relevant optional learned generator. The 2025 Nature Communications paper evaluates 22 PROTAC complexes using separate ligand-bound protein structures, samples 40 initial conformations, and ranks outputs by predicted aligned error. Its prominently reported mean DockQ of 0.65 should not be presented as rank-1 performance: the paper separately reports mean top-1 DockQ around 0.4. Training excludes known PROTACs and similar protein pairs under its stated split.[^10] Reproduction on independently curated TOPICS inputs is needed before transferring any performance claim.

The official repository supplies inference scripts, including a CPU route, and uses an Apache-2.0 code license. It also explicitly notes rigid chemical-handle constraints in conformational sampling.[^11] This makes a small isolated evaluation worthwhile. Exact checkpoint provenance, package compatibility, memory peaks and run time on an 8-GB laptop GPU remain unmeasured; published second-scale inference is not an end-to-end local guarantee.

PROflow explores learned iterative refinement with full PROTAC flexibility and synthetic training examples derived from protein complexes. It is a relevant 2024 workshop/preprint direction, but an official deployable checkpoint and reproducible local setup were not established in the available evidence. It should remain a research option rather than a first implementation dependency.[^12]

The apparent disagreement between AlphaFold comparisons is substantially explained by inputs and metrics:

- A 2025 Scientific Reports comparison used an AlphaFold server protocol described with protein sequences; it does not establish performance for local, explicit-PROTAC conditioning. Accessory chains also affect aggregate interface scores.[^13]
- Dunlop et al. used explicit ligands and compared AF3 with Boltz-1. Their 62-entry collection contains **48 ternary and 14 binary complexes**, not 62 independent ternary predictions. Ligand representation and release date affected results; performance was poorer on later structures. The associated PROTACFold pipeline and data are useful references.[^14]
- Riepenhausen et al. benchmarked AF3 and Boltz-2 in 2026 on 25 PROTAC and 15 molecular-glue complexes. Boltz-2 performed better in their evaluation, but the supporting tables explicitly use the best RMSD and best DockQ among five samples. Confidence-selected results must be compared separately.[^15]
- Chen et al.'s 2026 comparison of Chai-1, AF2, AF3 and Protenix identifies remaining errors in relative protein orientation and PROTAC placement despite favorable whole-complex results. Only the accessible abstract and data-availability description support that summary; no detailed leaderboard is inferred.[^16]

The recommendation is to benchmark DeepTernary for the existing binary-structure workflow and Boltz-2 for a sequence-plus-ligand workflow. Use AF3 and Chai-1/Protenix as additional external comparators if resources permit. Do not choose a universal winner by combining incompatible reported averages. A model supplied with two bound binary poses and a model supplied with sequences solve different input problems.

## 4. Proposed changes to the scientific pipeline

### Complete chemical head definitions and constrained linkers

Extend the input schema with explicit POI-head, E3-head and linker atom membership, attachment bonds, stereochemistry and local bound poses. Require complete, nonoverlapping atom maps covering the molecule; display ambiguous mappings for correction. Preserve schema-1 inputs as a clearly labeled legacy mode. Provide general SDF/SMILES input without requiring a PROTAC already present in the CCD.

For each head, preserve its intramolecular geometry and its relationship to the associated protein. The relative transform between the two protein–head units must remain an unknown to be sampled. In particular, applying coordinate constraints to both heads in their original ternary reference frame would leak the answer and is prohibited.

Compare two generators using identical chemical inputs and budgets. The first retains complete head geometry while sampling the linker to propose relative protein positions. The second proposes protein–head orientations and solves linker closure subject to bond, angle, torsion and chirality constraints. Share the feasibility checks and output format. Start with small batches and explicit failure reasons; increase sampling only when it adds distinct feasible solutions.

RDKit already supports constrained embedding and positional restraints during minimization. Its documentation explains that distance-based coordinate maps do not always preserve exact coordinates, and subsequent unconstrained minimization can move the constrained core.[^17] Check the complete heads after both embedding and refinement; do not assume that supplying `coordMap` completes the task. Retain TOPICS' chirality checks, including at newly interpreted attachment sites.

### Validity filtering and limited refinement

Replace the universal distance threshold as the sole physical criterion with typed nonbonded overlaps, appropriate exclusions, bond/angle plausibility, ring planarity, stereochemistry and complete-head pose preservation. Retain the old count for baseline comparison. Separate protein–protein, POI–PROTAC and E3–PROTAC checks so failures remain interpretable.

PoseBusters provides established chemical and geometric validation for predicted ligand poses.[^18] Evaluate it as an optional audit tool and explicitly define which checks apply to large bifunctional ligands and retained cofactors. PPI overlap checks remain necessary in addition to ligand validation. A failed checker or unparameterized atom must become a recorded failure, not a zero penalty.

Refine only diverse candidates that pass coarse feasibility checks. Begin with constrained ligand minimization, then add local side-chain relaxation while restraining the protein backbone and preserving each head's binding geometry. Use compatible protein and small-molecule force fields with explicit protonation, charges and metal treatment. Avoid arbitrary protein truncation without specifying boundary conditions; local flexibility can be restricted while retaining the full interaction environment.

OpenMM supplies bounded local energy minimization; an explicit iteration limit is needed because its default permits unlimited iterations.[^19] A local minimizer can relieve nearby steric problems but is not expected to repair a protein oriented tens of Å from the target pose. Long molecular dynamics or MM/GBSA across hundreds of invalid candidates should therefore be deferred. If MD is later used, evaluate predefined hypotheses on several seeds and report instability as well as stability.

### Ranking, clustering and confidence

Use physical feasibility as an explicit gate, then rank feasible candidates with a small documented set of features: ligand strain relative to a same-protomer reference ensemble, typed interactions, buried surface and unsatisfied polar groups, local binding-pose distortion and interface complementarity. Free-ligand MMFF energies from different molecules or protonation states are not directly comparable binding scores.

An initial consensus ranking or Pareto selection can avoid mixing unrelated units, but remains a heuristic requiring validation. Do not fit arbitrary weights on 5T35 and 6BOY. Introduce a learned ranker only after obtaining sufficiently diverse systems, with system-grouped validation. Preserve individual feature values and rejection reasons in outputs.

Cluster feasible complexes in a consistent protein frame using interface geometry or interaction fingerprints. Return representatives from distinct clusters rather than ten near-duplicates. ProLIF can encode typed molecular interactions and provides a reference for interpretable contact analysis.[^20] Contact fingerprints used for prediction may come from known binary poses or predicted candidates; native ternary contact recovery belongs exclusively to evaluation.

Represent uncertainty through ensemble diversity, independent-seed agreement, validity flags and measured calibration on held-out systems. Cluster population from biased conformer sampling is not a Boltzmann population. PAE, ipTM and docking scores are not automatically probabilities of structural correctness or degradation.

## 5. Benchmark design that can establish improvement

### Inputs, assemblies and exclusions

Freeze an annotated candidate inventory from the Rovers benchmark, DeepTernary input pairs, PROTACFold and the 2026 studies. These collections overlap; their sizes must not be added. Aim to curate 20–40 eligible, structurally diverse ternary cases initially, subject to actual input quality. Record the number retained after review rather than forcing a target count.

Keep 5T35 and 6BOY as development and regression controls. Separate matched bound-component reconstruction from prediction using independent binary structures, and separate both from apo or predicted-protein inputs. Identify biological assembly, author and label chain IDs, residue mapping, ligand atom correspondence, alternate conformers, occupancy, unresolved linker atoms and experimental quality. Explicitly resolve binary complexes and molecular glues into other tasks.

Preserve structural cofactors when biologically justified. Evaluate the POI–E3 interface independently of accessory chains; optionally report full-assembly compatibility as an additional endpoint. Cases without direct native PPI contacts require a documented metric fallback rather than automatic exclusion or a fabricated zero. The reported 8QJR arrangement is a useful curation candidate for this boundary; its PDB entry and associated experiment should be checked directly before eligibility is decided.[^9], [^21]

### Separation from development and training data

Group related structures by protein pair, target family and degrader scaffold/linker series before creating development, validation and final test partitions. Repeated crystal copies or near-identical analogues must remain in the same group. Report ligase-specific performance and include non-bromodomain targets. If data permit, add leave-family-out and temporal challenges rather than relying on a random PDB split.

For pretrained models, record checkpoint hashes, training cutoffs, available cluster exclusions and similarity to training structures. A later release date alone does not establish novelty: homologous proteins, related complexes and chemical dictionaries may carry overlapping information. Unknown training membership must be labeled unknown. Native-coordinate CCD features, idealized dictionary features and SMILES-generated coordinates also require separate provenance.

Decide the final test cases and success definitions before changing scientific scoring. After an algorithm or threshold is selected on validation data, freeze it and evaluate the final test once. A failed test is reported and becomes historical evidence; subsequent development needs a fresh protected assessment.

### Metrics and statistical unit

| Question | Required endpoint |
|---|---|
| Can the generator reach a good structure? | Oracle success among all candidates within a fixed budget; best DockQ and ligand RMSD, clearly reference-selected |
| Does ranking find it? | Success@1, @5 and @10 for ranked cluster representatives; conditional success when the ensemble contains a qualifying candidate |
| Is the proposed complex physically usable? | Valid-pose fraction, component-specific clashes, head preservation, chemistry and parameterization failures |
| Where is the error? | PPI DockQ, iRMSD, Fnat; role-defined partner RMSD; binding-site-aligned and symmetry-aware ligand RMSD |
| Does it generalize? | System-level distributions by ligase, target family, chemical series and input difficulty |
| Is it efficient? | End-to-end elapsed time, CPU time, RAM/VRAM peak, proposals and accepted clusters per budget |

Predefine a primary structural endpoint such as rank-1 PPI DockQ ≥0.23 and report the stricter ≥0.49 endpoint separately. Add a joint endpoint requiring chemical validity and a preregistered ligand pose tolerance; 2 Å is a useful initial convention for testing, not a proven universal threshold for full PROTACs. Report head and linker errors separately and use symmetry mappings that preserve chemical identity and stereochemistry.

Compute confidence intervals by resampling independent systems or chemical/protein-family groups, not the 96 correlated conformers of a single complex. Count timeouts, preparation failures and zero-feasible-pose outcomes in the eligible-system denominator. Present paired changes from the baseline, confidence intervals and failure cases. With only two controls, an overall accuracy percentage or a trained confidence model would be misleading.

## 6. Implementation sequence and decision gates

The following effort ranges are engineering estimates for one contributor, excluding scientific uncertainty, large downloads and external computations. They are not delivery promises or expected effect sizes. The operational checklist is in [ACCURACY_PLAN.md](ACCURACY_PLAN.md).

| Phase | Work | Estimated effort | Gate before advancing |
|---|---|---|---|
| P0 | Evaluation adapter, atom/residue mapping, frozen benchmark manifest, split and budget definitions | 3–5 working days | Native identity and rigid-transform invariance pass; invalid/empty cases report correctly; no reference enters prediction |
| P1 | Chemical head/linker partition, complete-head constraints, two constrained sampling variants | 5–8 days | Head geometry and stereochemistry survive each stage; feasible yield and oracle performance are measured on development/validation |
| P2 | Coarse PPI proposals, linker closure, selective local refinement and diversity clustering | 5–10 days | Improvement exceeds sampling noise under matched budget; chemistry and known binding modes are retained |
| P3 | Physical feature ranking and isolated DeepTernary comparison; optional Boltz-2 comparator | 4–7 days | Rank-1/5 gains reproduce on validation groups; resource ceilings and preparation failures are recorded |
| P4 | Frozen test evaluation, confidence reporting, publication figures and reproducibility package | 3–5 days | All eligible cases accounted for; paired uncertainty reported; source/version traceability complete |

Total sequential engineering scope is approximately 4–7 working weeks, with P3 model evaluation possible after P0 if it does not consume the resources needed for P1/P2. If P1 produces no feasible conformers, diagnose chemical constraints and closure before increasing scale. If feasible near-native structures exist but rank poorly, prioritize P3 scoring. If all generators fail one family, inspect input quality and protein flexibility before adding more score features.

Useful ablations are: unchanged baseline; complete-head constraints only; constraints plus alternative sampling; addition of local refinement; addition of clustering/ranking; learned generator with its native ranking; learned generator plus the same TOPICS validation/refinement. Keep inputs, evaluation and compute accounting fixed. Re-run the original baseline only when the measurement method changes, and retain both old and new metric versions.

## 7. Resource and integration strategy

Retain CPU operation as the default portable route. Use one or two compute threads, bounded batches and a wall-time limit for routine exploration. A proposed pilot budget is two minutes per control in the conservative mode and ten minutes in an explicitly selected extended mode, to be revised from measured completion rates. These are experimental budgets, not a claim that the full planned pipeline currently meets them.

The existing four-conformer parity check took 2.82 seconds on CPU and 2.84 seconds with CUDA distances. This small experiment shows numerical agreement, not a GPU speed advantage. Profile conformer generation, scoring, refinement and reporting separately before optimizing. CPU neighbor searches and cached static protein features may help short-range scoring more than transferring small distance blocks to a GPU. Do not change cutoff definitions when changing the search algorithm.

For learned models or OpenMM, use one GPU job at a time, batch size one initially, memory preflight, cancellation and restartable outputs. Measure peak VRAM on the actual 8-GB GPU before setting an admission limit, leaving headroom for the display and other processes. OpenMM exposes thread limits and CUDA synchronization/precision controls.[^22] Numerical precision affects computational reproducibility; it does not by itself correct an inaccurate physical model. Hardware temperature cannot be guaranteed by a fixed pause, so describe resource limits honestly and use telemetry only where available.

Install optional engines in isolated environments and call their documented interfaces. Boltz's current official repository states that code and weights are MIT-licensed and recommends a fresh environment.[^23] The AF3 repository currently labels its code Apache-2.0 and provides separate model-weight and output terms, superseding older papers' code-license descriptions.[^24] Check the exact versions and terms before redistribution. Do not bundle checkpoints, start Docker automatically, or assume a paper's cluster timing applies to a laptop.

The minimal near-term dependency change is an optional evaluation environment for DockQ/PoseBusters. Complete-head sampling can start with existing dependencies. Add OpenMM, a learned generator or an external docking engine only when the corresponding experiment reaches its gate.

## 8. Reporting and limits of interpretation

The interactive report should expose cluster representatives, the three interfaces, complete-head deviations, chemical validity, feature-level ranking and failed-stage counts. Reference overlays and native-contact recovery should appear only for evaluation runs. Ordinary predictions should show model confidence and ensemble disagreement without labeling unknown native accuracy as measured.

The publication package should include an input/split inventory, case-level CSV/JSON, failure table, success@k curves with system-level intervals, budget-versus-success curves, paired method comparisons, and representative structural overlays selected by a stated rule. Export vector SVG/PDF and high-resolution raster alternatives with units, accessible colors and independently reproducible source data. Preserve the legacy validation figures as the baseline rather than replacing them with an improved run under the same name.

Structural agreement, ternary binding, cooperativity, ubiquitination and cellular degradation remain separate endpoints. A geometrically plausible pose does not supply DC50, Dmax, binding free energy or a kinetic rate. Buried surface and proximity to an E2/ubiquitin assembly may be useful descriptors, but should not become universal activity thresholds. Any subsequent activity model needs matched assays, concentrations, cell lines, time points and negative examples, and must be evaluated separately from the structural benchmark.

## 9. Reproducible audit

The machine-readable evidence is [accuracy-audit.json](accuracy-audit.json). The accompanying [audit_accuracy.py](audit_accuracy.py) reads existing bound-component runs, checks identities and saved clash counts, evaluates selected poses, and records input/output hashes. It does not create predictions or change their ordering. Native self-comparisons validate the export and metric invocation.

Run the script with the TOPICS Python environment and a separate Python environment containing `DockQ==2.1.3`; supply paths appropriate to the installation:

```bash
python docs/audit_accuracy.py \
  --case-root runs/5t35 --case-root runs/6boy \
  --output <new-audit-file.json> \
  --dockq-python <evaluation-environment>/bin/python
```

On Windows, use the evaluation environment's `Scripts/python.exe`. Omitting `--dockq-python` runs the local geometry diagnostics only. Raw run folders remain outside the versioned evidence; the original controls can be regenerated using the instructions and inputs in `VALIDATION.md`. The six DockQ calculations here are an audit, not the complete evaluation adapter proposed in P0.

## Sources

[^1]: C. Mirabello and B. Wallner. “DockQ v2: improved automatic quality measure for protein multimers, nucleic acids, and small molecules.” *Bioinformatics* 40, btae586 (2024). [Paper](https://doi.org/10.1093/bioinformatics/btae586); [official implementation and score categories](https://github.com/wallnerlab/DockQ). Audit version: 2.1.3.
[^2]: E. Rovers and M. Schapira. “Benchmarking Methods for PROTAC Ternary Complex Structure Prediction.” *JCIM* 64, 6162–6173 (2024). [DOI](https://doi.org/10.1021/acs.jcim.4c00426). Comparative abstract and supporting-material description; no unverified runtime table reproduced.
[^3]: E. Rovers. “PROTAC ternary complex benchmark.” [Research repository](https://github.com/ERovers/PROTAC_ternary_complex_benchmark). Dataset provenance and preprocessing reference; accessed September 2026.
[^4]: D. Zaidman, J. Prilusky and N. London. “PRosettaC: Rosetta Based Modeling of PROTAC Mediated Ternary Complexes.” *JCIM* 60, 4894–4903 (2020). [DOI](https://doi.org/10.1021/acs.jcim.0c00589).
[^5]: LondonLab. “PRosettaC.” [Official repository, README and license](https://github.com/LondonLab/PRosettaC). Dependency and interface documentation; accessed September 2026.
[^6]: M. Ignatov et al. “High Accuracy Prediction of PROTAC Complex Structures.” *JACS* 145, 7123–7135 (2023). [DOI](https://doi.org/10.1021/jacs.2c09387); [author manuscript](https://pmc.ncbi.nlm.nih.gov/articles/PMC10240388/).
[^7]: A. Rao et al. “Bayesian optimization for ternary complex prediction (BOTCP).” *Artificial Intelligence in the Life Sciences* 3, 100072 (2023). [DOI](https://doi.org/10.1016/j.ailsci.2023.100072); [author publication record](https://hfooladi.github.io/publications/2023-12-01-botcp). Detailed deployment evidence remains incomplete.
[^8]: “MEGA PROTAC, MEGA DOCK-based PROTAC mediated ternary complex formation pipeline with sequential filtering and rank aggregation.” *Scientific Reports* (2025). [Paper](https://www.nature.com/articles/s41598-024-83558-2); [research code](https://github.com/yauz3/MEGA-PROTAC).
[^9]: S. Sueron et al. “COMPASS: A Computational Pipeline to Identify Linkers Predicting Ubiquitinable PROTAC-Induced Ternary Complexes.” *ChemMedChem* (15 July 2026). [Paper](https://doi.org/10.1002/cmdc.70385); [research code](https://github.com/ssueron/compass). Structural and activity endpoints are distinct.
[^10]: F. Xue et al. “SE(3)-equivariant ternary complex prediction towards target protein degradation.” *Nature Communications* 16, 5514 (1 July 2025). [Full paper](https://pmc.ncbi.nlm.nih.gov/articles/PMC12216337/). Sections on PROTAC evaluation and ranking distinguish the headline and top-1 results.
[^11]: DeepTernary authors. [Official repository](https://github.com/youqingxiaozhua/DeepTernary). Inference scripts, input-pair data description, chemical-handle constraints and code license; accessed September 2026.
[^12]: B. Qiang, W. Shi, Y. Song and M. Wu. “PROflow: An iterative refinement model for PROTAC-induced structure prediction.” GEM workshop, ICLR 2024 / preprint. [arXiv](https://arxiv.org/abs/2405.06654). Deployment not verified.
[^13]: “PRosettaC outperforms AlphaFold3 for modeling PROTAC ternary complexes.” *Scientific Reports* 15, 37620 (28 October 2025). [Paper and input protocol](https://www.nature.com/articles/s41598-025-21502-8).
[^14]: N. Dunlop, F. Erazo, F. Jalalypour and R. Mercado. “Predicting PROTAC-mediated ternary complexes with AlphaFold3 and Boltz-1.” *Digital Discovery* (27 October 2025). [DOI](https://doi.org/10.1039/D5DD00300H); [institutional full text](https://research.chalmers.se/publication/549755/file/549755_Fulltext.pdf); [PROTACFold](https://github.com/NilsDunlop/PROTACFold). Section 2.2 identifies 48 ternary and 14 binary entries.
[^15]: L. Riepenhausen, A.-C. Sarnow, D. Robaa and W. Sippl. “AI-Based Prediction of PROTAC- and Molecular Glue-Mediated Ternary Complexes: A Comparative Evaluation of AlphaFold 3 and Boltz-2.” *Archiv der Pharmazie* 359, e70225 (14 March 2026). [Paper](https://doi.org/10.1002/ardp.70225). Supporting-table description specifies best-of-five metrics.
[^16]: H. Chen et al. “Benchmarking Deep Learning for PROTAC Ternary Complex Prediction.” *Proteins* 94, 1329–1340 (2026; online 3 February). [Publisher abstract and data availability](https://doi.org/10.1002/prot.70117). Detailed full-text methods not used to make quantitative claims.
[^17]: G. Landrum. “More on constrained embedding.” RDKit developer tutorial (10 February 2023). [Technical explanation and examples](https://greglandrum.github.io/rdkit-blog/posts/2023-02-10-more-on-constrained-embedding.html).
[^18]: M. Buttenschoen, G. M. Morris and C. M. Deane. “PoseBusters: AI-based docking methods fail to generate physically valid poses or generalise to novel sequences.” *Chemical Science* 15, 3130–3139 (2024). [DOI](https://doi.org/10.1039/D3SC04185A); [author institutional manuscript](https://ora.ox.ac.uk/objects/uuid%3A898f0351-637a-454d-a51a-247277b8f4b4/files/sqj72p900h).
[^19]: OpenMM contributors. “LocalEnergyMinimizer.” [Official API documentation](https://docs.openmm.org/latest/api-python/generated/openmm.openmm.LocalEnergyMinimizer.html). Local optimization and iteration controls; accessed September 2026. Pin a stable release during implementation.
[^20]: C. Bouysset and S. Fiorucci. “ProLIF: a library to encode molecular interactions as fingerprints.” *Journal of Cheminformatics* 13, 72 (2021). [Paper](https://doi.org/10.1186/s13321-021-00548-6).
[^21]: RCSB PDB. “8QJR: BRG1 bromodomain in complex with VBC via compound 17.” [Entry and primary citation](https://www.rcsb.org/structure/8QJR). M. Berlin et al., *J. Med. Chem.* 67, 1262–1313 (2024), [DOI](https://doi.org/10.1021/acs.jmedchem.3c01781). Candidate boundary-case provenance, not an audited TOPICS result.
[^22]: OpenMM contributors. “Platform-Specific Properties.” [Official user guide](https://docs.openmm.org/latest/userguide/library/04_platform_specifics.html). CPU threads, CUDA precision, synchronization and determinism; accessed September 2026.
[^23]: Boltz contributors. [Official repository](https://github.com/jwohlwend/boltz). Current model interfaces, environment guidance and code/weight license statement; accessed September 2026.
[^24]: Google DeepMind. [AlphaFold 3 repository](https://github.com/google-deepmind/alphafold3). Current inference code and separate weight/output terms; accessed September 2026. Historical paper descriptions may reflect earlier distribution terms.
