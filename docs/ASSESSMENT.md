# Independent ensemble assessment

P0 implemented on 2026-09-12. Generation and ranking are unchanged. This adapter measures
every saved candidate and writes a new directory; it never rewrites the prediction or the
legacy `evaluate` outputs. See [the frozen inventory](../benchmarks/README.md).

## Install the optional tools

DockQ 2.1.3 requires NumPy below 2. Keep its environment separate from TOPICS:

```bash
uv venv .venv-evaluation --python 3.11
uv pip sync --python .venv-evaluation/bin/python requirements-evaluation.lock
uv run topics assess runs/5t35/run --input runs/5t35/inputs/input.json \
  --reference runs/5t35/inputs/reference.npz --out runs/5t35-assessment \
  --evaluation-python .venv-evaluation/bin/python
```

On Windows, substitute `.venv-evaluation/Scripts/python.exe`. The CLI also reads
`TOPICS_EVALUATION_PYTHON`; this enables external evaluation from interactive menu option 5.
Without an evaluation Python, local geometry and stereochemistry still run, while external
metrics remain unassessed. Missing tools and timeouts produce `partial` status, retained
candidate failures and a nonzero CLI exit. The external budget defaults to 1,800 seconds;
use `--timeout` to change it. Invalid inputs fail before successful results are reported.

The worker uses one thread, processes poses sequentially, and clears upstream DockQ caches
between poses. `mol_fast` avoids the expensive PoseBusters conformer-energy ensemble. The
evaluation lock records all resolved runtime packages; run provenance records actual versions.

## Identities and chemical head definitions

`--input` must have exactly the recorded content hashes. Saved output hashes, candidate
identities, coordinate counts, ordered SDF graphs and SDF/NPZ coordinate agreement are checked.
NPZ references use the existing `poi`, `e3`, `ligand`, `poi_keys`, `e3_keys`, `atom_names`
contract, allowing reordered atoms through complete identity lookup. No atom is dropped.

An optional `--mapping mapping.json` accepts these fields:

- `reference_selection`: required for a PDB/mmCIF reference; `poi` and `e3` contain author
  chain lists, with `ligand_chain`, `ligand`, and `ligand_residue` selecting one molecule.
- `protein_atom_map`: optional `poi`/`e3` arrays of `{"model": [...], "reference": [...]}`
  identity pairs. Each identity contains chain, residue number, insertion code, residue name
  and atom name, all as strings. Every atom must map once; incomplete structures need curated
  matched inputs. Automatic sequence alignment and partial coverage are not implemented.
- `ligand_atom_map`: complete model-to-reference CCD atom-name dictionary, when names differ.
  Chemical correspondence must be curated; this adapter cannot infer chemistry from arbitrary
  reference atom aliases.
- `head_groups`: complete, disjoint, connected heavy-atom name lists under `poi`, `e3`, and
  `linker`, accompanied by a `head_definition_source` citation. These are chemical annotations,
  not atom sets chosen by looking at reference RMSD or contacts.

Whole-head metrics remain null without those annotations. Schema-1 input selections are
reported separately as **anchor** metrics, even if they have historically been called heads.
DockQ currently accepts one selected protein chain per role. After complete identity matching,
temporary chains A/B and sequential residue numbers preserve the mapping across author offsets
and insertion codes. Multi-chain receptor assemblies require a future explicit interface adapter.

## Metric definitions

| Metric | Definition and interpretation |
|---|---|
| Legacy E3 Cα RMSD | E3 alpha carbons after fitting POI heavy atoms; original definition retained |
| Legacy ligand RMSD | Exact-name heavy atoms after the same POI fit, without a ligand fit |
| Symmetry ligand RMSD | Minimum over complete RDKit graph automorphisms with chirality enforced; at most 4,096 mappings, otherwise an explicit error |
| Local head RMSD | Each complete ligand head after independently fitting its own protein to the reference |
| PPI DockQ | Official DockQ v2 on explicitly mapped POI/E3 chains; also iRMSD, LRMSD and Fnat |
| Native contact recall (legacy) | Shared native-POI/predicted-E3 residue contacts within 5 Å divided by native contacts; use independent DockQ Fnat when POI conformation differs |
| Component clashes | POI/E3, POI/ligand and E3/ligand counts separately; retain the original 2 Å counts |
| Typed clashes | Heavy-atom distance below 0.75 times the sum of RDKit van der Waals radii; a steric heuristic |
| Ligand chemical validity | PoseBusters 0.6.5 `mol_fast` checks plus independently assigned 3D tetrahedral and double-bond stereochemistry |

**This is not full PoseBusters validity.** `mol_fast` checks loading, sanitization, InChI,
connectivity, radicals, bond lengths/angles, internal clashes, ring geometry and double-bond
flatness. Energy ratio, cofactors, waters, hydrogens and volume-overlap checks are absent.
Cross-component covalent complexes are outside this protocol.

PPI success requires DockQ ≥ 0.23. Ligand success requires symmetry RMSD ≤ 2 Å. Joint success
requires both, chemical validity, zero typed component clashes and both complete-head local
RMSDs ≤ 2 Å. It remains null when required assessments or chemical annotations are missing.
For a native with no PPI contacts, DockQ/recall are inapplicable; an exploratory E3 Cα RMSD
≤ 5 Å fallback is reported in a separate stratum. It is not a calibrated substitute for DockQ.

The summary retains every candidate in the denominator and reports coverage, success@1/5/10
and oracle success separately. Failures never count as passes. These are candidate ranks,
not distinct structural clusters. The benchmark denominator must also retain failed entire
runs; `assess` accepts only completed ensembles, while the frozen protocol defines those
case-level failures for the future benchmark runner.

## Evidence and output

`assessment.json` contains all results, native controls, versions, source hashes and criteria;
`assessment.csv` contains the full candidate table. `report.html` is self-contained and has
interactive metric plots plus source-data links. SVG/PDF and 600-dpi PNG share those data.
The original prediction report still provides the molecular 3D viewer.

The measured 192-candidate baseline is in [P0_RESULTS.md](P0_RESULTS.md). Tests cover native
identity, independent rigid transformations, reference atom permutations, explicit mmCIF
selection, whole-head metrics, missing PPI interfaces, symmetry, stereo inversion, corruption,
optional subprocess failure and protected dataset groups. CUDA is not required for this adapter.

Definitions follow the [DockQ implementation](https://github.com/wallnerlab/DockQ),
[PoseBusters API](https://posebusters.readthedocs.io/en/latest/api.html) and the pinned
[PoseBusters source](https://github.com/maabuu/posebusters); scientific context and papers are
cited in [ACCURACY_RESEARCH.md](ACCURACY_RESEARCH.md).
