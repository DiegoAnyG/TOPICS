# Initial structural controls

**The current baseline does not recover either crystal accurately.** These runs establish a reproducible starting point, not general predictive accuracy.

Both systems used 32 conformers for each of seeds 42, 43 and 44: 96 accepted candidates per structure. No parameters were tuned to improve these reported reference errors.

| PDB | Resolution (Å) | Top-ranked E3 Cα RMSD (Å) | Best sampled (Å) | Best sampled rank | Top native contact recall | Compute (s) |
|---|---:|---:|---:|---:|---:|---:|
| [5T35](https://www.rcsb.org/structure/5T35) | 2.70 | 51.18 | 11.73 | 24 | 0% | 70.70 |
| [6BOY](https://www.rcsb.org/structure/6BOY) | 3.33 | 48.88 | 14.74 | 6 | 0% | 42.81 |

E3 Cα RMSD is evaluated after aligning POI heavy atoms. The best sampled pose is selected using the reference only for diagnosis; it is not the method's chosen pose. Both top-ranked poses recover none of the native PPI residue contacts under the stated 5 Å definition.

## Interpretation

- Native-geometry positive controls give errors below 0.001 Å, verifying the atom mapping and rigid placement implementation. They are not predictions.
- Errors in the best sampled poses show inadequate conformational/orientational coverage at this budget. Errors in the selected poses also show limited ranking discrimination.
- Binary shapes and local ligand-head poses were extracted from each ternary crystal; this is retrospective bound-component assembly, with optimistic input information.
- 5T35 uses chains A/D and ligand 759 at D:301; Elongin B/C and the second copy are excluded. 6BOY uses C/B and RN6 at B:502; DDB1 and zinc are excluded.
- Two complexes, both containing BRD4 domains, cannot establish accuracy across proteins, ligases or chemistries. No confidence interval for general performance is claimed.

## Reproduce

```bash
uv sync --extra dev
uv run topics benchmark examples/5t35.json --out runs/5t35 --conformers 32 --seeds 42 43 44
uv run topics benchmark examples/6boy.json --out runs/6boy --conformers 32 --seeds 42 43 44
```

Runs used Python 3.11.16 on Linux/WSL with one CPU thread. Timings exclude report generation and are observations on this machine, not speed guarantees. The source and input hashes, versions and per-seed results are in [validation-results.json](validation-results.json).

Candidate metrics: [5T35 CSV](5t35-evaluation.csv), [6BOY CSV](6boy-evaluation.csv). Figures: [5T35 SVG](5t35-assessment.svg), [6BOY SVG](6boy-assessment.svg). Full offline reports and structures are generated locally under `runs/`.

CUDA validation on an RTX 4060 Laptop matched CPU distances within 1e-12 and four-conformer coordinates within 1e-10 Å, with identical candidate order and CSV metrics. Compute times were 2.82 s (CPU) and 2.84 s (CUDA): no acceleration benefit was established for this small workload.

## Next scientific work

Keep these structures as development controls. Establish a separate held-out set before adding more sampling, PPI-aware refinement, cofactor context or a better-supported selection method. Publish the original baseline alongside improvements rather than replacing its results.

See [full method definitions](METHODS.md) and the PDB entries above for the experimental structures and their primary papers.

## Software verification

- 11 automated checks passed with the local PDB data and CUDA runtime enabled.
- The wheel was built with `uv build --wheel` in an isolated build environment.
- A local Chromium check rendered the offline report, changed the candidate through the dropdown, and observed changed 3D coordinates without JavaScript errors.
- Both 96-candidate controls were repeated after finalizing the scientific implementation; their reported metrics were reproduced.
