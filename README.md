# TOPICS

An interactive terminal application for conformer-driven POI–PROTAC–E3 assembly,
crystallographic controls, standalone HTML reports and exportable scientific figures.

**Version 0.1 is a structural baseline, not a validated predictor.** Its initial controls show
substantial prediction errors. See [measured results](docs/VALIDATION.md).

The [accuracy research](docs/ACCURACY_RESEARCH.md) audits the baseline and compares published methods.
Read the [implementation plan](docs/ACCURACY_PLAN.md), [offline report](docs/ACCURACY_RESEARCH.html)
or [PDF](docs/ACCURACY_RESEARCH.pdf). P0 now provides [independent assessment](docs/ASSESSMENT.md)
and a [protected 32-case inventory](benchmarks/README.md); sampling and ranking improvements
remain pending. [Full-ensemble results](docs/P0_RESULTS.md) preserve the measured baseline.

## Install and start

Python 3.11+ is required. From the project directory:

```bash
uv sync
uv run topics
```

Without uv, create a virtual environment and run `python -m pip install -e .`, then `topics`.
The terminal menu runs a control, assembles custom inputs, evaluates a run or checks devices.
Option 5 assesses every saved candidate independently; optional DockQ/PoseBusters setup is
documented in [ASSESSMENT.md](docs/ASSESSMENT.md).

```bash
uv run topics doctor
uv run topics benchmark examples/5t35.json --out runs/5t35
uv run topics benchmark examples/6boy.json --out runs/6boy
```

Controls download the selected PDB structure and chemical component dictionary. Use
`--data-dir <data directory>` for local `<PDB ID>.cif` and `<CCD ID>.cif` files.
Compressed `pdb_0000<lowercase PDB ID>.cif.gz` coordinate files are also accepted.
Review [PDB usage and citation](https://www.rcsb.org/pages/usage-policy).

## Resource use and GPU

Defaults: **one CPU thread**, 16 conformers, one seed, CPU distances and a 50-ms pause between
candidates. Coordinate storage has an estimated 256 MiB budget. Split large jobs into runs.
No background workers or servers are launched.

CUDA accelerates distance blocks, not RDKit conformer generation. For a compatible NVIDIA GPU:

```bash
uv sync --extra gpu
uv run topics doctor --device cuda
uv run topics benchmark examples/5t35.json --out runs/5t35-gpu --device cuda
```

`--device auto` reports a CPU fallback if CUDA cannot run; explicit `cuda` reports an error.
Small jobs may not run faster on GPU: the initial four-conformer comparison took 2.82 s on
CPU and 2.84 s with GPU distances. No temperature guarantee is possible across hardware.
Use `--threads`, `--conformers`, `--seeds` and `--pause` to control workload.

## Custom inputs

```bash
uv run topics inspect protein.cif
uv run topics assemble input.json --out runs/custom --conformers 32 --seeds 42 43 44
uv run topics evaluate runs/custom --reference reference.npz
```

Input JSON defines two local protein–head poses, explicit CCD atom-name mappings and the
full PROTAC dictionary. Paths resolve relative to the JSON. See [methods](docs/METHODS.md)
and a generated `runs/5t35/inputs/input.json` for the input contract. Reference coordinates
are a separate evaluation input and are never used by the assembler or ranking function.

## Outputs

Open `runs/5t35/run/report.html` in a browser. It works offline, with rotatable 3D coordinates,
candidate selection, metric hover, a table and downloadable artifacts.

- `manifest.json`: versions, seeds, settings, input/output hashes and completion state.
- `candidates.csv`, `ensemble.sdf`, `ensemble.npz`, `best.cif`: ranking and structures.
- `evaluation.json`, `evaluation.csv`: reference metrics and per-seed results, when evaluated.
- `assessment.*`, `assembly.*`: SVG, PDF and 600-dpi PNG; captions in `figure_captions.txt`.
- `view.pml`: PyMOL script; open it from the run directory.

Protein figures use alpha-carbon points. Benchmark inputs include `positive_control.json`:
its near-zero error checks mapping and rigid transformations, not predictive accuracy.
Outputs require new directories; previous predictions are never silently overwritten.

## Development

```bash
uv sync --extra dev
uv run pytest -q
```

Set `TOPICS_TEST_DATA` to local PDB/CCD files to enable biological geometry checks. CUDA tests
skip without a runtime. See [AGENTS.md](AGENTS.md), [rules](_notes/RULES.md) and
[short changelog](CHANGELOG.md). No PoliScreen implementation was copied.
