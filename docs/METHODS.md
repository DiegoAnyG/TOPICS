# Methods and input contract

## Scope

TOPICS 0.1 accepts a complete, connected PROTAC and known local binding poses for both heads.
It generates free PROTAC conformers and places rigid protein partners on the corresponding
heads. It does not yet design linkers, dock unknown heads, refine protein interfaces or model
an intact ubiquitination assembly. Binding, cooperativity and degradation are not estimated.

Linker conformations and binary poses have precedent in [PRosettaC](https://pmc.ncbi.nlm.nih.gov/articles/PMC7592117/).
TOPICS is a separate, simpler baseline and does not claim that protocol or its performance.

## Input JSON, schema 1

```json
{
  "schema_version": 1,
  "title": "My assembly",
  "ligand_ccd": "ligand.cif",
  "poi": {
    "structure": "poi.cif",
    "head_atom_names": ["C1", "C2", "N1"],
    "head_coordinates": [[0, 0, 0], [1, 0, 0], [0, 1, 0]]
  },
  "e3": {
    "structure": "e3.cif",
    "head_atom_names": ["C20", "C21", "N5"],
    "head_coordinates": [[0, 0, 0], [1, 0, 0], [0, 1, 0]]
  }
}
```

Numbers illustrate the format, not a physical pose. Supply real head coordinates in the
frame of the corresponding protein. Atom names must exist in the CCD; each head needs at
least three unique, non-collinear atoms. Head selections must be disjoint. Prefer rigid head
atoms. Proteins may have multiple chains but exactly one model; curate missing residues first.

Only amino-acid heavy atoms are included. Hydrogens, waters and non-amino-acid cofactors are
excluded: a material limitation for metal/cofactor-dependent systems. Alternate atoms use
highest occupancy, then alternate-location label. Author chain, residue number, insertion code,
residue name and atom name are retained. Output chains gain a role prefix to avoid collisions.
Gemmi reads coordinates. Ligand bond orders, charges, atom names and chirality come from CCD.
Ideal CCD coordinates establish chirality and are removed before conformer generation.

## Sampling and placement

- RDKit ETKDGv3, random coordinates, 500 embedding iterations, enforced chirality, explicit
  seeds and threads; [RDKit documentation](https://www.rdkit.org/docs/RDKit_Book.html).
- MMFF94s, maximum 1500 iterations. Missing parameters cause an error; unconverged or
  chirality-inconsistent conformers are rejected. Counts are recorded.
- Proper-rotation Kabsch fits the generated POI head onto its binary pose. The E3 binary
  pose is fitted onto the other head of the transformed PROTAC.
- Rank lexicographically by total inter-component heavy-atom clash pairs below 2 Å, maximum
  head-fit RMSD, descending PPI residue-contact count, MMFF energy, then candidate ID.
  Geometric scores are heuristics; conformer energy is not binding free energy.
- A contact is a residue pair with any heavy-atom distance at most 5 Å. POI/E3, POI/PROTAC and
  E3/PROTAC are counted separately. These are proximity contacts, not typed hydrogen bonds.

There is no clustering or protein relaxation in 0.1. All accepted conformers are stored.
GPU distances use float64 CuPy operations in blocks of 128 source atoms; CPU uses SciPy cdist.
Blocks bound temporary memory. Transfer overhead can outweigh GPU benefits for small inputs;
see [CuPy guidance](https://docs.cupy.dev/en/stable/user_guide/performance.html).

## Evaluation

Control specifications name author chains, ligand residue and chemically explicit head atoms.
Binary partners are extracted from the ternary crystal and placed in separate local frames.
Bound protein shapes and local head poses are retained: this is retrospective bound-component
testing, not independent prediction from apo or separately determined binary structures.

`reference.npz` has arrays `poi`, `e3`, `ligand`, `poi_keys`, `e3_keys`, `atom_names` in the
prepared input's heavy-atom order, loaded without pickle. Evaluation fits predicted POI heavy
atoms onto the crystal, then computes E3 alpha-carbon RMSD and whole-ligand heavy-atom RMSD.
Exact atom names are used; symmetry permutations are not searched. Native PPI contact recall
and precision use the 5 Å residue-pair definition. These are not DockQ or uncertainty estimates.

Top-ranked and best sampled reference matches are separate. The latter is an oracle
diagnostic, never a prospective selection. Each seed has a summary; conformers are correlated
samples, not independent biological replicates.

The positive control assembles the binary partners using native ligand geometry. It checks
implementation consistency only. Automated tests also check translated-partner negatives and
independent frame rotations. General accuracy claims require a larger held-out benchmark.

## Provenance and resource limits

Manifests record raw-byte SHA256 inputs/outputs, versions, seeds, atom mappings, parameters
and computation time excluding reports. No personal paths or hostnames are recorded.
Each new run recalculates predictions. Incomplete runs keep a non-complete manifest state.

One CPU thread and pauses are defaults. A finite conformer budget, maximum eight explicitly
requested CPU threads and a 256 MiB estimate for stored coordinates constrain resource use.
These are not a thermal controller or a hard limit on all library allocations. Ctrl+C stops
the foreground job; no child workers or servers are launched.
