# Protected structural inventory

`inventory-v1.json` freezes 43 reviewed PDB entries, 32 eligible noncovalent X-ray ternaries
and 11 explicit exclusions. This permanent benchmark data is independent from runtime tool
discovery and prediction code. `topics benchmark-check <manifest>` validates its counts and
protected family/scaffold groups.

| Split | Eligible systems | Target groups |
|---|---:|---|
| Development | 16 | BET, SWI/SNF bromodomains |
| Validation | 10 | Protein kinases, WDR5 |
| Test | 6 | BCL2 family, FKBP, RAS |

Both existing controls stay in development. Each PDB entry contributes one representative
copy, selected by ligand completeness then author identifiers, without prediction scores.
All observed copies remain recorded. Target families and complete-ligand Bemis–Murcko
scaffolds cannot cross splits; a shared VHL recruiter alone does not merge otherwise different
PROTAC scaffolds. These are small, deliberately conservative groups, not proof of broad
generalization or an exhaustive homology/chemical-series leakage audit.

The inventory is drawn from the pinned
[Rovers benchmark](https://github.com/ERovers/PROTAC_ternary_complex_benchmark/tree/fe6312f22dd6c35b3fb85bc9e8a7547a72bd029a).
Coordinates and chemical dictionaries were retrieved directly from RCSB PDB, with individual
source links, primary citation DOIs and snapshot hashes. The `pdb-database` skill queried
experimental method, resolution, release date, title and primary citation for those 43 IDs.
mmCIF inspection established author-chain/ligand selections, UniProt role annotations and
identity-operation biological assemblies. Complete heavy-atom sets were compared to CCDs.

The pinned DeepTernary list records published test overlap, **not training exclusion**.
It contains 22 complex copies from 14 unique PDB entries. Every learned-engine training
overlap remains explicitly unknown until the actual checkpoint and training data are audited.
References were inspected for curation, but TOPICS predictions have been evaluated only for
the two pre-existing development controls. No protected test predictions have been run.

Exclusions distinguish preparation limitations from experimental scope:

- 6BN8, 6BN9 and 6BNB do not contain a complete resolved heterobifunctional component;
  two also exceed the 4 Å quality limit. They are not labeled molecular glues or binary controls.
- 8BDS/8BDX fail the current CCD valence reader, 8DSO has a CCD stereo/ideal-coordinate
  disagreement, and 8FY0/8FY1 lack required CCD ideal coordinates. These are software
  preparation exclusions; revisit them in P1 without changing molecular identity or split.
- 8QU8, 8R5H and 8RWZ are electron-microscopy structures, requiring a separate stratum.

No explicit molecular-glue or binary structure is included among eligible cases. 8QJR is
retained with no native PPI contacts and separately applicable metrics. Resolution alone
does not establish coordinate certainty; local density, occupancy and alternate conformers
still matter for interpreting later results.

## Frozen protocol and readiness

The manifest records input difficulty, seeds 42/43/44, 32 candidates per seed, one CPU thread,
1,800 seconds and a 4 GiB CPU benchmark budget, separate GPU comparisons, structural/chemical
success criteria, failure denominators and group-level analysis. These are benchmark rules;
the general `assemble` command does not yet enforce a process memory or wall-time limit.
Automated budget enforcement and a case-level runner are required before a full comparison.

**The structural inventory and split are frozen; prediction inputs are not yet ready.** P1
must curate complete chemical heads/linkers and prepare independently framed inputs. Binary,
apo and sequence-input comparisons need separately traced sources. Failed preparation and
timeouts must remain in the denominator once the execution inventory is fixed. A preparation
fix requires a versioned amendment, preserving groups and explaining changed eligibility.

## Recheck the evidence

Place PDB mmCIF files or `pdb_0000<id>.cif.gz` archives and CCD `<component>.cif` files in a
local data directory, then run:

```bash
uv run topics benchmark-check benchmarks/inventory-v1.json
uv run python docs/verify_benchmark.py --manifest benchmarks/inventory-v1.json \
  --data-dir <data directory> --out <new verification.json>
```

The offline verifier checks source hashes for every record and, for eligible cases, ligand
completeness, scaffold identity, selected biological assemblies and native contact counts.
`verification-v1.json` records the performed verification. Raw downloads are not versioned;
changed upstream files must not silently replace the frozen snapshot.
