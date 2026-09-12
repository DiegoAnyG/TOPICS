"""Batch commands and a small interactive terminal workflow."""

import argparse
import importlib.metadata
import json
import os
from pathlib import Path


def execution_options(parser):
    parser.add_argument("--conformers", type=int, default=16, help="Conformers per seed (default: 16)")
    parser.add_argument("--seeds", type=int, nargs="+", default=[42])
    parser.add_argument("--threads", type=int, default=1, help="CPU threads, default 1; maximum 8")
    parser.add_argument("--device", choices=["cpu", "auto", "cuda"], default="cpu", help="Distance backend; CPU minimizes GPU power use")
    parser.add_argument("--pause", type=float, default=0.05, help="Idle seconds between candidates (default: 0.05)")


def doctor(device="auto"):
    from .geometry import Distances
    backend = Distances(device)
    print(f"Available CPU cores: {os.cpu_count()}; default computation threads: 1")
    print(backend.reason)
    for package in ("topics-ternary", "rdkit", "numpy", "gemmi"):
        try:
            print(f"{package}: {importlib.metadata.version(package)}")
        except importlib.metadata.PackageNotFoundError:
            print(f"{package}: source checkout")


def interactive():
    print("\nTOPICS — ternary assembly and crystallographic controls")
    print("Default: 16 conformers, one CPU thread. Outputs are structural hypotheses.")
    while True:
        print("\n1  Run a crystallographic control\n2  Assemble from an input JSON\n3  Evaluate a saved run\n4  Check compute devices\n0  Exit")
        choice = input("Select: ").strip()
        if choice == "0":
            return 0
        if choice == "4":
            doctor(); continue
        if choice not in {"1", "2", "3"}:
            print("Choose 0, 1, 2, 3 or 4."); continue
        if choice == "3":
            args = ["evaluate", input("Run directory: ").strip(), "--reference", input("Reference NPZ: ").strip()]
        else:
            if choice == "1":
                available = sorted(Path("examples").glob("*.json"))
                if available:
                    print("Available controls: " + ", ".join(str(p) for p in available))
            config = input("Control specification JSON: " if choice == "1" else "Input JSON: ").strip()
            output = input("New output directory: ").strip()
            count = input("Conformers per seed [16]: ").strip() or "16"
            device = input("Distance device cpu/auto/cuda [cpu]: ").strip() or "cpu"
            args = ["benchmark" if choice == "1" else "assemble", config, "--out", output, "--conformers", count, "--device", device]
        try:
            main(args)
        except SystemExit:
            print("Invalid command options; try again.")


def main(argv=None):
    parser = argparse.ArgumentParser(prog="topics", description="Conformer-driven PROTAC assembly and transparent structural validation")
    commands = parser.add_subparsers(dest="command")
    commands.add_parser("interactive", help="Open the terminal menu")
    check = commands.add_parser("doctor", help="Check installed tools and GPU runtime")
    check.add_argument("--device", choices=["cpu", "auto", "cuda"], default="auto")
    inspect = commands.add_parser("inspect", help="List author chains and hetero residues")
    inspect.add_argument("structure")
    for command in ("assemble", "benchmark"):
        sub = commands.add_parser(command)
        sub.add_argument("config")
        sub.add_argument("--out", required=True)
        execution_options(sub)
        if command == "benchmark":
            sub.add_argument("--data-dir", help="Use local PDB/CCD downloads; no network")
    sub = commands.add_parser("evaluate")
    sub.add_argument("run")
    sub.add_argument("--reference", required=True)
    args = parser.parse_args(argv)
    try:
        if args.command in (None, "interactive"):
            return interactive()
        if args.command == "doctor":
            doctor(args.device)
        elif args.command == "inspect":
            import gemmi
            structure = gemmi.read_structure(args.structure)
            print(f"Models: {len(structure)}; resolution: {structure.resolution} Å")
            for chain in structure[0]:
                hetero = [(res.name, str(res.seqid)) for res in chain if res.het_flag == "H"]
                print(f"Author chain {chain.name}: {len(chain)} residues; hetero residues: {hetero}")
        elif args.command == "evaluate":
            from .benchmark import evaluate
            result = evaluate(args.run, args.reference)
            print(json.dumps({k: result[k] for k in ("type", "top_ranked", "best_sampled")}, indent=2))
        else:
            from .engine import run
            output = Path(args.out)
            config = args.config
            if args.command == "benchmark":
                from .benchmark import prepare
                if output.exists() and any(output.iterdir()):
                    raise ValueError("Choose a new benchmark output directory.")
                print("PDB usage and citation: https://www.rcsb.org/pages/usage-policy")
                config = prepare(config, output / "inputs", args.data_dir)
                output = output / "run"
            run(config, output, args.conformers, args.seeds, args.threads, args.device, args.pause)
            if args.command == "benchmark":
                from .benchmark import evaluate
                result = evaluate(output, Path(config).parent / "reference.npz")
                print(f"Retrospective control: top-ranked E3 Cα RMSD {result['top_ranked']['e3_ca_rmsd_A']:.2f} Å; best sampled {result['best_sampled']['e3_ca_rmsd_A']:.2f} Å.")
        return 0
    except (KeyboardInterrupt, EOFError):
        print("\nStopped. Partial output, if present, is not a completed run.")
        return 130
    except (ValueError, OSError, KeyError, RuntimeError) as exc:
        print(f"TOPICS: {exc}")
        return 1
