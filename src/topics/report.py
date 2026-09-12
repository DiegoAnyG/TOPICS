"""Standalone interactive reports and vector/600-dpi scientific figures."""

import csv
import html
import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import plotly.graph_objects as go


COLORS = {"POI": "#0072B2", "E3": "#D55E00", "PROTAC": "#009E73", "Crystal E3": "#777777"}


def render_assessment(root, result):
    """Export every assessed candidate, with missing metrics retained as gaps."""
    root = Path(root)
    rows = result["candidates"]
    ranks = [row["rank"] for row in rows]
    columns = [("dockq", "PPI DockQ (higher is better)"),
               ("ligand_symmetry_rmsd_A", "Ligand RMSD after POI fit (Å)"),
               ("e3_ca_rmsd_A", "E3 Cα RMSD after POI fit (Å)")]
    sections = []
    figure, axes = plt.subplots(1, len(columns), figsize=(12, 3.7), layout="constrained")
    for index, ((metric, label), axis) in enumerate(zip(columns, axes)):
        values = [row.get(metric) for row in rows]
        colors = ["#009E73" if row["steric_valid"] else "#D55E00" for row in rows]
        chart = go.Figure(go.Scatter(x=ranks, y=values, mode="markers", text=[row["candidate"] for row in rows],
                                    marker={"color": colors}, hovertemplate="%{text}<br>Rank %{x}<br>%{y:.4f}<extra></extra>"))
        chart.update_layout(template="plotly_white", xaxis_title="Original candidate rank", yaxis_title=label,
                            title=label, height=350)
        sections.append(chart.to_html(full_html=False, include_plotlyjs=True if index == 0 else False))
        mask = [i for i, value in enumerate(values) if value is not None]
        axis.scatter([ranks[i] for i in mask], [values[i] for i in mask], c=[colors[i] for i in mask], s=12)
        if metric == "dockq":
            axis.axhline(0.23, color="#555555", linestyle="--", linewidth=0.8)
            axis.set_ylim(0, 1)
        axis.set(xlabel="Original candidate rank", ylabel=label)
    figure.suptitle("Independent assessment; orange = typed steric violations, green = none")
    for suffix in ("svg", "pdf", "png"):
        figure.savefig(root / ("assessment." + suffix), dpi=600)
    plt.close(figure)
    fields = ["candidate", "rank", "dockq", "ligand_symmetry_rmsd_A", "e3_ca_rmsd_A", "chemical_valid", "steric_valid", "joint_success", "errors"]
    table = "<table><thead><tr>" + "".join("<th>" + html.escape(f) + "</th>" for f in fields) + "</tr></thead><tbody>"
    for row in rows:
        table += "<tr>" + "".join("<td>" + html.escape(str(row.get(f)) if row.get(f) is not None else "Not assessed") + "</td>" for f in fields) + "</tr>"
    table += "</tbody></table>"
    details = {k: v for k, v in result.items() if k != "candidates"}
    content = ("<!doctype html><html lang='en'><meta charset='utf-8'><meta name='viewport' content='width=device-width'>"
               "<title>TOPICS independent assessment</title><style>body{font:16px system-ui;margin:2rem;max-width:1200px}"
               "table{border-collapse:collapse;font-size:13px}td,th{padding:.4rem;border:1px solid #ddd}"
               "pre{white-space:pre-wrap;overflow-wrap:anywhere}.data{overflow:auto}</style><body>"
               "<h1>TOPICS independent assessment</h1><p>All candidates retain their original rank. "
               "Chemical validity is the mol_fast ligand-geometry subset plus defined stereochemistry. "
               "Joint success also requires explicit complete heads. Missing results are not passes.</p>"
               "<p><a href='assessment.csv'>Source data (CSV)</a> · <a href='assessment.json'>Provenance (JSON)</a> · "
               "<a href='assessment.svg'>SVG</a> · <a href='assessment.pdf'>PDF</a> · <a href='assessment.png'>600 dpi PNG</a></p>"
               + "".join(sections) + "<div class='data'>" + table + "</div><details><summary>Protocol and provenance</summary><pre>"
               + html.escape(json.dumps(details, indent=2)) + "</pre></details></body></html>")
    (root / "report.html").write_text(content)


def render(run_dir):
    root = Path(run_dir)
    manifest = json.loads((root / "manifest.json").read_text())
    display_manifest = {k: v for k, v in manifest.items() if k != "status"}
    rows = list(csv.DictReader((root / "candidates.csv").open()))
    arrays = np.load(root / "ensemble.npz", allow_pickle=False)
    evaluation = json.loads((root / "evaluation.json").read_text()) if (root / "evaluation.json").exists() else None
    rank = [int(r["rank"]) for r in rows]
    metric = ([r["e3_ca_rmsd_A"] for r in evaluation["candidates"]] if evaluation else [float(r["head_rmsd_A"]) for r in rows])
    ylabel = "E3 Cα RMSD after POI alignment (Å)" if evaluation else "Maximum binding-head fit RMSD (Å)"
    metric_fig = go.Figure(go.Scatter(x=rank, y=metric, mode="markers", text=[r["candidate"] for r in rows],
                                    marker={"color": [int(r["clash_pairs"]) for r in rows], "colorscale": "Cividis", "showscale": True,
                                            "colorbar": {"title": "Clash pairs"}},
                                    hovertemplate="%{text}<br>Rank %{x}<br>RMSD %{y:.3f} Å<extra></extra>"))
    metric_fig.update_layout(template="plotly_white", xaxis_title="Geometry rank", yaxis_title=ylabel,
                             title="Sampling and ranking assessment", height=430)
    viewer = go.Figure()
    poi_mask = np.array([k[-1] == "CA" for k in manifest["atom_keys"]["poi"]])
    e3_mask = np.array([k[-1] == "CA" for k in manifest["atom_keys"]["e3"]])
    for name, coords in [("POI", arrays["poi"][poi_mask]), ("E3", arrays["e3"][0][e3_mask]), ("PROTAC", arrays["ligand"][0])]:
        viewer.add_trace(go.Scatter3d(x=coords[:, 0], y=coords[:, 1], z=coords[:, 2], mode="markers", name=name,
                                     marker={"size": 3 if name == "PROTAC" else 2, "color": COLORS[name]},
                                     hovertemplate=f"{name}<br>%{{x:.2f}}, %{{y:.2f}}, %{{z:.2f}} Å<extra></extra>"))
    if evaluation:
        ref = np.load(root / "reference_view.npz", allow_pickle=False)["e3"][e3_mask]
        viewer.add_trace(go.Scatter3d(x=ref[:, 0], y=ref[:, 1], z=ref[:, 2], mode="markers", name="Crystal E3",
                                     marker={"size": 2, "color": COLORS["Crystal E3"], "opacity": 0.35}))
    buttons = []
    for i, row in enumerate(rows[:20]):
        coords = [arrays["e3"][i][e3_mask], arrays["ligand"][i]]
        buttons.append({"label": f"Rank {i + 1}: {row['candidate']}", "method": "restyle",
                        "args": [{axis: [c[:, j].tolist() for c in coords] for j, axis in enumerate("xyz")}, [1, 2]]})
    viewer.update_layout(template="plotly_white", height=540, title="Rotate, zoom and select a candidate (Cα and ligand atoms)",
                         scene={"aspectmode": "data", "xaxis_title": "x (Å)", "yaxis_title": "y (Å)", "zaxis_title": "z (Å)"},
                         updatemenus=[{"buttons": buttons, "direction": "down", "x": 0, "y": 1.1}])
    with plt.rc_context({"font.family": "DejaVu Sans", "font.size": 9, "svg.fonttype": "none", "pdf.fonttype": 42}):
        fig, ax = plt.subplots(figsize=(6.5, 3.8), layout="constrained")
        points = ax.scatter(rank, metric, c=[int(r["clash_pairs"]) for r in rows], cmap="cividis", s=25, edgecolors="black", linewidths=0.3)
        ax.set(xlabel="Geometry rank", ylabel=ylabel)
        ax.spines[["top", "right"]].set_visible(False)
        fig.colorbar(points, ax=ax, label="Heavy-atom clash pairs (<2 Å)")
        for suffix in ("svg", "pdf", "png"):
            fig.savefig(root / f"assessment.{suffix}", dpi=600)
        plt.close(fig)
        fig = plt.figure(figsize=(6.5, 4.5), layout="constrained")
        ax = fig.add_subplot(projection="3d")
        for name, coords in [("POI", arrays["poi"][poi_mask]), ("E3", arrays["e3"][0][e3_mask]), ("PROTAC", arrays["ligand"][0])]:
            ax.scatter(*coords.T, s=4, color=COLORS[name], label=name)
        if evaluation:
            ax.scatter(*ref.T, s=3, color=COLORS["Crystal E3"], alpha=0.35, label="Crystal E3")
        all_xyz = np.concatenate([arrays["poi"], arrays["e3"][0], arrays["ligand"][0]] + ([ref] if evaluation else []))
        ax.set_box_aspect(np.maximum(np.ptp(all_xyz, axis=0), 1))
        ax.set(xlabel="x (Å)", ylabel="y (Å)", zlabel="z (Å)")
        ax.legend(loc="upper left", fontsize=8)
        for suffix in ("svg", "pdf", "png"):
            fig.savefig(root / f"assembly.{suffix}", dpi=600)
        plt.close(fig)
    headers = list(rows[0])
    table = "<table><thead><tr>" + "".join(f"<th>{html.escape(k)}</th>" for k in headers) + "</tr></thead><tbody>"
    table += "".join("<tr>" + "".join(f"<td>{html.escape(str(row[k]))}</td>" for k in headers) + "</tr>" for row in rows) + "</tbody></table>"
    assessment = "No reference supplied. Predictive accuracy has not been measured."
    if evaluation:
        top, best = evaluation["top_ranked"], evaluation["best_sampled"]
        assessment = (f"Retrospective bound-component control. Top-ranked E3 Cα RMSD: {top['e3_ca_rmsd_A']:.2f} Å; "
                      f"best sampled: {best['e3_ca_rmsd_A']:.2f} Å (rank {best['rank']}). "
                      f"Top-ranked native PPI contact recall: {100 * top['native_contact_recall']:.1f}%. "
                      "Reference error is evaluated after ranking. This is not an independent blind benchmark.")
    source_url = manifest["provenance"].get("source", {}).get("url", "")
    source = f'<a href="{html.escape(source_url, quote=True)}">PDB source</a>' if source_url.startswith("https://www.rcsb.org/") else "User-supplied inputs"
    methods = f"{manifest['method']}. {manifest['parameters']}. {manifest['limitations']}"
    text = f"""<!doctype html><html lang="en"><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1">
<title>{html.escape(manifest['title'])} | TOPICS</title><style>
body{{font:16px system-ui,sans-serif;max-width:1200px;margin:30px auto;padding:0 20px;color:#183042;background:#fafcfd}}
h1{{font-size:30px}}.notice{{padding:20px;background:#eaf2f8;border-left:4px solid #0072b2}}table{{border-collapse:collapse;font-size:12px}}
td,th{{padding:8px;border-bottom:1px solid #ddd;text-align:right}}.table{{overflow:auto;max-height:420px}}a{{color:#006198}}pre{{white-space:pre-wrap}}
</style><h1>TOPICS · {html.escape(manifest['title'])}</h1><p>{source} · {len(rows)} candidates · {manifest['compute_seconds']:.2f} s compute · {html.escape(manifest['backend_note'])}</p>
<p class="notice">{html.escape(assessment)}</p>
{metric_fig.to_html(full_html=False, include_plotlyjs=True)}
{viewer.to_html(full_html=False, include_plotlyjs=False)}
<h2>Candidate table</h2><p>Ranking prioritizes fewer clashes, lower head-fit error, more PPI contacts, then lower conformer energy. Energy is not binding affinity.</p>
<div class="table">{table}</div><h2>Downloads</h2><p><a href="candidates.csv">Candidates CSV</a> · <a href="manifest.json">Manifest</a> · <a href="best.cif">Best assembly</a> · <a href="ensemble.sdf">Conformers SDF</a> · <a href="view.pml">PyMOL script</a></p>
<p><a href="assessment.svg">Assessment SVG</a> · <a href="assessment.pdf">PDF</a> · <a href="assessment.png">600-dpi PNG</a> · <a href="assembly.svg">Assembly SVG</a> · <a href="assembly.pdf">PDF</a> · <a href="assembly.png">600-dpi PNG</a></p>
<h2>Methods and limitations</h2><p>{html.escape(methods)}</p><details><summary>Reproducibility manifest</summary><pre>{html.escape(json.dumps(display_manifest, indent=2))}</pre></details>
<p>Figures are exportable scientific artifacts; interpretation and journal formatting require author review.</p></html>"""
    (root / "report.html").write_text(text)
    (root / "figure_captions.txt").write_text(
        f"Figure 1. {ylabel} versus geometry rank for {len(rows)} generated candidates. Color denotes heavy-atom pairs below 2 Å. "
        + ("Reference metrics were calculated after candidate selection. " if evaluation else "Head fit is an internal geometry measure, not predictive accuracy. ")
        + "\nFigure 2. Top-ranked rigid assembly: POI (blue), E3 (orange), PROTAC (green). Protein positions are represented by alpha carbons. "
        + ("The reference E3 is shown in gray after POI alignment. " if evaluation else "")
        + "\n" + assessment + "\n")
