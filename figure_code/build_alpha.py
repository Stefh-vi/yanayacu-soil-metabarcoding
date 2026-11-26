#!/usr/bin/env python3
"""Rebuild publication alpha figures from frozen sample data and test results.

Public entry point: python figure_code/reproduce_figures.py --output NEW_OUTPUT_DIRECTORY
Requires Python >= 3.10, NumPy and Matplotlib. No statistical fitting is done here.
"""
from __future__ import annotations
import argparse
from paths import QA
import csv
import hashlib
import json
import platform
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.ticker import FuncFormatter, MultipleLocator

ORDER = ["B1", "B2", "B3", "M1", "M2", "M3", "A1", "A2", "A3"]
CATEGORIES = ["Low degradation", "Moderate degradation", "High degradation"]
COLORS = {CATEGORIES[0]: "#287F83", CATEGORIES[1]: "#BD8B18", CATEGORIES[2]: "#904D6B"}
X_OFFSETS = [-0.19, 0, 0.19]
METRICS = [
    ("Observed ASVs", "Observed ASVs", "Observed ASVs", "Observed ASVs 2.5th percentile", "Observed ASVs 97.5th percentile"),
    ("Shannon diversity", "Shannon (ln, nats)", "Shannon diversity (ln, nats)", "Shannon 2.5th percentile", "Shannon 97.5th percentile"),
]

def read_tsv(path):
    with path.open(encoding="utf-8-sig", newline="") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))

def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()

def p_fmt(value):
    value = float(value)
    return "< 0.001" if value < .001 else f"= {value:.3f}"

def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--package", type=Path, required=True)
    args = ap.parse_args()
    package = args.package.resolve()
    out = QA / "alpha_render"
    out.mkdir(parents=True, exist_ok=True)
    inputs = {
        "alpha": package / "data/figure_inputs/alpha/alpha.tsv",
        "subsampling": package / "data/figure_inputs/alpha/subsampling.tsv",
        "tests": package / "data/figure_inputs/alpha/ALPHA_TESTS.tsv",
    }
    alpha = read_tsv(inputs["alpha"])
    technical = {(r["Marker"], r["Sample"]): r for r in read_tsv(inputs["subsampling"])}
    tests = {(r["Marker"], r["Endpoint"]): r for r in read_tsv(inputs["tests"])}
    points = {(r["Marker"], r["Sample"]): r for r in alpha}
    assert len(points) == 18 and len(alpha) == 18
    for marker in ("16S", "ITS"):
        assert {s for m, s in points if m == marker} == set(ORDER)
        for sample in ORDER:
            r = points[(marker, sample)]
            assert r["Locality"] in {"J", "SI"}
            assert int(r["Depth (sequences)"]) == (174661 if marker == "16S" else 450343)
        assert [sum(points[(marker, s)]["Category"] == c for s in ORDER) for c in CATEGORIES] == [3, 3, 3]
        for endpoint, column, _, lo, hi in METRICS:
            t = tests[(marker, endpoint)]
            assert int(t["n"]) == 9 and int(t["category_df"]) == 2 and int(t["residual_df"]) == 5
            assert float(t["P_raw"]) <= float(t["P_Holm"]) <= 1
            for s in ORDER:
                assert np.isfinite(float(points[(marker, s)][column]))

    plt.rcParams.update({
        "font.family": "DejaVu Sans", "font.size": 9,
        "axes.labelsize": 9, "xtick.labelsize": 8.4, "ytick.labelsize": 8.4,
        "axes.linewidth": .75, "pdf.fonttype": 42, "ps.fonttype": 42,
        "svg.fonttype": "none", "svg.hashsalt": "yanayacu-alpha-revision-20260921",
        "figure.facecolor": "white", "savefig.facecolor": "white",
        "text.color": "#242424", "axes.labelcolor": "#242424", "xtick.color": "#242424",
        "ytick.color": "#242424", "mathtext.default": "regular", "mathtext.fontset": "dejavusans",
    })
    source = []
    bounds_audit = []

    def panel(ax, marker, metric, letter):
        endpoint, column, ylabel, qlo, qhi = metric
        t = tests[(marker, endpoint)]
        # Header is deliberately outside the plotting region; all statistics are source driven.
        title = f"{letter}   {marker} | {'Observed ASVs' if endpoint == 'Observed ASVs' else 'Shannon diversity'}"
        ax.text(0, 1.225, title, transform=ax.transAxes, ha="left", va="bottom", fontsize=10, weight="bold")
        estimation = f"{int(t['Depth']):,} sequences; " + ("one rarefaction draw" if marker == "16S" else "mean of 100 draws")
        ax.text(0, 1.154, estimation, transform=ax.transAxes, fontsize=8.0, color="#545454")
        ax.text(0, 1.093, f"$F_{{2,5}}$ = {float(t['F']):.3f}; partial $R^2$ = {float(t['R2_partial']):.3f}",
                transform=ax.transAxes, fontsize=8.3)
        ax.text(0, 1.032, f"$P$ {p_fmt(t['P_raw'])}; Holm $P$ {p_fmt(t['P_Holm'])}", transform=ax.transAxes, fontsize=8.3)
        ax.set_xlim(-.43, 2.43)
        ax.set_xticks(range(3), ["Low\n(n = 3)", "Moderate\n(n = 3)", "High\n(n = 3)"])
        for label, cat in zip(ax.get_xticklabels(), CATEGORIES):
            label.set_color(COLORS[cat])
        ax.set_ylabel(ylabel, labelpad=6)
        ax.tick_params(axis="x", length=0, pad=7)
        ax.tick_params(axis="y", length=3, width=.7)
        ax.spines[["top", "right"]].set_visible(False)
        ax.spines[["left", "bottom"]].set_color("#696969")
        ax.grid(axis="y", color="#E7E7E7", linewidth=.65, zorder=0)
        if endpoint == "Observed ASVs":
            ax.set_ylim((2750, 5250) if marker == "16S" else (300, 1840))
            ax.yaxis.set_major_locator(MultipleLocator(500 if marker == "16S" else 300))
            ax.yaxis.set_major_formatter(FuncFormatter(lambda x, pos: f"{x:,.0f}"))
        else:
            ax.set_ylim((6.7, 8.09) if marker == "16S" else (3.75, 5.97))
            ax.yaxis.set_major_locator(MultipleLocator(.3 if marker == "16S" else .5))
            ax.yaxis.set_major_formatter(FuncFormatter(lambda x, pos: f"{x:.1f}"))
        labels = []
        for ci, cat in enumerate(CATEGORIES):
            cat_samples = [s for s in ORDER if points[(marker, s)]["Category"] == cat]
            median = float(np.median([float(points[(marker, s)][column]) for s in cat_samples]))
            ax.plot([ci-.30, ci+.30], [median, median], color="#5E5E5E", lw=1.45, zorder=2,
                    solid_capstyle="butt")
            for offset, sample in zip(X_OFFSETS, cat_samples):
                r = points[(marker, sample)]
                sr = technical[(marker, sample)]
                x, y = ci+offset, float(r[column])
                low, high = float(sr[qlo]), float(sr[qhi])
                # ITS technical percentile ranges are displayed on their measured scale;
                # ranges smaller than a symbol are not artificially enlarged.
                if marker == "ITS":
                    ax.vlines(x, low, high, color=COLORS[cat], linewidth=.75, zorder=2.5)
                    ax.hlines([low, high], x-.026, x+.026, color=COLORS[cat], linewidth=.65, zorder=2.5)
                ax.scatter([x], [y], c=[COLORS[cat]], marker="o" if r["Locality"] == "J" else "^",
                           s=40 if r["Locality"] == "J" else 48, edgecolors="white", linewidths=.7, zorder=4)
                # Fixed offsets reveal sample identities without changing the measured value.
                # For adjacent B1/B2 values, opposing label alignment avoids a collision.
                dx, dy, ha, va = 0, 7.0, "center", "bottom"
                if marker == "16S" and endpoint == "Observed ASVs" and sample == "B1":
                    dy, va = -7, "top"
                if marker == "ITS" and sample == "B1":
                    dx, dy, ha = -3.5, 7, "right"
                elif marker == "ITS" and sample == "B2":
                    dx, dy, ha = 3.5, 7, "left"
                label = ax.annotate(sample, (x, y), xytext=(dx, dy), textcoords="offset points",
                                    ha=ha, va=va, fontsize=8.0, color="#353535", zorder=5)
                labels.append(label)
                source.append({
                    "Marker": marker, "Endpoint": endpoint, "Sample": sample, "Locality_code": r["Locality"],
                    "Locality": "Jaén" if r["Locality"] == "J" else "San Ignacio", "Category": cat,
                    "Value": y, "Depth_sequences": t["Depth"], "Draws_for_estimator": t["Draws"],
                    "Technical_q025_displayed": low if marker == "ITS" else "",
                    "Technical_q975_displayed": high if marker == "ITS" else "",
                    "Category_median": median, "F_2_5": t["F"], "Partial_R2": t["R2_partial"],
                    "P_raw": t["P_raw"], "P_Holm": t["P_Holm"], "Holm_family": t["Family"],
                    "Horizontal_display_position": x,
                })
        return labels

    def footer(fig, standalone=False):
        handles = [
            Line2D([0], [0], marker="o", color="none", markerfacecolor="#555555", markeredgecolor="white", markersize=6, label="Jaén (J)"),
            Line2D([0], [0], marker="^", color="none", markerfacecolor="#555555", markeredgecolor="white", markersize=6.8, label="San Ignacio (SI)"),
            Line2D([0], [0], color="#5E5E5E", lw=1.45, label="Category median"),
        ]
        fig.legend(handles=handles, loc="lower center", bbox_to_anchor=(.53, .080 if not standalone else .125),
                   ncol=3, frameon=False, fontsize=8.1, handletextpad=.45, columnspacing=1.6, handlelength=1.6)
        fig.text(.53, .058 if not standalone else .085, "Locality-adjusted category tests; permutations restricted within locality.",
                 ha="center", fontsize=7.4, color="#494949")
        if standalone:
            text = "Ranges: 2.5th-97.5th percentiles of 100 rarefaction draws; technical variation only."
        else:
            text = "ITS ranges show technical rarefaction variation, not biological confidence intervals."
        fig.text(.53, .034 if not standalone else .045, text, ha="center", fontsize=7.4, color="#494949")

    def save(fig, stem, labels):
        fig.canvas.draw()
        renderer = fig.canvas.get_renderer()
        collisions = []
        for ax, lbls in labels:
            for i, a in enumerate(lbls):
                for b in lbls[i+1:]:
                    if a.get_window_extent(renderer).overlaps(b.get_window_extent(renderer)):
                        collisions.append([a.get_text(), b.get_text()])
                b = a.get_window_extent(renderer)
                assert fig.bbox.contains(b.x0, b.y0) and fig.bbox.contains(b.x1, b.y1), a.get_text()
        assert not collisions, f"Sample label collisions: {collisions}"
        bounds_audit.append({"figure": stem, "sample_label_overlaps": collisions, "labels_within_page": True,
                             "width_inches": fig.get_figwidth(), "height_inches": fig.get_figheight()})
        metadata = {"Title": stem.replace("Figure_03_", "Figure_02_").replace("_", " "), "Subject": "Alpha diversity of nine soil units; locality-adjusted permutation tests",
                    "CreationDate": None, "ModDate": None}
        fig.savefig(out / f"{stem}.pdf", metadata=metadata)
        fig.savefig(out / f"{stem}.svg", metadata={"Title": metadata["Title"], "Description": metadata["Subject"], "Date": None})
        fig.savefig(out / f"{stem}_600dpi.png", dpi=600)
        fig.savefig(out / f"{stem}_preview.png", dpi=170)
        plt.close(fig)

    fig, axes = plt.subplots(2, 2, figsize=(7.2, 7.45))
    fig.subplots_adjust(left=.105, right=.975, top=.883, bottom=.195, wspace=.38, hspace=.76)
    all_labels = []
    for ri, metric in enumerate(METRICS):
        for ci, marker in enumerate(("16S", "ITS")):
            ax = axes[ri, ci]
            all_labels.append((ax, panel(ax, marker, metric, "abcd"[ri*2+ci])))
    footer(fig)
    save(fig, "Figure_02_Alpha_Diversity_Revised", all_labels)
    # The table contains one row per marker x metric x sample; avoid standalone duplicate rows.
    source = source[:]
    joint_source = source[:]
    fig, axes = plt.subplots(1, 2, figsize=(7.2, 4.2))
    fig.subplots_adjust(left=.105, right=.975, top=.762, bottom=.29, wspace=.38)
    all_labels = [(ax, panel(ax, "ITS", metric, "ab"[i])) for i, (ax, metric) in enumerate(zip(axes, METRICS))]
    footer(fig, True)
    save(fig, "ITS_Alpha_Diversity_Revised", all_labels)
    with (out / "Alpha_Figure_Source_Data.tsv").open("w", encoding="utf-8", newline="") as fh:
        wr = csv.DictWriter(fh, delimiter="\t", fieldnames=list(joint_source[0]))
        wr.writeheader(); wr.writerows(joint_source)
    captions = """# Revised alpha-diversity figures

## Figure 2. Alpha diversity across soil degradation categories

Observed amplicon sequence variant (ASV) richness (a, b) and Shannon diversity using natural logarithms (c, d) for 16S (a, c) and ITS (b, d). Each point represents one of nine soil units (three per degradation category); sample labels retain the original sample identifiers. Colours indicate degradation category, and shapes indicate locality: circles, Jaén (J); triangles, San Ignacio (SI). Horizontal lines indicate category medians. The 16S estimates use one rarefaction realization at 174,661 sequences per sample; ITS estimates are means of 100 rarefaction draws at 450,343 sequences per sample. Thin ITS ranges represent the 2.5th and 97.5th percentiles across technical draws and may be narrower than the plotting symbols. These ranges are not confidence intervals for biological variation.

Category associations were evaluated with partial redundancy analysis of each univariate response, conditioning on locality. The displayed F statistics have 2 and 5 degrees of freedom; partial R² denotes the category fraction of variation remaining after locality adjustment. P values use all 4,319 non-identity permutations within locality, plus the observed configuration. Holm adjustment was applied within the three-endpoint 16S category family and the expanded four-endpoint ITS category family (two alpha and two beta endpoints). The ITS extension was documented before these new tests were executed, after the descriptive patterns had already been inspected. Neither alpha endpoint met a Holm-adjusted P < 0.05 criterion in either marker. Sampling categories and localities are unbalanced, and these associations do not establish an independent causal effect of degradation. Marker-specific richness estimates are not directly comparable measures of total organismal richness.

## ITS alpha-diversity panel

Observed ASV richness (a) and Shannon diversity (b) for the same nine ITS soil units. Values, encodings, rarefaction summaries, locality-adjusted tests and four-endpoint Holm correction are identical to the ITS panels of the revised Figure 2. This panel is an alternative presentation of those results, not an additional analysis or independent dataset.

## Reproduction and interpretation

- Plot inputs: `data/figure_inputs/alpha/alpha.tsv`, `inputs/subsampling.tsv`, and `results/ALPHA_TESTS.tsv`.
- Figure source data: `Alpha_Figure_Source_Data.tsv`, with one row per marker, endpoint and soil unit.
- Each PDF and SVG is vector artwork with text preserved. PNG files ending `_600dpi.png` are raster submission copies; `_preview.png` files are for convenient review.
- The horizontal offsets distinguish samples within a category and have no quantitative interpretation.
- Technical subsamples were not treated as biological replicates. The sample size remains nine units per marker.
- The expanded ITS family changes adjusted beta-diversity P values. Consult the revised category-test table rather than carrying forward adjustment from the historical two-test family.
"""
    (out / "Figure_Captions_EN.md").write_text(captions, encoding="utf-8")
    (out / "FIGURE_QA.json").write_text(json.dumps({"checks": bounds_audit,
        "source_rows": len(joint_source), "unique_marker_metric_sample_rows": len({(r['Marker'],r['Endpoint'],r['Sample']) for r in joint_source}),
        "statistics_read_from_result_file": True,
        "technical_ranges_not_rescaled": True,
        "visual_inspection": "Required after PDF rendering; see delivery audit."}, indent=2), encoding="utf-8")
    manifest = {
        "python": platform.python_version(), "numpy": np.__version__, "matplotlib": matplotlib.__version__,
        "inputs": {k: {"path": str(p.relative_to(package)), "sha256": digest(p)} for k, p in inputs.items()},
        "script_sha256": digest(Path(__file__)),
        "outputs": {p.name: digest(p) for p in sorted(out.iterdir()) if p.is_file() and
            (p.name.startswith(("Figure_02_Alpha_Diversity_Revised", "ITS_Alpha_Diversity_Revised")) or
             p.name in {"Alpha_Figure_Source_Data.tsv", "Figure_Captions_EN.md", "FIGURE_QA.json"})},
    }
    (out / "FIGURE_BUILD_MANIFEST.json").write_text(json.dumps(manifest, indent=2), encoding="utf-8")
    print(json.dumps({"output_directory": str(out), "source_rows": len(joint_source), "figures": 2, "sample_label_collisions": 0}))

if __name__ == "__main__":
    main()

