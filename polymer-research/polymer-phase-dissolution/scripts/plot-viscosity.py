"""Çalıştırma: python plot-viscosity.py --dist dist --output PS_Viskozite_Sinamasi.png"""

import argparse
import csv
import json
from pathlib import Path
import subprocess

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import numpy as np


SYSTEMS = {
    "ps-cyclohexane-theta": ("Siklohekzan · 34,5 °C", "#00796D"),
    "ps-benzene-25": ("Benzen · 25 °C", "#825199"),
}
INK = "#222D38"
GREY = "#536471"


def read_rows(directory, csv_name, script_name, global_name):
    for candidate in (directory, directory / "viscosity-inputs", directory / "intrinsic",
                      directory / "external", directory / "data"):
        path = candidate / csv_name
        if path.is_file():
            with path.open(encoding="utf-8-sig", newline="") as stream:
                return list(csv.DictReader(stream))
    matches = list(directory.rglob(script_name))
    if len(matches) != 1:
        raise ValueError(f"Expected one {script_name} or {csv_name} under {directory}")
    # Execute the two published, project-controlled browser data modules in Node.
    reader = ("global.window={}; require(process.argv[1]);"
              "process.stdout.write(JSON.stringify(window[process.argv[2]].rows))")
    result = subprocess.run(["node", "-e", reader, str(matches[0].resolve()), global_name],
                            check=True, capture_output=True, text=True)
    return json.loads(result.stdout)


def tr(value, digits=1):
    return f"{value:.{digits}f}".replace(".", ",")


def configure_axes(ax):
    ax.spines[["right", "top"]].set_visible(False)
    for side in ("bottom", "left"):
        ax.spines[side].set_color("#88949E")
        ax.spines[side].set_linewidth(.7)
    ax.tick_params(colors=GREY, labelsize=9.5, length=3)
    ax.set_axisbelow(True)


def make_figure(training_rows, external_rows):
    plt.rcParams.update({"font.family": "DejaVu Sans", "text.color": INK,
                         "axes.labelcolor": INK, "axes.unicode_minus": True})
    fig = plt.figure(figsize=(7.0, 8.7), facecolor="white")
    fig.text(.068, .965, "Polistirenin intrinsik viskozitesi", size=20, weight="bold")
    fig.text(.07, .929, "İki çözücüde Mark–Houwink uyumu ve bağımsız ölçümler",
             size=12, color=GREY)
    fig.text(.07, .865, "A   Viskozite ve mol kütlesi", size=13, weight="bold")
    top = fig.add_axes([.16, .535, .80, .285])
    bottom = fig.add_axes([.265, .126, .69, .232])
    configure_axes(top)
    configure_axes(bottom)
    top.grid(which="major", axis="both", linewidth=.6, color="#DBE2E5")
    top.set(xscale="log", yscale="log", xlabel=r"Kütle ortalamalı mol kütlesi $M_w$ (g/mol)",
            ylabel=r"İntrinsik viskozite $[\eta]$ (mL/g)")
    top.set_xlim(4.2e4, 8.3e7)
    top.set_ylim(13, 1.8e4)
    top.tick_params(which="minor", length=1.5, color="#AAB6BC")
    legend = []
    errors = []
    for system, (label, color) in SYSTEMS.items():
        observations = [r for r in training_rows
                        if r["system"] == system and str(r["training"]).lower() == "true"]
        external = [r for r in external_rows if r["system"] == system and r["role"] == "external"]
        mass = np.array([float(r["Mw_g_mol"]) for r in observations])
        eta = np.array([float(r["eta_mL_g"]) for r in observations])
        if len(mass) < 2 or not external:
            raise ValueError(f"Training or independent data missing for {system}")
        a, log_k = np.polyfit(np.log(mass), np.log(eta), 1)

        def predict(x):
            return np.exp(log_k) * x**a

        ext_mass = np.array([float(r["Mw_g_mol"]) for r in external])
        ext_eta = np.array([float(r["eta_mL_g"]) for r in external])
        lower, upper = min(mass), max(mass)
        x_in = np.geomspace(lower, upper, 80)
        top.plot(x_in, predict(x_in), c=color, lw=1.7)
        if min(ext_mass) < lower:
            x_low = np.geomspace(min(ext_mass), lower, 80)
            top.plot(x_low, predict(x_low), c=color, lw=1.4, ls=(0, (4, 3)))
        if max(ext_mass) > upper:
            x_high = np.geomspace(upper, max(ext_mass), 80)
            top.plot(x_high, predict(x_high), c=color, lw=1.4, ls=(0, (4, 3)))
        top.scatter(mass, eta, c=color, s=25, zorder=3)
        top.scatter(ext_mass, ext_eta, edgecolors=color, facecolors="white",
                    s=40, lw=1.5, marker="D", zorder=4)
        legend.append(Line2D([], [], color=color, lw=2, label=label))
        for row, m, observed in zip(external, ext_mass, ext_eta):
            errors.append((row["sample"], system, 100 * (predict(m)/observed - 1), lower <= m <= upper))
    top.legend(handles=legend, frameon=False, loc="upper left", fontsize=9.2,
               handlelength=1.6, borderpad=.2, labelcolor=INK)
    marker_key = [
        Line2D([], [], ls="", marker="o", color=GREY, markersize=5, label="Uyum ölçümü"),
        Line2D([], [], ls="", marker="D", color=GREY, markerfacecolor="white",
               markersize=5.5, label="Bağımsız ölçüm"),
        Line2D([], [], ls=(0, (4, 3)), color=GREY, label="Ekstrapolasyon")
    ]
    fig.legend(handles=marker_key, loc="center", bbox_to_anchor=(.53,.448), ncol=3,
               frameon=False, handlelength=1.8, columnspacing=1.4, fontsize=9.5)
    fig.text(.07, .395, "B   Bağımsız numunelerde tahmin hatası", size=13, weight="bold")
    bottom.spines["left"].set_visible(False)
    bottom.tick_params(axis="y", length=0)
    bottom.grid(axis="x", linewidth=.6, color="#DBE2E5")
    bottom.axvline(0, lw=1, color=GREY)
    errors.sort(key=lambda e: (e[1] == "ps-benzene-25", e[0]))
    ylabels = []
    for y, (sample, system, error, in_range) in enumerate(errors):
        _, color = SYSTEMS[system]
        bottom.barh(y, error, height=.56, color=color if in_range else "white",
                    edgecolor=color, lw=1.3, zorder=2)
        bottom.text(error + (.75 if error >= 0 else -.75), y,
                    ("+" if error >= 0 else "−") + tr(abs(error)),
                    fontsize=9.2, color=color, va="center",
                    ha="left" if error >= 0 else "right", weight="bold" if in_range else "normal")
        ylabels.append(sample + (" · benzen" if system == "ps-benzene-25" else ""))
    bottom.set_yticks(range(len(errors)), ylabels, fontsize=9.5)
    bottom.set_xlim(-10.5, 42)
    bottom.set_xticks([-10, 0, 10, 20, 30, 40])
    bottom.set_ylim(len(errors)-.45, -.6)
    bottom.set_xlabel("100 × (tahmin − ölçüm) / ölçüm  (%)", size=10.0)
    fig.text(.07, .053, "Dolu çubuk kalibrasyon aralığındaki numuneyi gösterir.", fontsize=9.3, color=GREY)
    fig.text(.07, .030, "Kaynaklar  Einaga vd. (1979), Baysal ve Uyanık (1992), NBS SRM 705 (1978).",
             fontsize=8.5, color=GREY)
    return fig


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--dist", type=Path, required=True, help="Dağıtım paketi veya veri dosyaları dizini")
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    fit = read_rows(args.dist, "viscosity-intrinsic.csv", "viscosity-intrinsic-data.js", "VIS_INTRINSIC")
    test = read_rows(args.dist, "external-intrinsic.csv", "viscosity-external-data.js", "VIS_EXTERNAL")
    figure = make_figure(fit, test)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    figure.savefig(args.output, dpi=200, facecolor="white")
    plt.close(figure)


if __name__ == "__main__":
    main()
