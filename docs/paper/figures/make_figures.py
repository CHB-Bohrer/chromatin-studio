"""Draw the paper's figures from docs/paper/data/. Only stale figures are redrawn.

    python docs/paper/figures/make_figures.py            # redraw what changed
    python docs/paper/figures/make_figures.py --force    # redraw everything
    python docs/paper/figures/make_figures.py val-wlc    # redraw named figures

Each figure is one function registered with @figure(name, inputs). A figure is
stale when its SVG is missing, or when the hash of its input files plus its own
drawing code (and the shared style) differs from the one recorded in
.figure-cache.json. Nothing here runs a simulation.

Style: three categorical hues at most (blue, orange, aqua: the set that stays
distinct under colour-vision deficiency), engine output as filled markers,
closed-form references as dashed grey lines, measurements as hollow diamonds,
accepted ranges as pale grey bands. One y-axis per panel.
"""

import hashlib
import inspect
import json
import math
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib.patches import FancyBboxPatch  # noqa: E402

HERE = Path(__file__).resolve().parent
DATA = HERE.parent / "data"
CACHE = HERE / ".figure-cache.json"

BLUE, ORANGE, AQUA = "#257de0", "#f47230", "#00b894"
HUES = [BLUE, ORANGE, AQUA]
INK, MUTED, GRID, REF, BAND = "#191919", "#666666", "#e6e6e6", "#7a7a7a", "#ececec"

plt.rcParams.update({
    "font.family": "DejaVu Sans", "font.size": 9, "axes.titlesize": 9.5, "axes.titleweight": "bold",
    "axes.labelsize": 9, "axes.edgecolor": "#bdbdbd", "axes.labelcolor": INK, "xtick.color": MUTED,
    "ytick.color": MUTED, "axes.spines.top": False, "axes.spines.right": False, "axes.grid": True,
    "grid.color": GRID, "grid.linewidth": 0.6, "legend.frameon": False, "legend.fontsize": 8,
    "lines.linewidth": 2, "lines.markersize": 5, "svg.fonttype": "path", "figure.dpi": 100,
})

REGISTRY = {}


def figure(name, inputs=()):
    def wrap(fn):
        REGISTRY[name] = (fn, list(inputs))
        return fn
    return wrap


def data(name):
    return json.loads((DATA / name).read_text())


def save(fig, name):
    fig.savefig(HERE / f"{name}.svg", bbox_inches="tight")
    plt.close(fig)


def model(ax, x, y, err=None, color=BLUE, label="Engine", marker="o"):
    ax.errorbar(x, y, yerr=err, fmt=marker, color=color, ms=7.5, mew=1.2, mec="white", elinewidth=1.2,
                capsize=0, label=label, zorder=3)


def reference(ax, x, y, label="Closed form", color=REF):
    ax.plot(x, y, "--", color=color, lw=1.6, label=label, zorder=2)


def measured(ax, x, y, err=None, label="Measured"):
    ax.errorbar(x, y, yerr=err, fmt="D", mfc="white", mec=INK, color=INK, ms=6, mew=1.3, elinewidth=1.1,
                capsize=2, label=label, zorder=4)


# ---------------------------------------------------------------- generic report renderer

def draw_report_figure(ax, spec, legend=True):
    """Render one figure description written by a component study (core.figure)."""
    groups = {}
    for s in spec["series"]:
        x = [v for v in s["x"]]
        y = s["y"]
        if s["role"] == "reference":
            ax.plot(x, y, "--", color=REF, lw=1.6, label=s["label"], zorder=2)
        elif s["role"] == "experiment":
            if s["style"] in ("line", "connected"):
                ax.plot(x, y, "-", color="#b5b5b5", lw=3, alpha=0.8, label=s["label"], zorder=1)
            else:
                measured(ax, x, y, label=s["label"])
        else:
            key = s.get("group") or s["label"]
            color = groups.setdefault(key, HUES[len(groups) % 3])
            err = None
            if s.get("low") and s.get("high"):
                err = [[max(v - lo, 0) for v, lo in zip(y, s["low"])], [max(hi - v, 0) for v, hi in zip(y, s["high"])]]
            if s["style"] in ("line", "step"):
                ax.plot(x, y, "-", color=color, lw=2, label=s["label"])
            else:
                model(ax, x, y, err, color=color, label=s["label"])
                if s["style"] == "connected":
                    ax.plot(x, y, "-", color=color, lw=1.2, alpha=0.6)
    for b in spec.get("bands") or []:
        if b.get("y"):
            lo, hi = b["y"]
            ax.axhspan(lo if lo is not None else ax.get_ylim()[0], hi if hi is not None else ax.get_ylim()[1],
                       color=BAND, zorder=0, label=b["label"])
    ax.set_xscale(spec.get("x_scale") or "linear")
    ax.set_yscale(spec.get("y_scale") or "linear")
    ax.set_xlabel(spec["x_label"])
    ax.set_ylabel(spec["y_label"])
    ax.set_title(spec["title"], loc="left")
    if legend:
        ax.legend(loc="best")


def report_panels(name, study, figure_ids, size):
    d = data(f"component-{study}.json")
    specs = {f["id"]: f for f in d.get("figures") or []}
    fig, axes = plt.subplots(1, len(figure_ids), figsize=size)
    axes = [axes] if len(figure_ids) == 1 else axes
    for ax, fid in zip(axes, figure_ids):
        draw_report_figure(ax, specs[fid])
    fig.tight_layout()
    save(fig, name)


def checks_panel(name, study, title, labels, size=(6.4, 2.6)):
    """Observed / expected for each analytical check, with the accepted window."""
    d = data(f"component-{study}.json")
    checks = [c for c in d["checks"] if c.get("expected")]
    fig, ax = plt.subplots(figsize=size)
    ys = list(range(len(checks)))[::-1]
    for y, c in zip(ys, checks):
        ratio = c["observed"] / c["expected"]
        err = (c.get("sem") or 0) / abs(c["expected"])
        allowed = (c.get("allowed") or 0) / abs(c["expected"])
        ax.add_patch(plt.Rectangle((1 - allowed, y - 0.3), 2 * allowed, 0.6, color=BAND, zorder=0, lw=0))
        ax.errorbar(ratio, y, xerr=err, fmt="o", color=BLUE if c["passed"] else ORANGE, ms=6, mec="white",
                    mew=1.5, elinewidth=1.4, zorder=3)
    ax.axvline(1, ls="--", color=REF, lw=1.4)
    ax.set_yticks(ys)
    ax.set_yticklabels([labels.get(c["id"], c["label"]) for c in checks])
    ax.set_xlabel("engine / closed form   (grey: accepted window)")
    ax.set_title(title, loc="left")
    ax.grid(axis="y", visible=False)
    fig.tight_layout()
    save(fig, name)


# ---------------------------------------------------------------- schematics

def box(ax, x, y, w, h, text, color="#f3f6fb", edge=BLUE, size=8.5, weight="normal"):
    ax.add_patch(FancyBboxPatch((x, y), w, h, boxstyle="round,pad=0.02,rounding_size=0.08", fc=color, ec=edge,
                                lw=1.3))
    ax.text(x + w / 2, y + h / 2, text, ha="center", va="center", fontsize=size, color=INK, weight=weight,
            wrap=True)


def arrow(ax, a, b, color=MUTED, style="-|>", rad=0.0):
    ax.annotate("", xy=b, xytext=a, arrowprops=dict(arrowstyle=style, color=color, lw=1.3,
                                                    connectionstyle=f"arc3,rad={rad}", shrinkA=2, shrinkB=2))


@figure("fig-feedback")
def fig_feedback():
    fig, ax = plt.subplots(figsize=(7.2, 3.4))
    ax.set_xlim(0, 10)
    ax.set_ylim(-0.1, 4.9)
    ax.axis("off")
    W, H = 1.9, 0.75
    nodes = {
        "marks": (0.2, 3.5, "Histone marks\n(§4)"), "nuc": (0.2, 2.1, "Nucleosomes &\nremodelers (§4)"),
        "tf": (2.7, 2.8, "Protein binding\n& recruitment (§3)"), "reg": (5.2, 2.8, "Promoter ON &\ninitiation (§7)"),
        "pol": (7.8, 2.8, "Polymerases\n& RNA (§7)"), "poly": (5.2, 0.5, "Chromatin folding\n& contacts (§2)"),
        "coh": (2.7, 0.5, "Cohesin loops\n(§5)"), "env": (7.8, 0.5, "Envelope, tethers,\nbodies (§6)"),
    }
    for key, (x, y, t) in nodes.items():
        box(ax, x, y, W, H, t, edge=ORANGE if key in ("marks", "nuc") else BLUE)
    c = {k: (x + W / 2, y + H / 2) for k, (x, y, _) in nodes.items()}
    e = lambda k, dx=0, dy=0: (c[k][0] + dx, c[k][1] + dy)  # noqa: E731
    arrow(ax, e("marks", W / 2), e("tf", -W / 2, 0.2))
    arrow(ax, e("nuc", W / 2), e("tf", -W / 2, -0.2))
    arrow(ax, e("tf", W / 2), e("reg", -W / 2))
    arrow(ax, e("reg", W / 2), e("pol", -W / 2))
    arrow(ax, e("tf", 0, H / 2), e("marks", W / 2, 0.2), rad=0.3)
    arrow(ax, e("pol", 0, -H / 2), e("poly", W / 2, 0.2), rad=-0.2)
    arrow(ax, e("pol", -0.3, H / 2), e("nuc", 0.3, H / 2), rad=0.35)
    arrow(ax, e("coh", W / 2), e("poly", -W / 2))
    arrow(ax, e("env", -W / 2), e("poly", W / 2, -0.15))
    arrow(ax, e("poly", 0, H / 2), e("reg", 0, -H / 2))
    arrow(ax, e("poly", -0.4, H / 2), e("tf", 0.4, -H / 2))
    ax.text(5.0, 0.05, "Orange: states on the genome    Blue: molecules, rates and 3D structure", ha="center",
            fontsize=8, color=MUTED)
    save(fig, "fig-feedback")


@figure("fig-step")
def fig_step():
    fig, ax = plt.subplots(figsize=(7.2, 1.9))
    ax.set_xlim(0, 10.6)
    ax.set_ylim(0, 2)
    ax.axis("off")
    box(ax, 0.1, 0.55, 1.6, 0.9, "Move particles\n(Brownian, §2)\n× sub-steps", color="#fff4ee", edge=ORANGE,
        size=7.6)
    steps = ["Bind &\nrelease", "Nucleosome\nturnover", "Sliding", "Marks", "Cohesin", "Tethers", "Genes\n& RNA"]
    x = 2.1
    for i, s in enumerate(steps):
        box(ax, x, 0.6, 1.0, 0.8, s, size=7.4)
        if i:
            arrow(ax, (x - 0.13, 1.0), (x, 1.0))
        x += 1.15
    arrow(ax, (1.72, 1.0), (2.1, 1.0))
    ax.annotate("", xy=(0.9, 0.5), xytext=(x - 0.6, 0.5),
                arrowprops=dict(arrowstyle="-|>", color=MUTED, lw=1.2, connectionstyle="arc3,rad=0.25"))
    ax.text(5.6, 0.02, "next step: t → t + Δt", ha="center", fontsize=8, color=MUTED)
    ax.text(6.1, 1.62, "reaction families, rates frozen over Δt", ha="center", fontsize=8, color=MUTED)
    save(fig, "fig-step")


@figure("fig-potentials")
def fig_potentials():
    import numpy as np
    fig, axes = plt.subplots(1, 3, figsize=(7.2, 2.2))
    r = np.linspace(0, 1.6, 200)
    eps, sigma = 8.0, 1.0
    axes[0].plot(r, np.where(r < sigma, eps / 2 * (1 - r ** 2 / sigma ** 2) ** 2, 0), color=BLUE)
    axes[0].set_title("Soft repulsion", loc="left")
    axes[0].set_xlabel("distance r / σ")
    axes[0].set_ylabel("energy (kT)")
    R, depth = 1.5, 1.0
    r2 = np.linspace(0, 2, 200)
    axes[1].plot(r2, np.where(r2 < R, -depth * (1 - r2 ** 2 / R ** 2) ** 2, 0), color=BLUE)
    axes[1].set_title("Mark attraction (optional)", loc="left")
    axes[1].set_xlabel("distance r / bond length")
    rr = np.linspace(0, 1.4, 200)
    axes[2].plot(rr, np.where(rr > 1, 0.5 * 30 * (rr - 1) ** 2, 0), color=BLUE)
    axes[2].set_title("Confining wall", loc="left")
    axes[2].set_xlabel("radius |x| / R_c  (R_c = 1)")
    fig.tight_layout()
    save(fig, "fig-potentials")


# ---------------------------------------------------------------- Section 2: polymer

@figure("val-rouse", ["benchmark-reference.json"])
def val_rouse():
    d = data("benchmark-reference.json")
    curves = {c["id"]: c for c in d["curves"]}
    fig, ax = plt.subplots(figsize=(4.6, 3.2))
    for cid, color, label in [("free", BLUE, "free particle"), ("rouse", ORANGE, "chain bead (internal)"),
                              ("rouse_com", AQUA, "chain centre of mass")]:
        pts = curves[cid]["points"]
        x = [p["lag"] for p in pts]
        model(ax, x, [p["observed"] for p in pts], [p["standard_error"] for p in pts], color=color, label=label)
        ax.plot(x, [p["discrete_exact"] for p in pts], "--", color=REF, lw=1.3)
    ax.plot([], [], "--", color=REF, label="exact (discrete time)")
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlabel("lag time (reduced)")
    ax.set_ylabel("mean-square displacement")
    ax.set_title("Brownian motion of particles and a Rouse chain", loc="left")
    ax.legend(loc="upper left")
    fig.tight_layout()
    save(fig, "val-rouse")


@figure("val-ideal-chain", ["extra-checks.json"])
def val_ideal_chain():
    d = data("extra-checks.json")["ideal_chain"]
    s = d["s"]
    fig, axes = plt.subplots(1, 2, figsize=(7.2, 2.9))
    ax = axes[0]
    model(ax, s, d["r2"], d["r2_sem"], label="Engine")
    reference(ax, s, d["r2_exact"], label="s ⟨b²⟩ (exact)")
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlabel("separation along chain, s (beads)")
    ax.set_ylabel("⟨R²(s)⟩ (reduced length²)")
    ax.set_title("Distance grows linearly with s", loc="left")
    ax.legend(loc="upper left")
    ax = axes[1]
    keep = [i for i, v in enumerate(d["contact"]) if v > 0]
    model(ax, [s[i] for i in keep], [d["contact"][i] for i in keep], [d["contact_sem"][i] for i in keep],
          label="Engine")
    reference(ax, s, d["contact_gaussian"], label="Gaussian chain")
    x0 = 16
    y0 = d["contact_gaussian"][x0 - 1]
    ax.plot([4, 63], [y0 * (4 / x0) ** -1.5, y0 * (63 / x0) ** -1.5], ":", color=MUTED, lw=1.2,
            label="slope −3/2")
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlabel("separation along chain, s (beads)")
    ax.set_ylabel("contact probability (r < 1)")
    ax.set_title("Contact probability P(s)", loc="left")
    ax.legend(loc="lower left")
    fig.tight_layout()
    save(fig, "val-ideal-chain")


@figure("val-wlc", ["component-polymer_wlc.json"])
def val_wlc():
    report_panels("val-wlc", "polymer_wlc", ["tangent_decay", "internal_distance"], (7.2, 3.0))


# ---------------------------------------------------------------- Section 3: proteins

@figure("val-diffusion", ["component-protein_diffusion.json"])
def val_diffusion():
    d = data("component-protein_diffusion.json")
    ids = [f["id"] for f in d.get("figures") or []][:2]
    report_panels("val-diffusion", "protein_diffusion", ids, (3.8 * len(ids), 3.0))


def relaxation(name, block, key, title, ylabel, size=(4.4, 3.0)):
    d = data("extra-checks.json")[block]
    fig, ax = plt.subplots(figsize=size)
    reference(ax, d["grid"], d["theory"])
    model(ax, d["times"], d[key], d[key + "_sem"], label="Engine, mean ± SEM")
    ax.set_xlabel("reduced time")
    ax.set_ylabel(ylabel)
    ax.set_title(title, loc="left")
    ax.legend(loc="best")
    fig.tight_layout()
    save(fig, name)


@figure("val-binding", ["extra-checks.json"])
def val_binding():
    relaxation("val-binding", "binding", "bound", "Binding relaxes to K/(K + k_off)", "fraction of molecules bound")


# ---------------------------------------------------------------- Section 4: marks and nucleosomes

@figure("val-turnover", ["component-nucleosome_turnover.json"])
def val_turnover():
    report_panels("val-turnover", "nucleosome_turnover", ["relaxation"], (4.4, 3.0))


@figure("val-marks", ["extra-checks.json"])
def val_marks():
    relaxation("val-marks", "marks", "marked", "Mark decay follows exp(−d t)", "fraction of beads marked")


# ---------------------------------------------------------------- Section 5: cohesin

@figure("val-cohesin", ["component-cohesin_extrusion.json"])
def val_cohesin():
    checks_panel("val-cohesin", "cohesin_extrusion", "Cohesin: engine against closed form", {})


# ---------------------------------------------------------------- Section 6: tethers

@figure("val-tethers", ["extra-checks.json"])
def val_tethers():
    relaxation("val-tethers", "tethers", "attached", "Tether attachment, eligible beads",
               "fraction of eligible beads attached")


# ---------------------------------------------------------------- Section 7: transcription

@figure("val-telegraph", ["component-telegraph_transcription.json"])
def val_telegraph():
    checks_panel("val-telegraph", "telegraph_transcription", "Bursting: engine against the telegraph model",
                 {"count_on_fraction": "ON fraction (count mode)", "count_mean": "mean RNA (count mode)",
                  "count_fano": "Fano factor (count mode)", "explicit_on_fraction": "ON fraction (individual RNA)",
                  "explicit_mean": "mean RNA (individual RNA)", "explicit_fano": "Fano factor (individual RNA)"})


# ---------------------------------------------------------------- Sections 10-11: experiments, calibration

@figure("exp-donovan", ["recorded-results.json"])
def exp_donovan():
    d = data("recorded-results.json")["donovan"]["factors"]
    fig, ax = plt.subplots(figsize=(5.6, 3.1))
    for i, f in enumerate(d):
        lo, hi = f["measured"]
        ax.add_patch(plt.Rectangle((i - 0.32, lo), 0.64, (hi or 1e4) - lo, color=BAND, zorder=0, lw=0))
        ax.plot([i - 0.32, i + 0.32], [lo, lo], color="#9a9a9a", lw=1)
        if hi:
            ax.plot([i - 0.32, i + 0.32], [hi, hi], color="#9a9a9a", lw=1)
        ax.errorbar(i - 0.1, f["occlusion"][0], yerr=f["occlusion"][1], fmt="o", color=BLUE, ms=7, mec="white",
                    mew=1.5, label="Blocking only" if i == 0 else None, zorder=3)
        if f["release"]:
            ax.errorbar(i + 0.1, f["release"][0], yerr=f["release"][1], fmt="s", color=ORANGE, ms=7,
                        mec="white", mew=1.5, label="Blocking + slower release" if i == 0 else None, zorder=3)
    ax.add_patch(plt.Rectangle((0, 0), 0, 0, color=BAND, label="Measured range"))
    ax.axhline(1, color=REF, ls="--", lw=1.2)
    ax.set_yscale("log")
    ax.set_ylim(0.1, 1000)
    ax.set_xlim(-0.6, len(d) - 0.4)
    ax.set_xticks(range(len(d)))
    ax.set_xticklabels([f["label"] for f in d])
    ax.set_ylabel("affinity ratio, nucleosome / DNA")
    ax.set_title("Nucleosome effect on affinity (Donovan et al.)", loc="left")
    ax.grid(axis="x", visible=False)
    ax.legend(loc="upper left", fontsize=7.5)
    fig.tight_layout()
    save(fig, "exp-donovan")


def _ratio(a, b):
    v = a[0] / b[0]
    return v, v * math.sqrt((a[1] / a[0]) ** 2 + (b[1] / b[0]) ** 2)


@figure("exp-kilic", ["recorded-results.json"])
def exp_kilic():
    d = data("recorded-results.json")["kilic"]
    mn = d["residence_measured_s"]["mn"]
    fig, axes = plt.subplots(1, 2, figsize=(7.2, 3.0))
    ax = axes[0]
    xs = [0, 25, 50, 75, 100]
    meas = [_ratio(d["residence_measured_s"][str(x)], mn) for x in xs]
    measured(ax, xs, [m[0] for m in meas], [m[1] for m in meas], label="Kilic et al. 2015")
    model(ax, xs, [d["residence_model_ratio"][str(x)] for x in xs], label="Model")
    ax.axhline(1, color=REF, ls="--", lw=1.2, label="single nucleosome")
    ax.set_xlabel("H3K9me3 on the array (%)")
    ax.set_ylabel("residence / mononucleosome residence")
    ax.set_title("Residence versus mark density", loc="left")
    ax.legend(loc="upper left", fontsize=7.5)
    ax = axes[1]
    concs = [50, 250, 1000]
    occ = [c / (c + d["kd_nM"]) for c in concs]
    meas = [_ratio(d["competition_measured_s"][str(c)], mn) for c in concs]
    measured(ax, occ, [m[0] for m in meas], [m[1] for m in meas], label="Kilic et al. 2015")
    model(ax, occ, [d["competition_model_ratio"][str(c)] for c in concs], label="Model")
    ax.axhline(1, color=REF, ls="--", lw=1.2, label="no-rebinding floor")
    for c, o in zip(concs, occ):
        ax.annotate(f"{c} nM", (o, 0.55), ha="center", fontsize=7, color=MUTED)
    ax.set_ylim(0.4, 2.6)
    ax.set_xlabel("competitor occupancy of sites")
    ax.set_title("Residence with competing HP1α (100% array)", loc="left")
    ax.legend(loc="upper right", fontsize=7.5)
    fig.tight_layout()
    save(fig, "exp-kilic")


@figure("cal-noise", ["recorded-results.json"])
def cal_noise():
    rows = data("recorded-results.json")["calibration_noise"]["rows"]
    s = [r["s"] for r in rows]
    fig, axes = plt.subplots(1, 2, figsize=(7.0, 2.6))
    model(axes[0], s, [r["time_s"] for r in rows], label="fitted clock")
    axes[0].plot(s, [r["time_s"] for r in rows], "-", color=BLUE, lw=1.2, alpha=0.6)
    axes[0].set_ylabel("seconds per reduced time unit")
    axes[0].set_title("Fitted clock", loc="left")
    model(axes[1], s, [r["log_rmse"] for r in rows], color=ORANGE, label="fit error")
    axes[1].plot(s, [r["log_rmse"] for r in rows], "-", color=ORANGE, lw=1.2, alpha=0.6)
    axes[1].set_ylabel("log RMSE of MSD fit")
    axes[1].set_title("Fit quality (best at s = 1)", loc="left")
    for ax in axes:
        ax.set_xlabel("noise offset scale s")
        ax.set_xticks(s)
    fig.tight_layout()
    save(fig, "cal-noise")


# ---------------------------------------------------------------- driver

STYLE_SOURCE = "".join(inspect.getsource(f) for f in (draw_report_figure, report_panels, checks_panel, box, arrow,
                                                       model, reference, measured, relaxation, save))


def fingerprint(name):
    fn, inputs = REGISTRY[name]
    h = hashlib.sha256(inspect.getsource(fn).encode() + STYLE_SOURCE.encode())
    for item in inputs:
        path = DATA / item
        h.update(path.read_bytes() if path.exists() else b"missing")
    return h.hexdigest()


def main(argv):
    force = "--force" in argv
    names = [a for a in argv if not a.startswith("--")] or list(REGISTRY)
    cache = json.loads(CACHE.read_text()) if CACHE.exists() else {}
    drawn, skipped, failed = [], [], []
    for name in names:
        _, inputs = REGISTRY[name]
        if any(not (DATA / i).exists() for i in inputs):
            failed.append(f"{name} (missing data: {', '.join(i for i in inputs if not (DATA / i).exists())})")
            continue
        fp = fingerprint(name)
        if not force and cache.get(name) == fp and (HERE / f"{name}.svg").exists():
            skipped.append(name)
            continue
        REGISTRY[name][0]()
        cache[name] = fp
        drawn.append(name)
    CACHE.write_text(json.dumps(cache, indent=1, sort_keys=True))
    print(f"drawn: {', '.join(drawn) or 'none'}")
    print(f"up to date: {len(skipped)}")
    if failed:
        print("not drawn: " + "; ".join(failed))


if __name__ == "__main__":
    main(sys.argv[1:])
