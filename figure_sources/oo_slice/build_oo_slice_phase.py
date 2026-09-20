"""Exact geometry of the OO diagonal and anti-diagonal coefficient slices.

Physical output: 15.6 cm x 9.5 cm; vector PDF and 300 dpi PNG preview.
No numerical connectedness test, sampled boundary, or point cloud is used.

For either iota(a)=(a,a) or iota(a)=(-a,a), 0<|a|<1,
  iota(a) in M_OO iff |a|>=1/sqrt(2) or (a real and |a|>=1/2).
The full-locus boundary intersection is exactly the critical circle plus
the real closed segments 1/2<=|a|<=1/sqrt(2).
"""
from pathlib import Path
import json
import math
import os

os.environ.setdefault("MPLCONFIGDIR", "/tmp/oo-slice-matplotlib")
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Circle, Rectangle


HERE = Path(__file__).resolve().parent
INK = "#203745"
GRAY = "#6D7D87"
LIGHT_GRAY = "#C5CDD2"
BOUNDARY = "#174E83"
CONNECTED = "#DCEAF8"
DUST = "#F6F2EB"

plt.rcParams.update({
    "font.family": "DejaVu Sans",
    "font.size": 9,
    "mathtext.fontset": "dejavusans",
    "axes.labelsize": 9,
    "xtick.labelsize": 8.5,
    "ytick.labelsize": 8.5,
    "text.color": INK,
    "axes.labelcolor": INK,
    "xtick.color": INK,
    "ytick.color": INK,
    "pdf.fonttype": 42,
    "ps.fonttype": 42,
    "savefig.facecolor": "white",
})

fig = plt.figure(figsize=(15.6/2.54, 9.5/2.54), facecolor="white")
ax = fig.add_axes([0.018, 0.135, 0.565, 0.765])
ax.set_aspect("equal")
ax.set_xlim(-1.22, 1.30)
ax.set_ylim(-1.13, 1.16)
ax.axis("off")

rcrit = 1/math.sqrt(2)
ax.add_patch(Circle((0, 0), 1, facecolor=CONNECTED, edgecolor="none", zorder=1))
ax.add_patch(Circle((0, 0), rcrit, facecolor=DUST, edgecolor="none", zorder=2))

# Thin coordinate axes; exact excluded and boundary sets are drawn above them.
arrow_kw = dict(arrowstyle="->", color=LIGHT_GRAY, lw=0.65, mutation_scale=7)
ax.annotate("", (1.20, 0), (-1.15, 0), arrowprops=arrow_kw, zorder=3)
ax.annotate("", (0, 1.13), (0, -1.10), arrowprops=arrow_kw, zorder=3)
ax.text(1.23, -0.035, r"$\mathrm{Re}\,a$", ha="left", va="center", fontsize=9)
ax.text(0.045, 1.10, r"$\mathrm{Im}\,a$", ha="left", va="bottom", fontsize=9)
for x, label in ((-1, r"$-1$"), (1, r"$1$")):
    ax.plot([x, x], [-0.024, 0.024], color=GRAY, lw=.65, zorder=4)
    ax.text(x + (.03 if x > 0 else -.03), -0.11, label,
            ha="left" if x > 0 else "right", va="top", fontsize=8.5)
ax.text(-.05, -.08, r"$0$", ha="right", va="top", fontsize=8.5)

ax.add_patch(Circle((0, 0), 1, facecolor="none", edgecolor=GRAY,
                    lw=.85, linestyle=(0, (3.2, 2.8)), zorder=4))
ax.add_patch(Circle((0, 0), rcrit, facecolor="none", edgecolor=BOUNDARY,
                    lw=1.45, zorder=5))
ax.plot([-rcrit, -.5], [0, 0], color=BOUNDARY, lw=2.8, solid_capstyle="round", zorder=6)
ax.plot([.5, rcrit], [0, 0], color=BOUNDARY, lw=2.8, solid_capstyle="round", zorder=6)
ax.plot([-.5, .5], [0, 0], linestyle="none", marker="o", ms=3.5,
        color=BOUNDARY, zorder=7)
ax.plot([0], [0], linestyle="none", marker="o", ms=4.5,
        markerfacecolor="white", markeredgecolor=GRAY, markeredgewidth=.8, zorder=7)
ax.text(.5, .105, r"$\frac{1}{2}$", ha="center", va="bottom", fontsize=9)
ax.text(-.5, .105, r"$-\frac{1}{2}$", ha="center", va="bottom", fontsize=9)

ax.text(0, .865, "connected", ha="center", va="center", fontsize=9)
ax.text(0, -.30, "disconnected", ha="center", va="center", fontsize=9)
ax.annotate(r"$|a|=1/\sqrt{2}$", xy=(.50, .50), xytext=(.81, .67),
            ha="left", va="center", fontsize=8.5, color=BOUNDARY,
            arrowprops=dict(arrowstyle="-", color=BOUNDARY, lw=.65,
                            connectionstyle="angle,angleA=0,angleB=45,rad=0"))

fig.text(.045, .966, "(a)  Coefficient plane", ha="left", va="top", fontsize=9.5, weight="bold")

# The angular graph retains isolated real-phase values and open lower endpoints.
ar = fig.add_axes([.687, .550, .280, .335])
ar.spines[["top", "right"]].set_visible(False)
ar.spines[["left", "bottom"]].set_color(LIGHT_GRAY)
ar.spines[["left", "bottom"]].set_linewidth(.7)
ar.set_xlim(-.18, 2*math.pi+.18)
ar.set_ylim(1.30, 2.12)
ar.set_xticks([0, math.pi, 2*math.pi], [r"$0$", r"$\pi$", r"$2\pi$"])
ar.set_yticks([math.sqrt(2), 2], [r"$\sqrt{2}$", r"$2$"])
ar.tick_params(length=3, width=.7, pad=3)
ar.set_xlabel(r"$\theta=\arg a$", labelpad=4)
ar.set_ylabel(r"$R(\theta)$", labelpad=5)
ar.plot([0, math.pi], [math.sqrt(2)]*2, color=BOUNDARY, lw=1.45)
ar.plot([math.pi, 2*math.pi], [math.sqrt(2)]*2, color=BOUNDARY, lw=1.45)
for theta in (0, math.pi, 2*math.pi):
    ar.plot(theta, math.sqrt(2), "o", ms=4.1, mfc="white", mec=BOUNDARY, mew=1.05, zorder=4)
    ar.plot(theta, 2, "o", ms=4.1, color=BOUNDARY, zorder=4)
ar.text(.52, 1.035, r"$\rho=1/|a|$", transform=ar.transAxes,
        ha="center", va="bottom", fontsize=8.5)
fig.text(.650, .966, "(b)  Radial ceiling", ha="left", va="top", fontsize=9.5, weight="bold")

# Compact vector legend. Text explains which ambient boundary is displayed.
legend_ax = fig.add_axes([.647, .130, .330, .285])
legend_ax.set_xlim(0, 1)
legend_ax.set_ylim(0, 1)
legend_ax.axis("off")
for y, color, label in ((.89, CONNECTED, "Connected"), (.70, DUST, "Disconnected")):
    legend_ax.add_patch(Rectangle((0, y-.04), .085, .08, facecolor=color, edgecolor="none"))
    legend_ax.text(.135, y, label, va="center", fontsize=8.5)
legend_ax.plot([0, .085], [.49, .49], color=BOUNDARY, lw=2)
legend_ax.text(.135, .49, "Boundary in the full 4D locus", va="center", fontsize=8.5)
legend_ax.plot([0, .085], [.28, .28], color=GRAY, lw=.85, linestyle=(0, (3.2, 2.8)))
legend_ax.text(.135, .28, r"Excluded: $|a|=1$", va="center", fontsize=8.5)
legend_ax.plot([.0425], [.08], "o", ms=4.5, mfc="white", mec=GRAY, mew=.8)
legend_ax.text(.135, .08, r"Excluded: $a=0$", va="center", fontsize=8.5)

fig.text(.050, .055, r"Both OO slices: $(a,b)=(a,a)$ and $(a,b)=(-a,a)$.",
         fontsize=9, ha="left", va="center")

HERE.mkdir(parents=True, exist_ok=True)
pdf = HERE / "OO_Exact_Slice_Phase.pdf"
png = HERE / "OO_Exact_Slice_Phase.png"
fig.savefig(pdf, metadata={"Title": "Exact opposite-opposite slice phase diagram",
                          "Author": "Bernat Espigule",
                          "Subject": "Exact connectedness, full-locus boundary and radial discontinuity"})
fig.savefig(png, dpi=300)
plt.close(fig)
(HERE / "OO_Exact_Slice_Phase.json").write_text(json.dumps({
    "physical_width_cm": 15.6,
    "physical_height_cm": 9.5,
    "font_sizes_pt": [8.5, 9.0, 9.5],
    "geometry": "Exact disks, closed real antenna segments and isolated radial values; no parameter classification by sampling.",
    "strict_coefficient_domain": "0 < |a| < 1",
    "connected_set": "|a| >= 1/sqrt(2), or a real and |a| >= 1/2",
    "full_locus_boundary_in_each_slice": "|a| = 1/sqrt(2), or a real and 1/2 <= |a| <= 1/sqrt(2)",
    "radial_ceiling": "R(theta)=2 for theta=0 mod pi; sqrt(2) otherwise",
    "png_dpi": 300,
}, indent=2)+"\n")
print(pdf)
print(png)
