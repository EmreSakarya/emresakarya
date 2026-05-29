"""
Methodology figure: the 2D nodal network and its control volumes (3x5 grid).

Produces ``node_diagram_final.png`` — a schematic of the plate geometry showing
the node positions, control-volume boundaries, corner heat sources (q'''), and
grid spacing.
"""

import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import matplotlib.lines as mlines

# Geometry in centimetres (schematic only).
nx = 3
ny = 5
Lx = 6.0
Ly = 16.0

x_nodes = [0.0, 3.0, 6.0]
y_nodes = [0.0, 4.0, 8.0, 12.0, 16.0]

x_cv_lines = [0.0, 1.5, 4.5, 6.0]
y_cv_lines = [0.0, 2.0, 6.0, 10.0, 14.0, 16.0]

CV_FILL = "#FDEBD0"
CV_EDGE = "#E67E22"
C_RED = "#C0392B"
C_BLUE = "#1A5276"
C_GREEN = "#1E8449"
C_TEXT = "#1C2833"


def node_color(i, j):
    corner = (i in (0, nx - 1)) and (j in (0, ny - 1))
    side = (i in (0, nx - 1)) or (j in (0, ny - 1))
    if corner:
        return C_RED, "Edge\n(Corner)"
    if side:
        return C_BLUE, "Side"
    return C_GREEN, "Interior"


q_vals = {(0, 0): 32, (2, 0): 35, (0, 4): 28, (2, 4): 37}

fig, ax = plt.subplots(figsize=(7.5, 12))
ax.set_xlim(-2.8, 9.0)
ax.set_ylim(-2.5, 20.2)
ax.set_aspect("equal")
ax.axis("off")
fig.patch.set_facecolor("white")

# title
ax.text(3.0, 19.6, "2D Nodal Network – Control Volumes  (Grid: 3×5)",
        ha="center", va="top", fontsize=11, fontweight="bold", color=C_TEXT)
ax.text(3.0, 19.05, "Lx = 6 cm,   Ly = 16 cm,   Lz = 1 cm",
        ha="center", va="top", fontsize=9.5, color=C_TEXT)

# plate + control volumes
plate = mpatches.Rectangle((0, 0), Lx, Ly, facecolor=CV_FILL,
                           edgecolor=CV_EDGE, linewidth=1.8, zorder=1)
ax.add_patch(plate)
for xb in x_cv_lines:
    ax.plot([xb, xb], [0, Ly], color=CV_EDGE, lw=1.0, zorder=2)
for yb in y_cv_lines:
    ax.plot([0, Lx], [yb, yb], color=CV_EDGE, lw=1.0, zorder=2)

# nodes
for i, xi in enumerate(x_nodes):
    for j, yj in enumerate(y_nodes):
        clr, lbl = node_color(i, j)
        ax.plot(xi, yj, "o", color=clr, markersize=10, zorder=5,
                markeredgecolor="white", markeredgewidth=0.6)
        if i == nx - 1:
            tx = xi + 0.45
            ha = "left"
        else:
            tx = xi - 0.45
            ha = "right"
        ax.text(tx, yj + 0.30, "(%d,%d)\n%s" % (i, j, lbl),
                fontsize=7, ha=ha, va="bottom", color=C_TEXT, linespacing=1.3)

# heat sources
q_pos = {
    (0, 4): (-1.0, 0.35, "center"),
    (2, 4): (1.0, 0.35, "center"),
    (0, 0): (-1.0, -1.15, "center"),
    (2, 0): (1.0, -1.15, "center"),
}
for (i, j), q in q_vals.items():
    xi = x_nodes[i]
    yj = y_nodes[j]
    ox, oy, ha = q_pos[(i, j)]
    ax.text(xi + ox, yj + oy, "q'''=%d kW/m³" % q,
            fontsize=8.5, ha=ha, color=C_RED, fontweight="bold",
            bbox=dict(boxstyle="round,pad=0.28", facecolor="#FDEDEC",
                      edgecolor=C_RED, linewidth=1.0))


def darrow(ax, x1, y1, x2, y2, label, lcolor="black",
           lx_off=0, ly_off=0.28, fs=8.5, rot=0):
    ax.annotate("", xy=(x2, y2), xytext=(x1, y1),
                arrowprops=dict(arrowstyle="<->", color=lcolor, lw=1.1, mutation_scale=12))
    mx = (x1 + x2) / 2 + lx_off
    my = (y1 + y2) / 2 + ly_off
    ax.text(mx, my, label, ha="center", va="bottom", fontsize=fs, color=lcolor, rotation=rot)


darrow(ax, 0, 17.5, 3, 17.5, "Δx = 3 cm")
darrow(ax, 3, 17.5, 6, 17.5, "Δx = 3 cm")
darrow(ax, -1.3, 4, -1.3, 8, "Δy = 4 cm", lx_off=-0.40, ly_off=0, fs=8, rot=90)
darrow(ax, -1.3, 0, -1.3, 4, "Δy = 4 cm", lx_off=-0.40, ly_off=0, fs=8, rot=90)
darrow(ax, 0.62, 0, 0.62, 2, "Δy/2", lcolor=CV_EDGE, lx_off=-0.42, ly_off=0, fs=8, rot=90)
darrow(ax, 0, -1.1, 1.5, -1.1, "Δx/2", lcolor=CV_EDGE, ly_off=0.22, fs=8)
darrow(ax, 0, -1.8, 6, -1.8, "Lx", ly_off=0.22, fs=9)
darrow(ax, 7.2, 0, 7.2, 16, "Ly", lx_off=0.42, ly_off=0, fs=9, rot=90)

info = ("CV Volumes:\n"
        "Corner   : (Δx/2)·(Δy/2)·Lz\n"
        "Side(s) : (Δx/2)·Δy·Lz\n"
        "Side(y) : Δx·(Δy/2)·Lz\n"
        "Interior: Δx·Δy·Lz")
ax.text(5.9, 9.5, info, fontsize=7, va="center", ha="right", linespacing=1.55,
        bbox=dict(boxstyle="round,pad=0.4", facecolor="#EBF5FB",
                  edgecolor="steelblue", linewidth=0.9))

leg_handles = [
    mlines.Line2D([], [], marker="o", linestyle="None", markerfacecolor=C_RED,
                  markeredgecolor="white", markersize=10, label="Edge/Corner  (4 nodes, q''' source)"),
    mlines.Line2D([], [], marker="o", linestyle="None", markerfacecolor=C_BLUE,
                  markeredgecolor="white", markersize=10, label="Side             (8 nodes)"),
    mlines.Line2D([], [], marker="o", linestyle="None", markerfacecolor=C_GREEN,
                  markeredgecolor="white", markersize=10, label="Interior         (3 nodes)"),
    mpatches.Patch(facecolor=CV_FILL, edgecolor=CV_EDGE, linewidth=1.2, label="Control Volume (CV)"),
]
ax.legend(handles=leg_handles, loc="lower center", bbox_to_anchor=(0.42, -0.11),
          fontsize=8.2, framealpha=0.95, ncol=1)

ax.text(3.0, -2.3, "X Dimension (cm)", ha="center", va="center", fontsize=10, color=C_TEXT)
ax.text(-2.6, 8.0, "Y Dimension (cm)", ha="center", va="center", fontsize=10, color=C_TEXT, rotation=90)

plt.tight_layout()
plt.savefig("node_diagram_final.png", dpi=160, bbox_inches="tight", facecolor="white")

if __name__ == "__main__":
    # plt.show()
    pass
