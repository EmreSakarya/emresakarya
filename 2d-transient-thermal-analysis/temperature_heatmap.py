"""
Temperature distribution map at the end of the heating phase.

Produces ``temperature_heatmap.png`` — a colour-coded map of the final nodal
temperatures (3x5 grid) once every node has exceeded the 1000 K target.
"""

import numpy as np
import matplotlib.pyplot as plt

# Final nodal temperatures at the end of the heating phase [K].
T_map = np.array([
    [1033.2, 1022.1, 1034.3],   # y = 0 (bottom nodes with clamps)
    [1009.4, 1008.2, 1009.1],   # y = 1
    [1002.3, 1001.8, 1002.1],   # y = 2 (centre interior nodes)
    [1008.5, 1007.1, 1008.3],   # y = 3
    [1029.0, 1021.0, 1034.3],   # y = 4 (top nodes with clamps)
])

fig, ax = plt.subplots(figsize=(6, 8))
im = ax.imshow(T_map, cmap="jet", origin="lower", aspect="auto", alpha=0.85)

cbar = fig.colorbar(im, ax=ax, fraction=0.046, pad=0.04)
cbar.set_label("Temperature [K]", rotation=270, labelpad=20, fontsize=12, fontweight="bold")

ax.set_xticks(np.arange(3))
ax.set_yticks(np.arange(5))
ax.set_xticklabels(["x = 0", "x = 1", "x = 2"], fontsize=11)
ax.set_yticklabels(["y = 0", "y = 1", "y = 2", "y = 3", "y = 4"], fontsize=11)
ax.set_xlabel("X Nodes (Horizontal)", fontsize=12, fontweight="bold", labelpad=10)
ax.set_ylabel("Y Nodes (Vertical)", fontsize=12, fontweight="bold", labelpad=10)
ax.set_title("2D Temperature Distribution Map\nEnd of Heating Phase",
             fontsize=13, fontweight="bold", pad=15)

for i in range(5):
    for j in range(3):
        text_color = "white" if (T_map[i, j] > 1025 or T_map[i, j] < 1005) else "black"
        ax.text(j, i, "%.1f K" % T_map[i, j], ha="center", va="center",
                color=text_color, fontweight="bold", fontsize=12)
        ax.text(j, i - 0.2, "Node(%d,%d)" % (j, i), ha="center", va="center",
                color=text_color, fontsize=9, alpha=0.9)

ax.set_xticks(np.arange(-0.5, 3, 1), minor=True)
ax.set_yticks(np.arange(-0.5, 5, 1), minor=True)
ax.grid(which="minor", color="black", linestyle="-", linewidth=2)
ax.tick_params(which="minor", bottom=False, left=False)

plt.tight_layout()
plt.savefig("temperature_heatmap.png", dpi=160, bbox_inches="tight")

if __name__ == "__main__":
    # plt.show()
    pass
