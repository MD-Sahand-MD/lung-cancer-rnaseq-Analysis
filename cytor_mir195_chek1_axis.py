"""Proposed CYTOR–miR-195-5p–CHEK1 regulatory axis in LUAD (three-node network figure)."""
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch

plt.rcParams.update({
    "font.family": "DejaVu Sans",
    "mathtext.fontset": "dejavusans",
    "pdf.fonttype": 42,      # keep text editable in Illustrator/Inkscape
    "svg.fonttype": "none",
})

# ---- palette (Okabe-Ito, colour-blind safe) -------------------------------
C_LNC, C_MIR, C_MRNA = "#E69F00", "#0072B2", "#009E73"
FILL = {C_LNC: "#FBEBC8", C_MIR: "#CFE4F2", C_MRNA: "#C9EBDF"}
C_BIND, C_INHIB, C_COEXP = "#333333", "#B5301F", "#6B6B6B"
TXT, SUB = "#1A1A1A", "#555555"

# ---- layout ---------------------------------------------------------------
W, H = 11.0, 7.8
fig = plt.figure(figsize=(W, H), dpi=300)
ax = fig.add_axes([0, 0, 1, 1])
ax.set_xlim(0, W); ax.set_ylim(0, H); ax.axis("off")

HW, HH = 1.2, 0.5   # node half-width / half-height
nodes = {
    "CYTOR":      dict(xy=np.array([2.0, 2.3]), col=C_LNC,  kind="lncRNA"),
    "miR-195-5p": dict(xy=np.array([5.5, 6.0]), col=C_MIR,  kind="miRNA"),
    "CHEK1":      dict(xy=np.array([9.0, 2.3]), col=C_MRNA, kind="mRNA"),
}

def boundary(a, b, gap=0.07):
    """Point on the border of node a's box along the line towards node b."""
    ca, cb = nodes[a]["xy"], nodes[b]["xy"]
    d = (cb - ca) / np.linalg.norm(cb - ca)
    tx = HW / abs(d[0]) if d[0] else np.inf
    ty = HH / abs(d[1]) if d[1] else np.inf
    return ca + d * (min(tx, ty) + gap), d

# ---- edges ----------------------------------------------------------------
# 1) CYTOR – miR-195-5p : binding (solid, round terminals)
p1, _ = boundary("CYTOR", "miR-195-5p")
p2, _ = boundary("miR-195-5p", "CYTOR")
ax.plot([p1[0], p2[0]], [p1[1], p2[1]], color=C_BIND, lw=2.4, zorder=1,
        solid_capstyle="round", marker="o", markersize=8, markevery=[0, 1])

# 2) miR-195-5p ⊣ CHEK1 : inhibition (solid, T-bar at target)
q1, _ = boundary("miR-195-5p", "CHEK1")
q2, u = boundary("CHEK1", "miR-195-5p", gap=0.16)
ax.plot([q1[0], q2[0]], [q1[1], q2[1]], color=C_INHIB, lw=2.4, zorder=1,
        solid_capstyle="butt")
perp = np.array([-u[1], u[0]])
bar = [q2 - perp * 0.22, q2 + perp * 0.22]
ax.plot([bar[0][0], bar[1][0]], [bar[0][1], bar[1][1]], color=C_INHIB,
        lw=3.4, solid_capstyle="butt", zorder=1)

# 3) CYTOR – CHEK1 : dashed, undirected
r1, _ = boundary("CYTOR", "CHEK1")
r2, _ = boundary("CHEK1", "CYTOR")
ax.plot([r1[0], r2[0]], [r1[1], r2[1]], color=C_COEXP, lw=2.2, zorder=1,
        linestyle=(0, (5, 3.5)))

# ---- nodes ----------------------------------------------------------------
for name, n in nodes.items():
    x, y = n["xy"]
    ax.add_patch(FancyBboxPatch(
        (x - HW, y - HH), 2 * HW, 2 * HH,
        boxstyle="round,pad=0.0,rounding_size=0.22",
        fc=FILL[n["col"]], ec=n["col"], lw=2.4, zorder=3))
    ax.text(x, y + 0.11, name, ha="center", va="center", fontsize=16,
            fontweight="bold", color=TXT, zorder=4)
    ax.text(x, y - 0.24, n["kind"], ha="center", va="center", fontsize=9.5,
            color=SUB, zorder=4)

# ---- edge annotations -----------------------------------------------------
def note(x, y, lines, ha, ec):
    ax.text(x, y, "\n".join(lines), ha=ha, va="center", fontsize=10,
            color=TXT, linespacing=1.45, zorder=5,
            bbox=dict(boxstyle="round,pad=0.45", fc="white", ec=ec, lw=1.2))

note(3.05, 4.15, [
    "Published interaction evidence",
    "(DIANA-LncBase; Zhang and Li)",
    "LUAD expression correlation:",
    r"$r$ = 0.005, $P$ = 0.908 (ns)",
], ha="right", ec=C_BIND)

note(7.95, 4.15, [
    "Published direct targeting evidence",
    "LUAD negative correlation:",
    r"$r$ = −0.310, $P$ = 6.93 × 10$^{-13}$",
], ha="left", ec=C_INHIB)

ax.text(5.5, 2.3, "positive co-expression", ha="center", va="center",
        fontsize=10.5, style="italic", color=TXT, zorder=5,
        bbox=dict(boxstyle="round,pad=0.3", fc="white", ec="none"))
note(5.5, 1.35, [
    r"LUAD Spearman $\rho$ = 0.35, $P$ = 3.9 × 10$^{-15}$",
    "Indirect regulation remains hypothetical",
], ha="center", ec=C_COEXP)

# ---- title ----------------------------------------------------------------
ax.text(W / 2, 7.4, "Proposed CYTOR–miR-195-5p–CHEK1 regulatory axis in LUAD",
        ha="center", va="center", fontsize=16.5, fontweight="bold", color=TXT)

# ---- legend (edge types), drawn by hand so the terminals match the edges ---
ly = 0.68
def legend_item(x0, label, color, style):
    x1 = x0 + 0.75
    if style == "bind":
        ax.plot([x0, x1], [ly, ly], color=color, lw=2.4, marker="o", markersize=7,
                solid_capstyle="round")
    elif style == "inhib":
        ax.plot([x0, x1], [ly, ly], color=color, lw=2.4, solid_capstyle="butt")
        ax.plot([x1, x1], [ly - 0.12, ly + 0.12], color=color, lw=3.4,
                solid_capstyle="butt")
    else:
        ax.plot([x0, x1], [ly, ly], color=color, lw=2.2, linestyle=(0, (5, 3.5)))
    ax.text(x1 + (0.02 if style == "dash" else 0.18), ly, label, ha="left", va="center", fontsize=9.5, color=TXT)

legend_item(0.95, "Binding (published interaction)", C_BIND, "bind")
legend_item(4.20, "Inhibition (published targeting)", C_INHIB, "inhib")
legend_item(7.45, "Co-expression only (undirected)", C_COEXP, "dash")
ax.text(W / 2, 0.22, "ns, not significant.  Solid edges: published molecular "
        "interactions; dashed edge: correlation only.",
        ha="center", va="center", fontsize=8.5, color=SUB)

# ---- save: ask the user where to put the file -------------------------------
import tkinter as tk
from tkinter import filedialog

root = tk.Tk()
root.withdraw()                       # hide the empty main window
root.attributes("-topmost", True)     # keep the dialog in front
save_path = filedialog.asksaveasfilename(
    parent=root,
    title="Save figure as...",
    initialfile="cytor_mir195_chek1_axis",
    defaultextension=".png",
    filetypes=[("PNG image", "*.png"), ("PDF (vector)", "*.pdf"),
               ("SVG (vector)", "*.svg"), ("TIFF image", "*.tiff")],
)
root.destroy()

if save_path:                         # empty string means the user cancelled
    fig.savefig(save_path, dpi=300, facecolor="white")   # format from extension
    print(f"Saved to {save_path}")
else:
    print("Save cancelled - figure not written.")
