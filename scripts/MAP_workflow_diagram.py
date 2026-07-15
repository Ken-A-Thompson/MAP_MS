"""
MAP_workflow_diagram.py
Reproduces the MAP pipeline flowchart (left: per-run processing,
right: per-marker × amplicon steps) as a publication-quality PNG.
"""

import os
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch

# ── Colours ───────────────────────────────────────────────────────────────────
C_GREY   = "#EFEFEF"    # neutral process box
C_BLUE   = "#C5DCF0"    # differs by read type
C_YELLOW = "#FFF0C0"    # conditional on marker / length
C_EDGE   = "#9AABB8"
C_ARROW  = "#555555"
C_SUB    = "#666666"    # sub-label colour

# ── Font sizes ────────────────────────────────────────────────────────────────
FS_MAIN  = 17           # box main label
FS_SUB   = 14           # box sub-label
FS_ANNOT = 14           # branch annotations ("PE = Yes" etc.)
FS_LEG   = 15           # legend
FS_HDR   = 14           # "for each marker × amplicon" header

# ── Layout constants ────────────────────────────────────────────────────────────
BW   = 3.9               # standard box width
BH   = 0.78              # standard box height (single-line label)
BHT  = 1.0               # tall box height (main + sub-label)
BW_B = 3.4               # branch (blue) box width

LCX = 4.0                # left column centre x
L1X = 1.6                # left branch x
L2X = 6.4                # right branch x

RCX = 11.9               # right column centre x

DIV_X = 8.4              # (unused — divider removed)
CHAN_X = 9.5             # vertical channel used to route the connector up


# ── Helpers ─────────────────────────────────────────────────────────────────────

def draw_box(ax, cx, cy, w, h, label, sublabel=None,
             color=C_GREY, fs=FS_MAIN):
    rect = FancyBboxPatch(
        (cx - w / 2, cy - h / 2), w, h,
        boxstyle="round,pad=0.10",
        facecolor=color, edgecolor=C_EDGE, linewidth=1.3, zorder=3
    )
    ax.add_patch(rect)
    if sublabel:
        ax.text(cx, cy + h * 0.17, label,
                ha="center", va="center",
                fontsize=fs, fontweight="bold", zorder=4)
        ax.text(cx, cy - h * 0.22, sublabel,
                ha="center", va="center",
                fontsize=FS_SUB, color=C_SUB, zorder=4)
    else:
        ax.text(cx, cy, label,
                ha="center", va="center",
                fontsize=fs, zorder=4)


def arrow(ax, x1, y1, x2, y2):
    ax.annotate(
        "", xy=(x2, y2), xytext=(x1, y1),
        arrowprops=dict(
            arrowstyle="-|>",
            color=C_ARROW,
            lw=1.6,
            mutation_scale=14,
        ),
        zorder=2,
    )


def hline(ax, x1, x2, y):
    ax.plot([x1, x2], [y, y], color=C_ARROW, lw=1.6, zorder=2)


def vline(ax, x, y1, y2):
    ax.plot([x, x], [y1, y2], color=C_ARROW, lw=1.6, zorder=2)


def branch_label(ax, x, y, text):
    ax.text(x, y, text, ha="center", va="center",
            fontsize=FS_ANNOT, color=C_SUB, zorder=5)


def draw_branch(ax, cx, top_y, l1x, l2x, bottom_y_l1, bottom_y_l2, label1, label2):
    """Draw a split from a single point (cx, top_y) down to two children,
    with an hline placed just above the child boxes and the PE labels
    placed in the clear space between the parent box and that hline."""
    hline_y = max(bottom_y_l1, bottom_y_l2) + 0.30
    vline(ax, cx, top_y, hline_y)
    hline(ax, l1x, l2x, hline_y)
    arrow(ax, l1x, hline_y, l1x, bottom_y_l1)
    arrow(ax, l2x, hline_y, l2x, bottom_y_l2)
    label_y = (top_y + hline_y) / 2
    branch_label(ax, l1x, label_y, label1)
    branch_label(ax, l2x, label_y, label2)


# ── Figure setup ──────────────────────────────────────────────────────────────
fig, ax = plt.subplots(figsize=(17.5, 16))
ax.set_xlim(-0.7, 16.1)
ax.set_ylim(0, 14.6)
ax.axis("off")

# ─────────────────────────────────────────────────────────────────────────────
# LEFT PANEL — per-run processing
# ─────────────────────────────────────────────────────────────────────────────

# Row y positions (top → bottom), evenly spaced with room for branch arrows
y0 = 13.3   # Collect run parameters
y1 = 11.6   # Merge PE / Merge raw fastq
y2 = 9.95   # Filter & trim reads
y3 = 8.3    # Demultiplex
y4 = 6.5    # Trim primers (×2)
y5 = 4.85   # Screen chimeras
y6 = 3.2    # Denoise / Cluster

# ── Collect run parameters ────────────────────────────────────────────────────
draw_box(ax, LCX, y0, BW, BH, "Collect run parameters")

# Branch point between collect and merge rows
draw_branch(ax, LCX, y0 - BH / 2, L1X, L2X,
            y1 + BHT / 2, y1 + BHT / 2, "PE = Yes", "PE = No")

# ── Merge row (blue) ──────────────────────────────────────────────────────────
draw_box(ax, L1X, y1, BW_B, BHT, "Merge PE reads",    color=C_BLUE)
draw_box(ax, L2X, y1, BW_B, BHT, "Merge raw fastq",   color=C_BLUE)

# Converge from merge to filter
jp1 = (y1 - BHT / 2 + y2 + BH / 2) / 2
vline(ax, L1X, y1 - BHT / 2, jp1)
vline(ax, L2X, y1 - BHT / 2, jp1)
hline(ax, L1X, L2X, jp1)
arrow(ax, LCX, jp1, LCX, y2 + BH / 2)

# ── Filter & trim reads ───────────────────────────────────────────────────────
draw_box(ax, LCX, y2, BW, BH, "Filter & trim reads")
arrow(ax, LCX, y2 - BH / 2, LCX, y3 + BH / 2)

# ── Demultiplex ───────────────────────────────────────────────────────────────
draw_box(ax, LCX, y3, BW, BH, "Demultiplex")

draw_branch(ax, LCX, y3 - BH / 2, L1X, L2X,
            y4 + BHT / 2, y4 + BHT / 2, "PE = Yes", "PE = No")

# ── Trim primers (blue, with sub-labels) ──────────────────────────────────────
draw_box(ax, L1X, y4, BW_B, BHT, "Trim primers",
         sublabel="5′ fwd · 3′ rev", color=C_BLUE)
draw_box(ax, L2X, y4, BW_B, BHT, "Trim primers",
         sublabel="both-end fwd & rev", color=C_BLUE)

# Converge to screen chimeras
jp2 = (y4 - BHT / 2 + y5 + BHT / 2) / 2
vline(ax, L1X, y4 - BHT / 2, jp2)
vline(ax, L2X, y4 - BHT / 2, jp2)
hline(ax, L1X, L2X, jp2)
arrow(ax, LCX, jp2, LCX, y5 + BHT / 2)

# ── Screen chimeras ───────────────────────────────────────────────────────────
draw_box(ax, LCX, y5, BW, BHT, "Screen chimeras",
         sublabel="per sample")

# Branch to denoise / cluster
bp3 = (y5 - BHT / 2 + y6 + BHT / 2) / 2
vline(ax, LCX, y5 - BHT / 2, bp3)
hline(ax, L1X, L2X, bp3)
arrow(ax, L1X, bp3, L1X, y6 + BHT / 2)
arrow(ax, L2X, bp3, L2X, y6 + BHT / 2)

# ── Denoise / Cluster (blue) ──────────────────────────────────────────────────
draw_box(ax, L1X, y6, BW_B, BHT, "Denoise → OTUs",
         sublabel="paired-end (unoise)", color=C_BLUE)
draw_box(ax, L2X, y6, BW_B, BHT, "Cluster → OTUs",
         sublabel="long-read (fast)", color=C_BLUE)

# ─────────────────────────────────────────────────────────────────────────────
# RIGHT PANEL — per-marker × amplicon
# ─────────────────────────────────────────────────────────────────────────────

r1 = 12.5   # Cluster run-wide OTUs
r2 = 10.9   # Chimera check (yellow)
r3 = 9.3    # Remove NUMTs & errors (yellow)
r4 = 7.7    # Assign taxonomy
r5 = 6.1    # Remove contaminants
r6 = 4.5    # BIN match (yellow)
r7 = 2.9    # Reporting

draw_box(ax, RCX, r1, BW, BH,  "Cluster run-wide OTUs")
arrow(ax, RCX, r1 - BH / 2,   RCX, r2 + BHT / 2)

draw_box(ax, RCX, r2, BW, BHT, "Chimera check",
         sublabel="non-PE · ≥ 200 bp", color=C_YELLOW)
arrow(ax, RCX, r2 - BHT / 2,  RCX, r3 + BHT / 2)

draw_box(ax, RCX, r3, BW, BHT, "Remove NUMTs & errors",
         sublabel="COI-5P · ≥ 500 bp", color=C_YELLOW)
arrow(ax, RCX, r3 - BHT / 2,  RCX, r4 + BH / 2)

draw_box(ax, RCX, r4, BW, BH,  "Assign taxonomy")
arrow(ax, RCX, r4 - BH / 2,   RCX, r5 + BH / 2)

draw_box(ax, RCX, r5, BW, BH,  "Remove contaminants")
arrow(ax, RCX, r5 - BH / 2,   RCX, r6 + BHT / 2)

draw_box(ax, RCX, r6, BW, BHT, "BIN match",
         sublabel="COI-5P · ≥ 300 bp", color=C_YELLOW)
arrow(ax, RCX, r6 - BHT / 2,  RCX, r7 + BH / 2)

draw_box(ax, RCX, r7, BW, BH,  "Reporting")

# ─────────────────────────────────────────────────────────────────────────────
# Connector: routes from bottom of the left column, down, across, and up into
# the top of the right column (illustrates the two panels are one continuous
# pipeline, split for layout only).
# ─────────────────────────────────────────────────────────────────────────────
y_merge   = y6 - BHT / 2 - 0.35   # merge point below the two bottom-left boxes
y_channel = 1.55                  # horizontal channel beneath everything
y_top     = r1 + BH / 2 + 0.55    # horizontal jog above the right column

vline(ax, L1X, y6 - BHT / 2, y_merge)
vline(ax, L2X, y6 - BHT / 2, y_merge)
hline(ax, L1X, L2X, y_merge)
vline(ax, LCX, y_merge, y_channel)
hline(ax, LCX, CHAN_X, y_channel)
vline(ax, CHAN_X, y_channel, y_top)
hline(ax, CHAN_X, RCX, y_top)
arrow(ax, RCX, y_top, RCX, r1 + BH / 2)

ax.text((CHAN_X + RCX) / 2, y_top + 0.22, "for each marker × amplicon",
        ha="center", va="bottom",
        fontsize=FS_HDR, color=C_SUB, style="italic")

# ── Legend ────────────────────────────────────────────────────────────────────
leg_y = 0.7
lx = 0.35
for color, label in [
    (C_BLUE,   "Differs by read type"),
    (C_YELLOW, "Conditional on marker / length"),
]:
    rect = FancyBboxPatch(
        (lx, leg_y - 0.22), 0.42, 0.42,
        boxstyle="round,pad=0.06",
        facecolor=color, edgecolor=C_EDGE, linewidth=1.0, zorder=3
    )
    ax.add_patch(rect)
    ax.text(lx + 0.55, leg_y, label,
            ha="left", va="center", fontsize=FS_LEG)
    lx += 3.8

# ── Save ──────────────────────────────────────────────────────────────────────
script_dir = os.path.dirname(os.path.abspath(__file__))
out = os.path.join(script_dir, "..", "figs_tables", "MAP_workflow_diagram.png")
plt.savefig(out, dpi=250, bbox_inches="tight", facecolor="white")
print(f"Saved: {out}")
