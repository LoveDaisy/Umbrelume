"""Shared matplotlib style constants: the repo's single display-style definition.

This module is a zero-import leaf: figure tooling outside the solver
(tools/render_mainline.py, validation/*) must be able to consume the shared
style without importing any solver code.
"""

# Keys shared by every figure-producing script in the repo (eclipse.report,
# tools/render_mainline, validation/*). Scripts with additional display
# preferences layer them on top of BASE_RC instead of redefining shared keys.
BASE_RC = {
    "font.family": "DejaVu Sans", "font.size": 10,
    "axes.spines.top": False, "axes.spines.right": False,
    "figure.facecolor": "white", "axes.facecolor": "white",
    "savefig.facecolor": "white",
}
HOUSE_STYLE = {**BASE_RC, "axes.labelcolor": "#30343a", "text.color": "#252a31",
               "axes.edgecolor": "#bdc2c6", "xtick.color": "#545d66",
               "ytick.color": "#545d66"}
