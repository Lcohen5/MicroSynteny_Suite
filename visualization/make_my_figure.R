# ══════════════════════════════════════════════════════════════════════════════
# MAKE MY FIGURE
#
# Run this file ONE LINE AT A TIME: click on a line, then press
#   Ctrl + Enter   (Windows)      Cmd + Return   (Mac)
# After each line, read what appears in the Console pane below.
#
# Lines that start with # are notes for people. R skips them.
# ══════════════════════════════════════════════════════════════════════════════

# 1. Load the tool. This draws nothing yet; it teaches R a new command,
#    plot_synteny(). You only need to do this once each time you open RStudio.
source("comparative_synteny_plot.R")

# 2. Draw every species in the sample table. Tables live in the output folder,
#    beside this one; the tool looks there, and saves the figures there too.
#    When it finishes, the Console names the three files it saved:
#    a PDF, an SVG and a PNG of the same figure.
plot_synteny("ImpL2_Sample.tsv")

# 3. Draw only one subgenus.
plot_synteny("ImpL2_Sample.tsv", level = "subgenus", name = "Sophophora")

# 4. Draw only one species group.
#    Try it: change "obscura" to "repleta" and run the line again.
plot_synteny("ImpL2_Sample.tsv", level = "group", name = "obscura")

# 5. Now your own table. The merge saved it in the output folder. Change the
#    file name inside the quotes to match it exactly, and run the line.
plot_synteny("ImpL2_merged.tsv")
