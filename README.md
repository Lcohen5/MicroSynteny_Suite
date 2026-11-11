# MicroSynteny Suite

Tools for recording the genes around a target gene in several *Drosophila*
species and drawing the result as a comparative synteny figure.

Created by Logan Cohen for the Genomics Education Partnership.

Contact: logancohen@thegep.org

Repository: https://github.com/Lcohen5/MicroSynteny_Suite

## Getting it

On the repository page, click the green **Code** button, then **Download ZIP**,
and unzip it. The folder is called `MicroSynteny_Suite-main`. In Posit Cloud,
choose **New Project from Git Repository** and paste the repository address.

## What is here

| Folder | What it does | Runs in |
|---|---|---|
| `table_builder/` | Interactive program that builds a table (TSV) of a target gene and its neighbours, one species at a time | Python, in JupyterLab (Anaconda Navigator) |
| `merge/` | Combines the tables for one target gene into one, and writes a report of what it kept and what to check | Python, in JupyterLab |
| `visualization/` | Draws the comparative figure from a merged table, as a PDF, an SVG and a PNG | R |
| `output/` | Where every table and figure is saved. Ships with one finished table, `ImpL2_Sample.tsv`, to try the visualization on | |

## Order of work

1. Build a table for your species with the Table Builder. It is saved in `output/`.
2. Merge the tables for one target gene with `merge/merge_microsynteny.py`.
   It reads the tables in `output/` and saves the merged table there.
3. Draw the figure with `visualization/comparative_synteny_plot.R`. It finds
   the merged table in `output/` and saves the figure there.

Keep the four folders together. The programs find `output/` by looking beside
their own folder.

Steps 1 and 2 run from one JupyterLab notebook kept in the suite's main
folder, one cell each:

    %run table_builder/MicroSynteny_Table_Builder.py
    %run merge/merge_microsynteny.py

Step 3 runs in RStudio: open `visualization/visualization.Rproj`, then
`make_my_figure.R`.

Fuller instructions for the Table Builder and the merge are at the top of each
file. Instructions for plotting are in `visualization/README.md`.

Name table files `gene_lastname.tsv`, for example `ImpL2_Smith.tsv`.

## Requirements

- Python 3 (Anaconda is the tested route)
- R 4.0 or newer, with the `tidyverse`, `ggnewscale` and `svglite` packages:
  `install.packages(c("tidyverse", "ggnewscale", "svglite"))`

## Acknowledgements

With thanks to the GEP faculty and students whose feedback made this possible:

- the students and Summer 2026 Reconcilers of the Reed Lab at the University of Alabama
- the GEP Advanced Research Lab Faculty Implementation Group
- participants in the Summer 2026 Faculty Reconciliation Trainings

## How this was made

MicroSynteny Suite was a hybrid effort. It was started by hand and refined with
AI assistance (Claude, from Anthropic). The author is responsible for the code
and its documentation.

## License

MIT. See `LICENSE`.
