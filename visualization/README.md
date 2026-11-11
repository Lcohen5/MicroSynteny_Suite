# MicroSynteny Visualization

Part of **MicroSynteny Suite** (see the README in the main folder for the whole).

One R script that turns a merged MicroSynteny table into a comparative synteny figure suitable for a poster, paper, or lab meeting.

- `comparative_synteny_plot.R` — defines `plot_synteny()`, which renders a merged table as an arrow-and-ribbon figure
- `make_my_figure.R` — a short starter script to run line by line
- `visualization.Rproj` — the RStudio project file; open this to start
- `species_taxonomy.tsv` — species order and grouping (see `species_taxonomy_SOURCES.md`)

Merging tables is done before this step, in Python: see `merge/merge_microsynteny.py` in the main folder.

---

## About the developer

Created by Logan Cohen for the Genomics Education Partnership.

- Contact: logancohen@thegep.org
- ORCID: 0000-0002-2961-0334
- Repository: https://github.com/Lcohen5/MicroSynteny_Suite

## License

Released under the MIT License. See `LICENSE` in the main folder for full terms.

## Citation

If you use these scripts in published work, please cite the suite:

> Cohen, Logan. (2026). *MicroSynteny Suite* [Computer software]. https://github.com/Lcohen5/MicroSynteny_Suite

---

## Section 1 — How to use

### What the tool does

The visualization half of the suite takes a merged MicroSynteny table and produces a **slot-scale figure**: each species-defined position is drawn as a uniform column, so the plot reads at a glance as "who's next to what," across all species side by side.

One call produces one figure, saved three ways: a PDF, an SVG and a PNG of the same plot. Which species appear in it is up to you: the whole phylogeny, or any single subgenus, species group, or subgroup.

The visual grammar:

- Genes are drawn as pentagon-shaped arrows pointing in the direction of transcription
- Gene identity is encoded by fill color, so an ortholog is the same color across every species
- Nested genes tile a subrow directly above their container, at the same arrow height
- Ribbons connect matching orthologs between species, in four states:
    - **Teal** — same slot, same orientation. Fully conserved.
    - **Blue** — different slot, same orientation. The gene moved, but did not invert.
    - **Magenta** — different orientation. Drawn as a twist, so an inversion is visible from the shape as well as the color.
    - **Dashed, faint** — drift: the gene is absent from every species between the two it connects. Drawn behind everything else.
- Target-relative orientation is enforced: the target gene is always drawn pointing right in every species, with upstream on the left and downstream on the right, regardless of the target's actual genomic strand. Species that had to be mirrored are marked *(- strand)* above the target.

The two orientation colors are drawn from the Okabe-Ito colorblind-safe palette. Gene colors come from R's Spectral palette, chosen for perceptual band separation.

> **Note on orphan marks.** Earlier versions drew a grey X where a gene was present in one species and absent in its neighbour. That rendering is off by default — the absence of a ribbon already carries the information, and the X marks added visual noise. The classification is still computed if you switch `DRAW_ORPHAN_MARKS` to `TRUE` near the bottom of the plot script.

### What you'll need before you start

- R version 4.0 or newer
- The `tidyverse`, `ggnewscale` and `svglite` R packages
    - Install once with `install.packages(c("tidyverse", "ggnewscale", "svglite"))`
- A merged MicroSynteny table at schema v7, in the suite's `output` folder
- `species_taxonomy.tsv` in the same folder as the plot script

That last one is not optional. It supplies the top-to-bottom species order for every figure, and the subgenus/group/subgroup memberships used for subset plots.

If a species in your table is not listed in the taxonomy table, it is not drawn, and the script names it in a note. Add it to the table rather than editing the script.

### Plotting in RStudio

1. Open `visualization.Rproj` (double-click it). RStudio opens in this folder.
2. Open `make_my_figure.R` and run it one line at a time.

The lines in that file are all there is to it:

```r
source("comparative_synteny_plot.R")     # load the tool; draws nothing yet

plot_synteny("ImpL2_Sample.tsv")         # every species in the table
```

`source()` teaches R one new command, `plot_synteny()`. Each call to it stands alone: plotting a second table needs no clean-up after the first.

A table named without a path is looked for in the suite's `output` folder, which is also where the figures are saved. The target gene is read from the table's `target_gene` column; pass `target = "ImpL2"` to state it yourself.

### Plotting from a terminal

```
Rscript comparative_synteny_plot.R TSV [TARGET] [LEVEL] [NAME]

Rscript comparative_synteny_plot.R ImpL2_merged.tsv
Rscript comparative_synteny_plot.R ImpL2_merged.tsv ImpL2 group obscura
```

From the suite's main folder, the same commands start with `Rscript visualization/comparative_synteny_plot.R`.

### Plotting one subgenus, group, or subgroup

Add a level and a name:

```r
plot_synteny("chico_merged.tsv", level = "subgenus", name = "Sophophora")
plot_synteny("chico_merged.tsv", level = "group",    name = "repleta")
plot_synteny("chico_merged.tsv", level = "subgroup", name = "melanogaster")
```

Valid levels are `full`, `subgenus`, `group`, and `subgroup`. The name must match a value in the corresponding column of `species_taxonomy.tsv` — pass a name it doesn't recognise and it will tell you what the valid options are.

**A subset needs at least two species present in your data.** One species has nothing to compare against, so the script stops with a message naming the species it did find, rather than drawing a figure with no ribbons in it.

**Adjacency is recomputed within the subset.** Two species that are neighbours in the filtered plot get a ribbon even if unshown species separate them in the full phylogeny. This is deliberate — a subset figure is meant to stand alone — but it does mean a subset ribbon asserts less about phylogenetic distance than a full-phylogeny ribbon does. Every subset figure says so in its subtitle.

### What `plot_synteny()` takes

| Argument | Effect |
|---|---|
| `tsv` | The merged table: a file name in quotes |
| `target` | Target gene name; read from the table if left out |
| `level` / `name` | Which species to draw; `level` defaults to `"full"` |
| `taxonomy` | A taxonomy table other than the one beside the script |
| `output_prefix` | An output file name stem other than the default below |

The function also returns the plot, so `p <- plot_synteny(...)` keeps it for further changes with ggplot2.

### Output

Three files are saved in the suite's `output` folder. For a full-phylogeny run:

```
<Target>_synteny.pdf
<Target>_synteny.svg
<Target>_synteny.png
```

For a subset run:

```
<Target>_<level>_<name>_synteny.pdf     e.g. chico_subgroup_melanogaster_synteny.pdf
<Target>_<level>_<name>_synteny.svg
<Target>_<level>_<name>_synteny.png
```

The PDF and SVG are vector: sharp at any size, and editable in Inkscape or Illustrator. The SVG keeps its labels as text and opens in a web browser. The PNG is raster, convenient for slides, posters and review.

If the script has been moved away from the suite and there is no `output` folder beside its own, tables are looked for and figures are saved in the working directory instead.

Figure width is `W_INCHES` (default 16). Height scales with species count, with a floor that differs for full and subset plots — a two-species subgroup would otherwise float in a great deal of whitespace. Both, plus `DPI`, are constants near the top of `plot_synteny()`.

---

## Section 2 — How it works

### Ribbon geometry

Each ribbon is drawn as a filled polygon whose top and bottom edges are cubic Bezier curves. The curves leave and arrive with horizontal tangents, so the ribbon reads as flowing straight down out of one species and straight up into the next, curving smoothly in between rather than veering off at an angle.

For **preserved orientation** ribbons, the top edge connects upper-left to lower-left, and the bottom edge connects upper-right to lower-right. No crossing, no twist.

For **flipped orientation** ribbons, the destination x-coordinates swap. The two curves cross in the middle, creating a natural twist — the ribbon flips 180° between species, which is exactly the topology of a genomic inversion. Flipped ribbons carry two signals, color *and* twist, so a reader who cannot distinguish the two colors still sees the inversion clearly from the shape.

### The four ribbon states

Two questions are asked of every gene shared by two species: is it in the same slot, and is it pointing the same way. Those give three states — conserved, moved, flipped — and each is drawn in its own color.

The fourth is different in kind. A **drift** ribbon connects two species that share a gene which is absent from *every* species between them. It is drawn faint and dashed, behind the solid ribbons, because it represents a gene that recurs across the phylogeny without being stably positioned in it.

### Target-relative orientation

Comparing microsynteny across species is easier when every species is displayed with its target pointing the same way. If Sik3 is on the plus strand in D. mel and the minus strand in D. hydei, the raw genomic layout would put "upstream of Sik3" on opposite sides of the two species, and reading synteny across the phylogeny becomes an exercise in mental flipping.

The plot solves this by displaying every species with the target arrow always pointing right. For any species whose target is on the minus strand, every gene's *displayed* strand is the opposite of its raw strand, and the species is marked *(- strand)* above the target.

No x-coordinate mirroring is involved, and none is needed: `slot` is already biological in the v7 schema — negative always means upstream in the direction of transcription — so the horizontal layout is target-relative before the script starts.

The stored data is never modified. The `strand`, `slot`, `gene_start` and `gene_end` columns in the TSV remain the true genomic values.

Ribbon orientation comparison happens on the *displayed* strand, because that's the biologically meaningful comparison. Two orthologs that both point away from the target in their respective species share their relative orientation, even if their raw genomic strands differ.

### Nested genes

Nested genes are drawn as full-height arrows on a subrow directly above their container, tiled at uniform width across the container's slot span, in biological order relative to the container's displayed orientation. The container itself is drawn as a normal arrow on the main row.

This lets you see both genes with equal visual weight, and read the containment relationship from the vertical stacking rather than from a shrunken inset arrow. Drawing nested genes at their true genomic width produced slivers too small to label.

A target gene that sits inside a larger gene is handled by the same machinery: the container relationship is resolved into an ordinary parent-child link when the data loads, so nothing downstream needs a special case.

### Palette and readability

Gene identity uses R's Spectral palette, which gives distinct perceptual bands across the color wheel. With many unique genes (15+), the palette necessarily uses shades whose luminance varies substantially — some fills come out dark, some light.

The orientation ribbon colors are drawn from the Okabe-Ito palette, which is designed to remain distinguishable under the three common forms of color blindness.

### Species order and grouping

Species order and taxonomic grouping both come from `species_taxonomy.tsv`. Its `order` column determines the top-to-bottom arrangement; its `subgenus`, `group` and `subgroup` columns determine which species a subset plot includes.

Species present in the TSV but absent from the table are not drawn; the script names them in a note. Species in the table but absent from the TSV are simply skipped in the layout.

Adapting to a different phylogeny is a matter of editing that table. No code changes are required.

See `species_taxonomy_SOURCES.md` for where the classifications came from and which conventions were chosen where the literature disagrees.

---

## Section 3 — Developer notes

### File layout

The plot script is standalone. It needs the `tidyverse`, `ggnewscale` and `svglite` packages and `species_taxonomy.tsv` beside it, and nothing else from the suite.

### The plot script

`comparative_synteny_plot.R` defines one function, `plot_synteny()`, organized top-down: input resolution, taxonomy load and subset selection, output dimensions, palette, data load and validation, container-of-target resolution, arrow-polygon generation, ribbon generation, plot assembly, file save.

The `make_arrow()` helper converts a rectangle spec into a 5-vertex pentagon polygon suitable for `geom_polygon`. The `bezier_ribbon()` helper generates a curved ribbon polygon between two species rows. Both are pure functions and can be reused for other rendering experiments.

The script does not mutate its input, and its output files are overwritten cleanly on each run.

### Adapting to other clades

- **Species order and grouping.** Edit `species_taxonomy.tsv`. This is the only change required.
- **Arrow proportions.** `ARROW_HEIGHT`, `NESTED_ROW_OFFSET`, `SLOT_HALF_WIDTH`, and the `tip_frac` argument to `make_arrow()` control arrow shape and spacing.
- **Ribbon curvature.** The `ctrl_offset` calculation inside `bezier_ribbon()` (currently 55% of vertical distance) controls how much the Bezier curve bulges.
- **Colors.** `gene_palette` is set from `hcl.colors(..., "Spectral")`. `COLOR_PRESERVED`, `COLOR_MOVED`, `COLOR_FLIPPED` and `COLOR_ORPHAN` are constants.
- **Figure size.** `W_INCHES`, `H_INCHES`, `H_INCHES_MIN_SUBSET`, `H_PER_SPECIES`, `H_HEADER`, `DPI`.

### Contributing

I encourage you to report bugs, suggest improvements, and ask questions. I can be contacted at logancohen@thegep.org.

### Known limitations

- Very large gene counts (>30 unique genes) exceed the perceptual capacity of the Spectral palette; several fills become hard to tell apart. Subset plots substantially reduce the gene count and are the current workaround. A future release will consider a hierarchical fill approach.
- In a full-phylogeny plot, solid ribbons are only drawn between phylogenetically-adjacent species pairs. Non-adjacent conservation is present in the data and partly visible through drift ribbons, but is not otherwise drawn. This is a deliberate choice to avoid ribbon-soup at scale, and I recognize it hides some real signal.
- Subset plots recompute adjacency within the subset, so a ribbon in a subset figure and a ribbon in the full figure do not assert the same thing. The subtitle says so, but it is worth keeping in mind when the two are shown side by side.
- Gene labels are always drawn in black. On the darkest fills of the palette they are harder to read.
