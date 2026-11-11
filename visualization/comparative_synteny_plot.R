# ══════════════════════════════════════════════════════════════════════════════
# COMPARATIVE SYNTENY PLOT (v16)
#
# Draws one comparative synteny figure (as a PDF, an SVG and a PNG) from a merged
# MicroSynteny table (schema v7).
#
# ── HOW TO USE IT IN RSTUDIO ─────────────────────────────────────────────────
#   1. Open visualization.Rproj (double-click it). RStudio opens in this folder.
#   2. Open make_my_figure.R and run it one line at a time.
#
#   The lines in that file are all there is to it:
#
#       source("comparative_synteny_plot.R")       # load the tool (once)
#       plot_synteny("ImpL2_merged.tsv")           # every species
#       plot_synteny("ImpL2_merged.tsv", level = "subgenus", name = "Sophophora")
#       plot_synteny("ImpL2_merged.tsv", level = "group",    name = "obscura")
#       plot_synteny("ImpL2_merged.tsv", level = "subgroup", name = "melanogaster")
#
# ── HOW TO USE IT FROM A TERMINAL ────────────────────────────────────────────
#       Rscript comparative_synteny_plot.R TSV [TARGET] [LEVEL] [NAME]
#       Rscript comparative_synteny_plot.R ImpL2_merged.tsv
#       Rscript comparative_synteny_plot.R ImpL2_merged.tsv ImpL2 group obscura
#
# ── WHAT plot_synteny() TAKES ────────────────────────────────────────────────
#   tsv            the merged table (a file name in quotes)
#   target         the target gene; read from the table if you leave it out
#   level, name    which species to draw: level is "full" (default),
#                  "subgenus", "group" or "subgroup"; name is one of the names
#                  in that column of species_taxonomy.tsv
#   taxonomy       a taxonomy table other than the one beside this script
#   output_prefix  a file name stem (with a path, if wanted) other than the
#                  default; a stem without a path saves in the working directory
#
#   The suite keeps tables and figures in one folder, output/, beside this
#   script's folder. A table named without a path is looked for there, and the
#   figure is saved there, as <target>_synteny.pdf, .svg and .png
#   (or <target>_<level>_<name>_synteny for a subset). The plot is also
#   returned, so   p <- plot_synteny(...)   keeps it for further changes.
#
# What v16 changed from v15:
#   - The script now defines a function, plot_synteny(), instead of reading
#     variables left in the R session. Each call stands alone, so plotting a
#     second file can no longer reuse the first one's target or file name.
#   - species_taxonomy.tsv is found beside this script wherever it is run from.
#   - A target gene that is not in the table stops with a message instead of
#     drawing an unoriented figure.
#   - Species missing from species_taxonomy.tsv are named instead of dropped
#     silently.
#   - The PNG always has a white background.
#   - An SVG is saved alongside the PDF and PNG (needs the svglite package).
#   The drawing itself is unchanged from v15.
#
# Conventions:
#   - Genes drawn as pentagon-shaped arrows, pointing in transcription direction
#     within the target-relative frame
#   - Container genes drawn as full-size arrows on the main row; their
#     nested contents tile a uniform-width subrow above
#   - Ribbons connect matching gene_names between species:
#       teal    = same slot AND same relative orientation (fully conserved)
#       blue    = different slot, same orientation (positional shift, no inversion)
#       magenta = different orientation (inversion; drawn as twisted ribbon)
#       dashed  = drift: gene absent from every intervening species (any state)
# ══════════════════════════════════════════════════════════════════════════════

# ── Packages ─────────────────────────────────────────────────────────────────
.synteny_needs <- c("tidyverse", "ggnewscale", "svglite")
.synteny_missing <- .synteny_needs[
  !vapply(.synteny_needs, requireNamespace, logical(1), quietly = TRUE)]
if (length(.synteny_missing) > 0) {
  stop(sprintf(paste0(
    "This script needs R package(s) that are not installed: %s\n",
    "Run this line once, then try again:\n",
    "  install.packages(c(\"tidyverse\", \"ggnewscale\", \"svglite\"))"),
    paste(.synteny_missing, collapse = ", ")), call. = FALSE)
}
suppressPackageStartupMessages({
  library(tidyverse)
  library(ggnewscale)
})

# ── Where this script lives ──────────────────────────────────────────────────
# Recorded when the script is loaded, so plot_synteny() can find
# species_taxonomy.tsv beside it whatever the working directory is.
.SYNTENY_SCRIPT_DIR <- local({
  from_rscript <- grep("^--file=", commandArgs(trailingOnly = FALSE),
                       value = TRUE)
  from_source <- NULL
  for (i in seq_len(sys.nframe())) {
    f <- sys.frame(i)$ofile
    if (!is.null(f)) from_source <- f
  }
  p <- if (!is.null(from_source)) {
    from_source
  } else if (length(from_rscript) > 0) {
    # Rscript writes spaces in the path as "~+~".
    gsub("~+~", " ", sub("^--file=", "", from_rscript[1]), fixed = TRUE)
  } else {
    NULL
  }
  if (is.null(p)) getwd() else dirname(normalizePath(p, mustWork = FALSE))
})

plot_synteny <- function(tsv, target = NULL, level = "full", name = NULL,
                         taxonomy = NULL, output_prefix = NULL) {

  if (missing(tsv) || length(tsv) == 0 || !is.character(tsv)) {
    stop("Give the table to plot, in quotes: plot_synteny(\"ImpL2_merged.tsv\")",
         call. = FALSE)
  }
  TSV_PATH              <- tsv
  TARGET_GENE_CANONICAL <- target
  SUBSET_LEVEL          <- level
  SUBSET_NAME           <- if (is.null(name)) NA_character_ else name

  # The suite keeps tables and figures in one shared folder, output/, beside
  # this script's folder. A table named without a path is looked for in the
  # working directory first and then there; figures are saved there. If the
  # script has been moved away from the suite, the working directory is used.
  .shared_dir <- file.path(dirname(.SYNTENY_SCRIPT_DIR), "output")
  OUTPUT_DIR <- if (dir.exists(.shared_dir)) {
    normalizePath(.shared_dir)
  } else {
    getwd()
  }
  .in_shared <- !file.exists(TSV_PATH) &
    file.exists(file.path(OUTPUT_DIR, TSV_PATH))
  TSV_PATH[.in_shared] <- file.path(OUTPUT_DIR, TSV_PATH[.in_shared])

  .not_found <- TSV_PATH[!file.exists(TSV_PATH)]
  if (length(.not_found) > 0) {
    .tables <- function(dir) {
      f <- list.files(dir, pattern = "\\.tsv$", ignore.case = TRUE)
      f <- setdiff(f, "species_taxonomy.tsv")
      if (length(f) == 0) "none" else paste(f, collapse = ", ")
    }
    stop(sprintf(paste0(
      "Cannot find the file: %s\n",
      "Tables in %s:\n  %s\n",
      "Check the spelling, including capital letters and the .tsv ending."),
      paste(.not_found, collapse = ", "), OUTPUT_DIR, .tables(OUTPUT_DIR)),
      call. = FALSE)
  }

  # If TSV_PATH is a single directory, discover .tsv files automatically.
  if (length(TSV_PATH) == 1 && dir.exists(TSV_PATH)) {
    TSV_PATH <- list.files(TSV_PATH, pattern = ".*\\.tsv$",
                           full.names = TRUE)
  }

  # ── Species taxonomy, display order, and subsetting ──────────────────────────
  # Species order AND the subgenus/group/subgroup memberships used for subset
  # plots both come from species_taxonomy.tsv. Adapting the suite to a different
  # clade is a matter of editing that table -- no code change is required.
  #
  # The table is looked for, in order:
  #   1. the `taxonomy` argument, if given
  #   2. next to this script
  #   3. the working directory
  #
  # Adjacency is recomputed WITHIN the subset: two species that are neighbours in
  # the filtered plot get a ribbon even if unshown species separate them in the
  # full phylogeny. That is deliberate -- a subset plot is meant to stand alone --
  # but it does mean a subset ribbon asserts less about phylogenetic distance
  # than a full-phylogeny ribbon does.
  if (!is.null(taxonomy)) {
    if (!file.exists(taxonomy)) {
      stop(sprintf("Cannot find the taxonomy table: %s", taxonomy), call. = FALSE)
    }
    TAXONOMY_PATH <- taxonomy
  } else {
    .candidates <- unique(c(file.path(.SYNTENY_SCRIPT_DIR, "species_taxonomy.tsv"),
                            file.path(getwd(), "species_taxonomy.tsv")))
    .found <- .candidates[file.exists(.candidates)]
    if (length(.found) == 0) {
      stop(sprintf(paste0(
        "Cannot find species_taxonomy.tsv. Looked in:\n  %s\n",
        "Keep it in the same folder as comparative_synteny_plot.R."),
        paste(.candidates, collapse = "\n  ")), call. = FALSE)
    }
    TAXONOMY_PATH <- .found[1]
  }

  SPECIES_TAXONOMY <- readr::read_tsv(TAXONOMY_PATH, show_col_types = FALSE,
                                      progress = FALSE)

  .required_cols <- c("species", "subgenus", "group", "subgroup", "order")
  .missing_cols <- setdiff(.required_cols, names(SPECIES_TAXONOMY))
  if (length(.missing_cols) > 0) {
    stop(sprintf("%s is missing required column(s): %s",
                 TAXONOMY_PATH, paste(.missing_cols, collapse = ", ")),
         call. = FALSE)
  }

  SPECIES_TAXONOMY <- SPECIES_TAXONOMY |> arrange(order)
  SPECIES_ORDER <- SPECIES_TAXONOMY$species

  if (anyDuplicated(SPECIES_ORDER) > 0) {
    stop(sprintf("%s lists these species more than once: %s", TAXONOMY_PATH,
                 paste(unique(SPECIES_ORDER[duplicated(SPECIES_ORDER)]),
                       collapse = ", ")), call. = FALSE)
  }


  # Resolve the requested subset.
  SUBSET_LEVEL <- tolower(as.character(SUBSET_LEVEL))

  .valid_levels <- c("full", "subgenus", "group", "subgroup")
  if (!SUBSET_LEVEL %in% .valid_levels) {
    stop(sprintf("level must be one of %s (got '%s')",
                 paste(sQuote(.valid_levels), collapse = ", "), SUBSET_LEVEL),
         call. = FALSE)
  }

  if (SUBSET_LEVEL == "full") {
    SUBSET_NAME <- NA_character_
  } else {
    .known <- unique(SPECIES_TAXONOMY[[SUBSET_LEVEL]])
    if (is.na(SUBSET_NAME) || SUBSET_NAME == "") {
      stop(sprintf("level = \"%s\" also needs a name, like name = \"%s\".\nAvailable: %s",
                   SUBSET_LEVEL, sort(.known)[1],
                   paste(sort(.known), collapse = ", ")), call. = FALSE)
    }
    if (!SUBSET_NAME %in% .known) {
      stop(sprintf("No %s named '%s' in %s.\nAvailable: %s",
                   SUBSET_LEVEL, SUBSET_NAME, basename(TAXONOMY_PATH),
                   paste(sort(.known), collapse = ", ")), call. = FALSE)
    }
    SPECIES_ORDER <- SPECIES_TAXONOMY$species[
      SPECIES_TAXONOMY[[SUBSET_LEVEL]] == SUBSET_NAME]
    message(sprintf("Subset: %s '%s' (%d species in the taxonomy)",
                    SUBSET_LEVEL, SUBSET_NAME, length(SPECIES_ORDER)))
  }

  # ── Output dimensions ─────────────────────────────────────────────────────────
  W_INCHES <- 16
  H_INCHES <- 8    # height floor for a full-phylogeny plot
  DPI      <- 300
  # Height is driven by species count. The floor above suits a 20+ species figure
  # but leaves a two-species subset floating in whitespace, so subset plots get a
  # smaller floor. The per-species rate is the same either way.
  H_INCHES_MIN_SUBSET <- 3.5
  H_PER_SPECIES       <- 1.2
  # Headroom for the title block, which is now three subtitle lines tall.
  H_HEADER            <- 1.6

  # ── Orientation colors (Okabe-Ito colorblind-safe) ────────────────────────────
  COLOR_PRESERVED <- "#009E73"   # teal — same slot, same orientation
  COLOR_MOVED     <- "#0072B2"   # blue — different slot, same orientation
  COLOR_FLIPPED   <- "#CC79A7"   # magenta — different orientation (with or without slot change)
  COLOR_ORPHAN    <- "#7F7F7F"   # grey (for X-mark synteny breaks)

  # ── Target-relative orientation ───────────────────────────────────────────────
  # When TRUE, every species is displayed so the target gene points RIGHT.
  # For species where the target is on the '-' strand, this means:
  #   - Every gene's arrow direction is flipped in the display (strand-of-display)
  #   - Ribbon "preserved vs flipped" reflects strand relative to target rather
  #     than raw genomic strand
  # The underlying data is unchanged; this is purely a rendering choice.
  # When FALSE, arrows show raw strand as recorded in the TSV.
  ORIENT_TO_TARGET <- TRUE

  # ══════════════════════════════════════════════════════════════════════════════
  # LOAD
  # ══════════════════════════════════════════════════════════════════════════════
  message("Reading TSV(s)...")

  # Read every file and rbind. Target_gene normalization happens AFTER
  # reading so we can infer TARGET_GENE_CANONICAL from the file itself
  # when it wasn't passed on the command line.
  d_list <- list()
  for (path in TSV_PATH) {
    if (!file.exists(path)) {
      message(sprintf("  Skipping (not found): %s", path))
      next
    }
    message(sprintf("  Reading: %s", path))
    d_one <- read_tsv(path, show_col_types = FALSE)
    d_list[[path]] <- d_one
  }
  if (length(d_list) == 0) stop("No TSVs could be read.")
  d <- bind_rows(d_list, .id = "source_file")

  # ── Schema gate ──────────────────────────────────────────────────────────────
  # This script is written against schema v7. Handing it an older file used to
  # fail deep inside a dplyr verb ("object 'target_nesting_gene' not found"),
  # which tells a student nothing about what is actually wrong or how to fix it.
  # Check up front and say exactly what to run.
  REQUIRED_SCHEMA <- "v7"
  V7_ONLY_COLUMNS <- c("target_nesting_gene", "species_number")

  .schema_seen <- if ("schema_version" %in% names(d)) {
    unique(as.character(d$schema_version[!is.na(d$schema_version) &
                                         d$schema_version != ""]))
  } else character(0)
  .missing_v7 <- setdiff(V7_ONLY_COLUMNS, names(d))

  if (length(.missing_v7) > 0 || !all(.schema_seen == REQUIRED_SCHEMA)) {
    .found <- if (length(.schema_seen) == 0) "none recorded"
              else paste(.schema_seen, collapse = ", ")
    .files <- paste(basename(names(d_list)), collapse = ", ")
    stop(sprintf(paste0(
      "This file is not on schema %s, so it cannot be plotted yet.\n",
      "  file(s):          %s\n",
      "  schema_version:   %s\n",
      "  missing column(s): %s\n\n",
      "It was made with an older Table Builder and needs converting to\n",
      "schema %s before it can be merged or plotted."),
      REQUIRED_SCHEMA, .files, .found,
      if (length(.missing_v7) == 0) "none" else paste(.missing_v7, collapse = ", "),
      REQUIRED_SCHEMA), call. = FALSE)
  }

  # Infer canonical target_gene from the file if none was supplied on the
  # command line. The most-common target_gene value wins (which handles
  # both single-file and merged-with-normalized-labels cases).
  if (is.null(TARGET_GENE_CANONICAL)) {
    target_counts <- d |>
      filter(!is.na(target_gene), target_gene != "") |>
      count(target_gene, sort = TRUE)
    if (nrow(target_counts) == 0) {
      stop("No target_gene values found in any file. Pass the target as the second CLI arg.")
    }
    TARGET_GENE_CANONICAL <- target_counts$target_gene[1]
    message(sprintf("  Target gene: %s (read from the table)", TARGET_GENE_CANONICAL))
    if (nrow(target_counts) > 1) {
      other_targets <- paste(target_counts$target_gene[-1], collapse = ", ")
      message(sprintf(
        "  Note: file also contains target_gene values [%s]. Rows with those values will still be normalized to '%s'. Pass an explicit target on the CLI to override.",
        other_targets, TARGET_GENE_CANONICAL))
    }
  }

  # Now normalize every row's target_gene to the canonical.
  d$target_gene <- TARGET_GENE_CANONICAL

  # Multi-file merge: each faculty file typically contains D. mel + one other
  # species. To avoid duplicating D. mel across files (each with slightly
  # different gene name variants like "tra_Schwartz"), for each species we
  # keep rows from a single source file — the one contributing the most rows
  # for that species. Prevents the color palette from treating "tra" and
  # "tra_Schwartz" as different genes.
  species_source_counts <- d |>
    group_by(species, source_file) |>
    summarise(n = n(), .groups = "drop") |>
    group_by(species) |>
    slice_max(n, n = 1, with_ties = FALSE) |>
    ungroup() |>
    select(species, canonical_source = source_file)

  d <- d |>
    inner_join(species_source_counts, by = "species") |>
    filter(source_file == canonical_source) |>
    select(-source_file, -canonical_source)

  # Species the taxonomy table does not list cannot be placed in the figure.
  # Say so by name: a species that silently vanishes looks like missing data.
  .unplaced <- setdiff(unique(d$species), SPECIES_TAXONOMY$species)
  if (length(.unplaced) > 0) {
    message(sprintf(
      "  Note: not listed in %s, so not drawn: %s",
      basename(TAXONOMY_PATH), paste(sort(.unplaced), collapse = ", ")))
  }

  d <- d |> filter(species %in% SPECIES_ORDER)

  # Explicitly coerce numeric columns. If read_tsv detected a column as
  # character (because a row was empty, or had "NA", or a stray non-numeric
  # character), later math like `slot - SLOT_HALF_WIDTH` will fail with
  # "non-numeric argument to binary operator" and R will fall through to the
  # methods::slot() function producing very confusing error messages.
  d <- d |>
    mutate(
      slot        = suppressWarnings(as.numeric(slot)),
      gene_start  = suppressWarnings(as.numeric(gene_start)),
      gene_end    = suppressWarnings(as.numeric(gene_end)),
      strand      = as.character(strand),
      gene_name   = as.character(gene_name),
      locid       = as.character(locid),
      feature_uid = as.character(feature_uid),
      parent_uid  = as.character(parent_uid)
    )

  # Report any rows dropped due to malformed slot/coords
  n_before <- nrow(d)
  d <- d |> filter(!is.na(slot), !is.na(gene_start), !is.na(gene_end))
  n_dropped <- n_before - nrow(d)
  if (n_dropped > 0) {
    message(sprintf("  Warning: dropped %d rows with non-numeric slot/coords", n_dropped))
  }

  species_present <- intersect(SPECIES_ORDER, unique(d$species))
  message(sprintf("  Found %d rows across %d species", nrow(d), length(species_present)))

  # A comparative figure needs something to compare. One species yields a
  # neighbourhood diagram with no ribbons, which is not what this script is for;
  # the batch runner pre-filters these, so reaching here means a direct call.
  if (length(species_present) < 2) {
    stop(sprintf(
      "Need at least 2 species to draw a comparative figure; %s has %d (%s).",
      if (SUBSET_LEVEL == "full") "this dataset"
        else sprintf("%s '%s'", SUBSET_LEVEL, SUBSET_NAME),
      length(species_present),
      if (length(species_present) == 0) "none" else paste(species_present, collapse = ", ")),
      call. = FALSE)
  }

  d$row_y <- setNames(seq_along(species_present),
                      rev(species_present))[as.character(d$species)]

  # Validate parent_uid: any value that isn't a real feature_uid in the file
  # is treated as top-level (empty) with a warning. This catches data-entry
  # errors where a gene name was typed into the parent_uid column instead of
  # the UID hash — those rows would otherwise silently fail to render.
  valid_uids <- unique(d$feature_uid)
  bad_parents <- d |>
    filter(!is.na(parent_uid) & parent_uid != "" & !(parent_uid %in% valid_uids))
  if (nrow(bad_parents) > 0) {
    message(sprintf(
      "  Warning: %d row(s) have parent_uid values that are not valid feature_uids.",
      nrow(bad_parents)))
    message("  These rows will be treated as top-level. Rows affected:")
    bad_summary <- bad_parents |>
      select(species, gene_name, parent_uid) |>
      distinct()
    for (i in seq_len(nrow(bad_summary))) {
      message(sprintf("    - %s / %s: parent_uid = %s",
                      bad_summary$species[i],
                      bad_summary$gene_name[i],
                      bad_summary$parent_uid[i]))
    }
    d <- d |>
      mutate(parent_uid = ifelse(parent_uid %in% valid_uids, parent_uid, ""))
  }

  # ══════════════════════════════════════════════════════════════════════════════
  # TARGET-RELATIVE ORIENTATION
  # For every species whose target sits on the '-' strand, flip the displayed
  # direction of every gene so the target points right. Concretely: those
  # species get a `display_strand` that is the opposite of their raw strand.
  #
  # No x-coordinate mirroring happens here, and none is needed: `slot` is already
  # biological in the v7 schema (negative = upstream in the direction of
  # transcription), so the horizontal layout is target-relative before we start.
  # Coordinate mirroring belonged to the genomic-scale plot, which was retired.
  #
  # This is display-only. The stored `strand`, `slot`, `gene_start` and `gene_end`
  # columns are never modified.
  # ══════════════════════════════════════════════════════════════════════════════
  # Find the target's raw strand per species (from the row at slot 0 whose
  # gene_name matches target_gene AND has no parent — i.e. the actual target,
  # not a same-slot container/nested)
  target_strand_by_species <- d |>
    filter(slot == 0, gene_name == TARGET_GENE_CANONICAL,
           is.na(parent_uid) | parent_uid == "") |>
    select(species, target_strand = strand) |>
    distinct()

  # The target has to be in the table under exactly this name, or no species
  # can be oriented and the figure would be wrong without looking wrong.
  if (!any(d$gene_name == TARGET_GENE_CANONICAL)) {
    .at_zero <- sort(unique(d$gene_name[d$slot == 0]))
    stop(sprintf(paste0(
      "No gene named '%s' in this table, so the figure cannot be oriented.\n",
      "Genes at the target position are: %s\n",
      "Check the spelling of the target, including capital letters."),
      TARGET_GENE_CANONICAL,
      if (length(.at_zero) == 0) "none" else paste(.at_zero, collapse = ", ")),
      call. = FALSE)
  }
  .no_target <- setdiff(species_present, target_strand_by_species$species)
  if (length(.no_target) > 0) {
    message(sprintf(paste0(
      "  Warning: no top-level '%s' at the target position in: %s\n",
      "  These species are drawn as recorded, without target orientation."),
      TARGET_GENE_CANONICAL, paste(.no_target, collapse = ", ")))
  }

  d <- d |> left_join(target_strand_by_species, by = "species")

  # For any species where we couldn't identify a target strand, default to '+'
  d <- d |> mutate(
    target_strand = ifelse(is.na(target_strand), "+", target_strand),
    # If ORIENT_TO_TARGET is on and target is on '-', flip every gene's
    # display strand within this species
    needs_flip = (ORIENT_TO_TARGET & target_strand == "-"),
    display_strand = ifelse(needs_flip,
                            ifelse(strand == "+", "-", "+"),
                            strand)
  )

  # ══════════════════════════════════════════════════════════════════════════════
  # COLOR PALETTE — spectral rainbow with distinct bands (jet-like but perceptual)
  # ══════════════════════════════════════════════════════════════════════════════
  unique_genes <- unique(d$gene_name)
  gene_palette <- setNames(
    hcl.colors(length(unique_genes), palette = "Spectral"),
    unique_genes
  )

  # ══════════════════════════════════════════════════════════════════════════════
  # CONTAINER-OF-TARGET RESOLUTION
  # The TSV schema records a special case: when the target gene is nested
  # inside a larger gene at the same locus, both rows sit at slot 0 with empty
  # parent_uid. The container is marked with target_nesting_gene == "y". This
  # is semantically a nesting relationship, but the parent_uid column doesn't
  # encode it because the container itself isn't nested inside anything.
  #
  # For rendering, we want the container to be a normal slot-0 arrow with the
  # target (and any siblings) tiled above it — exactly like any other
  # container/nested case. To make that happen, we synthetically link the
  # target's parent_uid to the container's feature_uid, per species.
  container_links <- d |>
    filter(slot == 0, target_nesting_gene == "y") |>
    select(species, container_uid = feature_uid) |>
    distinct(species, .keep_all = TRUE)

  if (nrow(container_links) > 0) {
    d <- d |>
      left_join(container_links, by = "species") |>
      mutate(
        # Assign the container as parent for the target and for any other
        # slot-0 gene that lacks a parent and isn't itself the container.
        # This catches Ilp4 (the target) plus any incidental slot-0 siblings.
        parent_uid = ifelse(
          !is.na(container_uid) &
            slot == 0 &
            (is.na(target_nesting_gene) | target_nesting_gene != "y") &
            (is.na(parent_uid) | parent_uid == ""),
          container_uid,
          parent_uid
        )
      ) |>
      select(-container_uid)

    n_relinked <- sum(d$slot == 0 &
                      d$target_nesting_gene != "y" &
                      d$parent_uid != "" &
                      d$parent_uid %in% container_links$container_uid,
                      na.rm = TRUE)
    if (n_relinked > 0) {
      message(sprintf(
        "  Linked %d slot-0 gene(s) to their container(s) for tiled rendering.",
        n_relinked))
    }
  }

  # ══════════════════════════════════════════════════════════════════════════════
  # SEPARATE TOP-LEVEL FROM NESTED
  # ══════════════════════════════════════════════════════════════════════════════
  top   <- d |> filter(is.na(parent_uid) | parent_uid == "")
  nest  <- d |> filter(!is.na(parent_uid) & parent_uid != "")

  nest <- nest |>
    left_join(top |> select(parent_uid = feature_uid,
                            parent_slot           = slot,
                            parent_start          = gene_start,
                            parent_end            = gene_end,
                            parent_row_y          = row_y,
                            parent_raw_strand     = strand,
                            parent_display_strand = display_strand),
              by = "parent_uid") |>
    mutate(
      # child_transcription_offset: how far into the container's transcription
      # this child sits. Computed as distance from the container's transcription
      # start (which is gene_start for '+' containers, gene_end for '-'). This
      # gives an ortholog-consistent ordering across species regardless of
      # whether the container is on + or - strand in that species.
      child_transcription_offset = ifelse(
        parent_raw_strand == "+",
        gene_start - parent_start,
        parent_end - gene_start
      )
    )

  # ══════════════════════════════════════════════════════════════════════════════
  # HELPER: build an arrow-shaped polygon
  # ══════════════════════════════════════════════════════════════════════════════
  make_arrow <- function(xmin, xmax, y_center, strand, row_height,
                         tip_frac = 0.15) {
    h  <- row_height / 2
    width <- xmax - xmin
    tip <- width * tip_frac
    if (strand == "+") {
      data.frame(
        x = c(xmin, xmax - tip, xmax,        xmax - tip, xmin),
        y = c(y_center + h, y_center + h, y_center, y_center - h, y_center - h)
      )
    } else {
      data.frame(
        x = c(xmax, xmin + tip, xmin,        xmin + tip, xmax),
        y = c(y_center + h, y_center + h, y_center, y_center - h, y_center - h)
      )
    }
  }

  # ══════════════════════════════════════════════════════════════════════════════
  # PLOT A: SLOT-SCALE
  # ══════════════════════════════════════════════════════════════════════════════
  message("Building slot-scale plot...")

  SLOT_HALF_WIDTH   <- 0.42
  ARROW_HEIGHT      <- 0.20  # every arrow, nested or not, is this tall
  NESTED_ROW_OFFSET <- 0.30  # nested arrows sit this far above their container

  top_slot_arrows <- top |>
    rowwise() |>
    mutate(
      arrow_polygon = list(make_arrow(
        xmin = slot - SLOT_HALF_WIDTH,
        xmax = slot + SLOT_HALF_WIDTH,
        y_center = row_y,
        strand = display_strand,
        row_height = ARROW_HEIGHT
      ))
    ) |>
    ungroup()

  top_slot_polys <- top_slot_arrows |>
    select(feature_uid, gene_name, arrow_polygon) |>
    unnest(cols = arrow_polygon)

  if (nrow(nest) > 0) {
    # Nested arrows all get UNIFORM width (Option A) and are TILED across the
    # container's slot span (Option D). Order is biological relative to the
    # container's display direction: for a container that's display-flipped
    # (parent_display_strand == "-"), the leftmost tile is the child that sits
    # HIGHEST in genomic coordinates (i.e. what's "upstream" in the container's
    # target-frame). For a container displaying "+", leftmost is lowest genomic.
    #
    # The container's own arrow remains on the main row at its true size; the
    # nested subrow above holds all children as evenly-sized tiles.
    nest_slot <- nest |>
      # Order tiles by transcription offset within the container. This gives
      # an ortholog-consistent ordering across species: the child that's
      # transcribed first (nearest the container's transcription start)
      # gets the lowest offset regardless of which strand the container is
      # on genomically.
      #
      # Direction of the tile row within the container:
      #   parent_display "+" (container points RIGHT in display) →
      #     leftmost tile = transcribed FIRST (ascending offset)
      #   parent_display "-" (container points LEFT in display) →
      #     leftmost tile = transcribed LAST (descending offset)
      arrange(species, parent_uid,
              if_else(parent_display_strand == "-",
                      desc(child_transcription_offset),
                      as.numeric(child_transcription_offset))) |>
      group_by(species, parent_uid) |>
      mutate(
        n_children = n(),
        child_idx  = row_number(),
        # Tile width per child, filling the container's slot span
        tile_width = (2 * SLOT_HALF_WIDTH) / n_children,
        child_xmin = parent_slot - SLOT_HALF_WIDTH +
                     (child_idx - 1) * tile_width,
        child_xmax = parent_slot - SLOT_HALF_WIDTH +
                      child_idx      * tile_width,
        nested_row_y = parent_row_y + NESTED_ROW_OFFSET
      ) |>
      ungroup() |>
      rowwise() |>
      mutate(
        arrow_polygon = list(make_arrow(
          xmin = child_xmin,
          xmax = child_xmax,
          y_center = nested_row_y,
          strand = display_strand,
          row_height = ARROW_HEIGHT
        ))
      ) |>
      ungroup()

    nest_slot_polys <- nest_slot |>
      select(feature_uid, gene_name, arrow_polygon) |>
      unnest(cols = arrow_polygon)
  } else {
    nest_slot     <- nest |> mutate(child_xmin = numeric(0), child_xmax = numeric(0))
    nest_slot_polys <- data.frame()
  }

  # Containers are drawn as full arrows (not brackets), so they stay in the
  # top-level arrow set; nested arrows are drawn on top of them.

  # ══════════════════════════════════════════════════════════════════════════════
  # UNIFIED RENDERED ARROWS
  # For ribbon and orphan calculations we care about every drawn arrow,
  # whether it's a top-level gene at its slot column or a nested gene tiled
  # above its container. This frame provides the shared shape:
  #
  #   gene_name, species, y_center, x_left, x_right, x_center, display_strand
  #
  # so the ribbon builder can connect matching genes wherever they render.
  # ══════════════════════════════════════════════════════════════════════════════
  top_rendered <- top_slot_arrows |>
    transmute(
      gene_name, species, display_strand,
      y_center = row_y,
      x_left   = slot - SLOT_HALF_WIDTH,
      x_right  = slot + SLOT_HALF_WIDTH,
      x_center = slot
    )

  if (nrow(nest) > 0) {
    nest_rendered <- nest_slot |>
      transmute(
        gene_name, species, display_strand,
        y_center = nested_row_y,
        x_left   = child_xmin,
        x_right  = child_xmax,
        x_center = (child_xmin + child_xmax) / 2
      )
    rendered_arrows_meta <- bind_rows(top_rendered, nest_rendered)
  } else {
    rendered_arrows_meta <- top_rendered
  }

  # Ribbons
  # ══════════════════════════════════════════════════════════════════════════════
  # BEZIER RIBBON HELPER
  #
  # Builds a curved ribbon polygon between two species-rows. Each ribbon has
  # a top curve and a bottom curve, sampled from cubic Bezier splines.
  #
  # For preserved orientation: top edge connects (x1_left → x2_left),
  # bottom edge connects (x1_right → x2_right). No crossing — flowing cloth.
  #
  # For flipped orientation: destination x-coordinates swap. Top edge
  # connects (x1_left → x2_right), bottom edge (x1_right → x2_left). The two
  # edges naturally cross at the midpoint, creating a twist that visually
  # encodes the inversion.
  #
  # The Bezier control points sit on horizontal tangents (control points
  # displaced vertically only), so the ribbon leaves species 1 heading
  # straight down and arrives at species 2 heading straight down. Classic
  # synteny-plot look.
  # ══════════════════════════════════════════════════════════════════════════════
  bezier_ribbon <- function(x1_left, x1_right, y1_edge,
                            x2_left, x2_right, y2_edge,
                            orientation = "preserved",
                            n_steps = 40) {
    # For flipped ribbons, swap x2's endpoints so the top edge crosses to
    # the opposite side at species 2.
    if (orientation == "flipped") {
      tmp <- x2_left
      x2_left <- x2_right
      x2_right <- tmp
    }
    # Control point vertical displacement: how far the tangent extends before
    # curving back. A larger value = smoother, longer curves.
    # Control-point vertical displacement: how far the tangent extends before
    # curving back. Larger = smoother, longer curves. Signed, so it points from
    # the start edge toward the destination edge.
    ctrl_offset <- (y2_edge - y1_edge) * 0.55
    t <- seq(0, 1, length.out = n_steps)
    bez <- function(p0, p1, p2, p3, t) {
      (1 - t)^3 * p0 + 3 * (1 - t)^2 * t * p1 +
        3 * (1 - t) * t^2 * p2 + t^3 * p3
    }
    # Both edges share the same y-curve: they differ only in x. Computing it
    # once is not just tidier, it guarantees the two edges stay in step.
    edge_y <- bez(y1_edge, y1_edge + ctrl_offset,
                  y2_edge - ctrl_offset, y2_edge, t)
    top_x  <- bez(x1_left,  x1_left,  x2_left,  x2_left,  t)
    bot_x  <- bez(x1_right, x1_right, x2_right, x2_right, t)
    # Bottom edge is reversed so the polygon closes cleanly.
    data.frame(
      x = c(top_x, rev(bot_x)),
      y = c(edge_y, rev(edge_y))
    )
  }

  build_ribbons_slot <- function(meta) {
    ribbons <- list()
    # Ribbons connect the same gene name across species. Two flavors:
    #
    # SOLID ribbons — between phylogenetically adjacent species (row_y
    #   differs by exactly 1). Full opacity, solid line.
    #
    # DRIFT ribbons — between non-adjacent species, but only where the gene
    #   is absent from every species between them. Represents "drifting in
    #   and out of the neighborhood" — the gene is nearby biologically but
    #   not stably positioned. Lower opacity, dashed line. Rendered behind
    #   solid ribbons.
    #
    # Three ribbon states are distinguished (same for solid and drift):
    #   "stayed"  — same slot, same relative orientation
    #   "moved"   — different slot, same relative orientation
    #   "flipped" — different relative orientation (drawn twisted)
    #
    # `meta` includes both top-level and nested rendered arrows, so nested
    # target genes get ribboned correctly. Species adjacency is determined
    # by each species' base row_y; nested subrows inherit their container's
    # species row_y so cross-species adjacency is preserved.
    species_row <- meta |>
      group_by(species) |>
      summarise(species_row_y = min(y_center), .groups = "drop")
    reps <- meta |>
      left_join(species_row, by = "species") |>
      select(gene_name, species, y_center, species_row_y,
             x_left, x_right, x_center, display_strand)

    # For drift detection we need the set of species row_ys that a given
    # gene appears in, so we can tell whether intervening species are
    # "empty" for this gene (they contain other genes but not this one).
    for (gene in unique(reps$gene_name)) {
      rows <- reps |> filter(gene_name == gene)
      if (nrow(rows) < 2) next
      gene_species_rows <- sort(unique(rows$species_row_y))

      for (i in seq_len(nrow(rows) - 1)) {
        for (j in (i + 1):nrow(rows)) {
          sy1 <- rows$species_row_y[i]; sy2 <- rows$species_row_y[j]
          row_gap <- abs(sy1 - sy2)
          if (row_gap == 0) next  # same species — shouldn't happen but defensive

          continuity <- if (row_gap == 1) {
            "solid"
          } else {
            # Drift candidate. Only keep it if every intervening species
            # row is truly empty for this gene. If any intermediate row
            # holds a copy of this gene, that copy will get its own
            # adjacent-solid ribbon and this longer jump is redundant.
            lower <- min(sy1, sy2)
            upper <- max(sy1, sy2)
            intermediate_rows <- seq(lower + 1, upper - 1)
            any_intermediate_has_gene <- any(intermediate_rows %in% gene_species_rows)
            if (any_intermediate_has_gene) next
            "drift"
          }

          y_upper_edge <- max(rows$y_center[i], rows$y_center[j]) - ARROW_HEIGHT / 2
          y_lower_edge <- min(rows$y_center[i], rows$y_center[j]) + ARROW_HEIGHT / 2
          same_strand <- rows$display_strand[i] == rows$display_strand[j]
          same_slot   <- rows$x_center[i] == rows$x_center[j]
          state <- if (!same_strand) {
            "flipped"
          } else if (same_slot) {
            "stayed"
          } else {
            "moved"
          }
          # Twist geometry only applies to flipped ribbons. Moved-preserved
          # ribbons are drawn with the same "preserved" ribbon topology
          # (no crossing), but the endpoints are at different x-positions
          # which produces the diagonal look naturally.
          upper_i <- if (rows$y_center[i] > rows$y_center[j]) i else j
          lower_i <- if (rows$y_center[i] > rows$y_center[j]) j else i
          ribbon_df <- bezier_ribbon(
            x1_left  = rows$x_left[upper_i],
            x1_right = rows$x_right[upper_i],
            y1_edge  = y_upper_edge,
            x2_left  = rows$x_left[lower_i],
            x2_right = rows$x_right[lower_i],
            y2_edge  = y_lower_edge,
            orientation = state
          )
          ribbon_df$gene <- gene
          ribbon_df$state <- state
          ribbon_df$continuity <- continuity
          ribbon_df$gid <- paste0(gene, "_", i, "_", j)
          ribbons[[length(ribbons) + 1]] <- ribbon_df
        }
      }
    }
    if (length(ribbons) == 0) return(data.frame())
    bind_rows(ribbons)
  }

  # Orphan X marks
  build_orphans_slot <- function(meta) {
    x_marks <- list()
    species_row <- meta |>
      group_by(species) |>
      summarise(species_row_y = min(y_center), .groups = "drop")
    reps <- meta |>
      left_join(species_row, by = "species") |>
      select(gene_name, species, species_row_y, x_center)
    for (r in seq_len(nrow(reps))) {
      gene_here <- reps$gene_name[r]
      species_here <- reps$species[r]
      sy_here <- reps$species_row_y[r]
      x_here <- reps$x_center[r]
      for (dy in c(-1, 1)) {
        sy_adj <- sy_here + dy
        if (!(sy_adj %in% reps$species_row_y)) next
        exists_adj <- any(reps$species_row_y == sy_adj & reps$gene_name == gene_here)
        if (exists_adj) next
        halfway_y <- sy_here + dy * 0.5
        x_hw <- 0.10
        y_hw <- 0.10
        x_marks[[length(x_marks) + 1]] <- data.frame(
          x_center = x_here,
          y_center = halfway_y,
          seg = c(1, 1, 2, 2),
          xend_offset = c(-x_hw, x_hw, -x_hw, x_hw),
          yend_offset = c(-y_hw, y_hw, y_hw, -y_hw),
          gene = gene_here,
          source_species = species_here,
          gid = paste0(gene_here, "_", species_here, "_dy", dy)
        )
      }
    }
    if (length(x_marks) == 0) return(data.frame())
    bind_rows(x_marks) |>
      mutate(
        x = x_center + xend_offset,
        y = y_center + yend_offset
      )
  }

  ribbons_slot <- build_ribbons_slot(rendered_arrows_meta)

  # ── Assemble slot-scale plot ──────────────────────────────────────────────────
  species_bands <- data.frame(
    ymin = seq_along(species_present) - 0.5,
    ymax = seq_along(species_present) + 0.5,
    band = seq_along(species_present) %% 2
  )

  p_slot <- ggplot()

  p_slot <- p_slot +
    geom_rect(
      data = species_bands,
      aes(xmin = -Inf, xmax = Inf, ymin = ymin, ymax = ymax,
          alpha = factor(band)),
      fill = "grey50", inherit.aes = FALSE
    ) +
    scale_alpha_manual(values = c("0" = 0, "1" = 0.05), guide = "none")

  if (nrow(ribbons_slot) > 0) {
    # Split into drift (background) and solid (foreground) so drift always
    # renders behind, and each layer has its own styling. Drift ribbons have
    # lower opacity and a dashed outline; solid ribbons are opaque fill with
    # no border.
    #
    # Both layers share the same color scale (by state), so the legend is a
    # single "Ortholog relationship" scale — the drift-vs-solid distinction
    # is called out in the subtitle instead of a second legend.
    ribbons_drift <- ribbons_slot |> filter(continuity == "drift")
    ribbons_solid <- ribbons_slot |> filter(continuity == "solid")

    if (nrow(ribbons_drift) > 0) {
      p_slot <- p_slot +
        geom_polygon(
          data = ribbons_drift,
          aes(x = x, y = y, group = gid, fill = state),
          alpha = 0.12, color = "grey40", linetype = "dashed",
          linewidth = 0.25, show.legend = TRUE
        )
    }

    if (nrow(ribbons_solid) > 0) {
      p_slot <- p_slot +
        geom_polygon(
          data = ribbons_solid,
          aes(x = x, y = y, group = gid, fill = state),
          alpha = 0.30, color = NA, show.legend = TRUE
        )
    }

    p_slot <- p_slot +
      scale_fill_manual(
        values = c("stayed"  = COLOR_PRESERVED,
                   "moved"   = COLOR_MOVED,
                   "flipped" = COLOR_FLIPPED),
        name = "Ortholog relationship",
        labels = c("stayed"  = "Same slot, same orientation",
                   "moved"   = "Different slot, same orientation",
                   "flipped" = "Different orientation")
      ) +
      ggnewscale::new_scale_fill()
  }

  # Orphan X marks are computed above but not rendered by default. The
  # absence of a ribbon between two adjacent species already communicates
  # "gene missing here." An explicit X mark added visual noise without
  # adding much information. The classification is preserved (in
  # `orphans_slot`) so a future extension could re-enable the rendering
  # or export the orphan list programmatically.
  DRAW_ORPHAN_MARKS <- FALSE
  if (DRAW_ORPHAN_MARKS) {
    orphans_slot <- build_orphans_slot(rendered_arrows_meta)
  }
  if (DRAW_ORPHAN_MARKS && nrow(orphans_slot) > 0) {
    p_slot <- p_slot +
      geom_path(
        data = orphans_slot,
        aes(x = x, y = y, group = paste0(gid, "_seg", seg)),
        color = COLOR_ORPHAN,
        linewidth = 0.9,
        alpha = 0.8,
        inherit.aes = FALSE
      )
  }

  p_slot <- p_slot +
    geom_polygon(
      data = top_slot_polys,
      aes(x = x, y = y, group = feature_uid, fill = gene_name),
      color = "black", linewidth = 0.3
    ) +
    scale_fill_manual(values = gene_palette,
                      name = "Gene",
                      guide = guide_legend(ncol = 1)) +
    ggnewscale::new_scale_fill()

  if (nrow(nest_slot_polys) > 0) {
    p_slot <- p_slot +
      geom_polygon(
        data = nest_slot_polys,
        aes(x = x, y = y, group = feature_uid, fill = gene_name),
        color = "black", linewidth = 0.25
      ) +
      scale_fill_manual(values = gene_palette, guide = "none") +
      ggnewscale::new_scale_fill()
  }

  # (Container brackets removed in this iteration: containers now render as
  # full arrows, and nested genes are drawn overlaid at their true positions
  # within the container's genomic span.)

  p_slot <- p_slot +
    geom_text(
      data = top_slot_arrows,
      aes(x = slot, y = row_y, label = gene_name),
      size = 2.0, color = "black", fontface = "italic", vjust = 0.5
    )
    # LOCID under-label disabled (uncomment for high-detail figures):
    # +
    # geom_text(
    #   data = top_slot_arrows |> filter(!is.na(locid) & locid != ""),
    #   aes(x = slot, y = row_y + LABEL_LINE2_Y_OFFSET, label = locid),
    #   size = 1.7, color = "grey30", vjust = 0.5
    # )

  # Strand indicator above the target arrow, but ONLY when the target's raw
  # strand differs from its displayed strand (i.e. the species was mirrored
  # so that the target could always point right). This tells the reader
  # "we're showing this species flipped; the target is actually on the -
  # strand in the genome."
  #
  # The target may live either in the top-level arrows (species where no
  # container-of-target exists) or in the nested arrows (species where the
  # target has been promoted to a child of its container). Look in both.
  target_name <- TARGET_GENE_CANONICAL

  # Title and scope note. A subset figure has to say what it is a subset of, or a
  # reader cannot tell it from a full plot that happens to be missing species.
  plot_title <- if (SUBSET_LEVEL == "full") {
    sprintf("Comparative microsynteny around %s", target_name)
  } else {
    sprintf("Comparative microsynteny around %s: %s %s",
            target_name, SUBSET_NAME, SUBSET_LEVEL)
  }
  subtitle_scope <- if (SUBSET_LEVEL == "full") {
    sprintf("Target: %s. All species present in the data.", target_name)
  } else {
    sprintf(paste0("Target: %s. Restricted to the %s %s; ribbon adjacency is ",
                   "recomputed within the subset."),
            target_name, SUBSET_NAME, SUBSET_LEVEL)
  }
  # Position the strand indicator above the tallest thing at the target's
  # rendered position. Because the target may have nested genes tiled above
  # it (nested subrow at parent_row_y + NESTED_ROW_OFFSET), the indicator
  # needs to clear that subrow to avoid overlapping the nested arrows.
  strand_indicator_y_offset <- NESTED_ROW_OFFSET + ARROW_HEIGHT / 2 + 0.15
  target_top_flipped <- top_slot_arrows |>
    filter(gene_name == target_name, strand != display_strand) |>
    transmute(x = slot,
              y = row_y + strand_indicator_y_offset)
  target_nest_flipped <- if (nrow(nest) > 0) {
    # Rare case: the target is itself a nested child of a container.
    # Position the indicator above the nested subrow.
    nest_slot |>
      filter(gene_name == target_name, strand != display_strand) |>
      transmute(x = (child_xmin + child_xmax) / 2,
                y = nested_row_y + ARROW_HEIGHT / 2 + 0.15)
  } else {
    data.frame(x = numeric(0), y = numeric(0))
  }
  target_flipped <- bind_rows(target_top_flipped, target_nest_flipped)
  if (nrow(target_flipped) > 0) {
    p_slot <- p_slot +
      geom_text(
        data = target_flipped,
        aes(x = x, y = y, label = "(- strand)"),
        size = 2.4, color = "grey30", fontface = "italic", vjust = 0
      )
  }

  if (nrow(nest_slot) > 0) {
    p_slot <- p_slot +
      geom_text(
        data = nest_slot,
        aes(x = (child_xmin + child_xmax) / 2,
            y = nested_row_y, label = gene_name),
        size = 1.6, color = "black", fontface = "italic"
      )
    # LOCID under-label for nested (also disabled):
    # +
    # geom_text(
    #   data = nest_slot |> filter(!is.na(locid) & locid != ""),
    #   aes(x = (child_xmin + child_xmax) / 2,
    #       y = nested_row_y - 0.06, label = locid),
    #   size = 1.4, color = "grey30"
    # )
  }

  # Target asterisk removed — the target's bold border and unique fill color
  # from the legend are enough to identify it.

  p_slot <- p_slot +
    scale_y_continuous(
      breaks = seq_along(species_present),
      labels = rev(species_present),
      expand = expansion(add = 0.9)
    ) +
    scale_x_continuous(
      breaks = sort(unique(top$slot)),
      labels = function(x) ifelse(
        x == 0, "Target",
        ifelse(x < 0, sprintf("Upstream %d", -x),
                      sprintf("Downstream %d", x))
      ),
      expand = expansion(add = 0.3)
    ) +
    labs(
      title    = plot_title,
      # Hard-wrapped onto three lines. As one long string this was silently
      # clipped at the right edge of the device, losing its last clause.
      subtitle = paste0(
        "Target-relative orientation: upstream on the left, downstream on the right. ",
        "Nested genes tile above their container.\n",
        "Ribbons: teal = same slot & orientation, blue = shifted position, ",
        "magenta = flipped (twisted). Dashed = gene drifts across intervening species.\n",
        subtitle_scope
      ),
      x        = NULL,
      y        = NULL
    ) +
    theme_minimal(base_size = 11) +
    theme(
      panel.grid         = element_blank(),
      axis.text.y        = element_text(face = "italic", size = 11),
      axis.text.x        = element_text(angle = 45, hjust = 1, size = 9),
      plot.title.position = "plot",
      plot.margin        = margin(15, 15, 20, 15),
      legend.position    = "right",
      legend.key.size    = unit(0.35, "cm"),
      legend.text        = element_text(size = 8),
      legend.title       = element_text(size = 10, face = "bold")
    )


  # ══════════════════════════════════════════════════════════════════════════════
  # SAVE
  #
  # Output filenames are prefixed with the target gene so that batch runs
  # produce one distinct plot per target rather than overwriting the same
  # file each time. Passing OUTPUT_PREFIX before source()ing overrides the
  # default, if you want a specific project label.
  # ══════════════════════════════════════════════════════════════════════════════
  h_floor <- if (SUBSET_LEVEL == "full") H_INCHES else H_INCHES_MIN_SUBSET
  h_final <- max(h_floor, length(species_present) * H_PER_SPECIES + H_HEADER)

  if (is.null(output_prefix)) {
    OUTPUT_PREFIX <- if (SUBSET_LEVEL == "full") {
      paste0(TARGET_GENE_CANONICAL, "_synteny")
    } else {
      # e.g. Ilp7_subgroup_melanogaster_synteny
      paste0(TARGET_GENE_CANONICAL, "_", SUBSET_LEVEL, "_",
             gsub("[^A-Za-z0-9._-]+", "-", SUBSET_NAME), "_synteny")
    }
  } else {
    OUTPUT_PREFIX <- output_prefix
  }

  # Default names are saved in the shared output folder.
  if (is.null(output_prefix)) {
    OUTPUT_PREFIX <- file.path(OUTPUT_DIR, OUTPUT_PREFIX)
  }

  pdf_path <- paste0(OUTPUT_PREFIX, ".pdf")
  png_path <- paste0(OUTPUT_PREFIX, ".png")
  svg_path <- paste0(OUTPUT_PREFIX, ".svg")

  ggsave(pdf_path, p_slot,
         width = W_INCHES, height = h_final, units = "in",
         limitsize = FALSE)
  ggsave(png_path, p_slot,
         width = W_INCHES, height = h_final, units = "in", dpi = DPI,
         bg = "white", limitsize = FALSE)
  # SVG is written by the svglite package. Text stays as text, so labels can
  # be edited in Inkscape or Illustrator.
  ggsave(svg_path, p_slot,
         width = W_INCHES, height = h_final, units = "in",
         bg = "white", limitsize = FALSE)

  message(sprintf("Saved three files in %s", normalizePath(dirname(pdf_path))))
  message(sprintf("  %s   (vector, for print and zooming)", basename(pdf_path)))
  message(sprintf("  %s   (vector, for the web and for editing)", basename(svg_path)))
  message(sprintf("  %s   (picture, for slides and posters)", basename(png_path)))

  invisible(p_slot)
}   # end of plot_synteny()

# ── Run from a terminal, or just load ────────────────────────────────────────
# Rscript runs this file at the top level and is not interactive, so the
# arguments on the command line are used. In RStudio (Source, or source()),
# the function is loaded and nothing is drawn until plot_synteny() is called.
if (sys.nframe() == 0L && !interactive()) {
  .args <- commandArgs(trailingOnly = TRUE)
  .flags <- .args[startsWith(.args, "--")]
  .args <- .args[!startsWith(.args, "--")]
  if (length(.flags) > 0) {
    stop(sprintf("This script takes no options; remove: %s",
                 paste(.flags, collapse = " ")), call. = FALSE)
  }
  if (length(.args) < 1) {
    cat("Usage: Rscript comparative_synteny_plot.R TSV [TARGET] [LEVEL] [NAME]\n")
    quit(save = "no", status = 1)
  }
  plot_synteny(
    tsv    = .args[1],
    target = if (length(.args) >= 2) .args[2] else NULL,
    level  = if (length(.args) >= 3) .args[3] else "full",
    name   = if (length(.args) >= 4) .args[4] else NULL
  )
} else {
  message("plot_synteny() is loaded. Example: plot_synteny(\"ImpL2_merged.tsv\")")
}
