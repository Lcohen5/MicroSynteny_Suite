#!/usr/bin/env python3
r"""
MicroSynteny Merge 1.0

Combines the MicroSynteny Table Builder files (.tsv, schema v7) for ONE target
gene into a single table for the visualizer, and writes a report of what it
did.

It uses only Python's standard library, so there is nothing to install.


HOW TO RUN IT IN JUPYTERLAB (Anaconda Navigator)
------------------------------------------------
1. Put every table for one target gene in the suite's output folder. Tables
   you built with the Table Builder are already there; add the ones other
   people built. The suggested file name is gene_lastname.tsv (for example
   ImpL2_Smith.tsv).
2. In JupyterLab's file list, open the suite's main folder (the one that
   holds table_builder, merge, visualization and output). In a notebook
   there (the same one you used for the Table Builder is fine), run a cell
   containing:

       %run merge/merge_microsynteny.py

   It shows the folder and the files it found, and asks before merging.
   (If this file has been moved away from the suite, it looks in the
   notebook's own folder instead.)

Without the questions, from a notebook or another script:

       from merge_microsynteny import merge
       merge("C:/path/to/folder")

From a terminal:

       python merge_microsynteny.py path/to/folder


WHAT COMES OUT
--------------
Two files, written into the same folder:

    <gene>_merged.tsv         the merged table, ready for the visualizer and
                              still readable by the Table Builder
    <gene>_merge_report.txt   what was read, what was kept, what to check

Your input files are never changed. Running the merge again replaces the two
output files. Files ending in _merged.tsv or _Sample.tsv are never read as
input, so the sample table that ships with the suite stays out of your merge.


WHAT THE MERGE DECIDES, AND HOW
-------------------------------
One target gene per merge.
    If the files do not all name the same target gene, the merge stops and
    lists which file names which gene. Differences in capital letters only
    (ImpL2 / impl2) are treated as the same gene; the most common spelling is
    used and the difference is reported.

One copy of each species.
    Every file normally contains D. melanogaster, and two people may have
    built the same species. For each species found in more than one file, the
    copies are compared gene by gene on: name, slot, strand, scaffold, start,
    end, and the gene it is nested inside.
      - All copies agree: one is kept.
      - Copies disagree: the version found in the most files is kept (ties go
        to the copy with more genes, then to the first file name). Every
        difference is listed in the report by file and gene. A disagreement
        is either an error or a finding; the report is there so a person can
        decide which.

Gene names are reported, never changed.
    Names that differ only in capital letters (TipE / tipE), and names that
    occur in a single species, are listed so the spelling can be checked on
    FlyBase and corrected in the Table Builder. The visualizer matches genes
    by exact name, so a misspelled gene loses its colour and its ribbon.
"""

import csv
import datetime
import io
import os
import sys
import unicodedata
from collections import Counter, OrderedDict
from typing import Dict, List, Optional, Tuple

MERGE_VERSION = "1.0"
SCHEMA_VERSION = "v7"

HEADER = [
    "schema_version", "target_gene", "species", "species_number",
    "feature_uid", "parent_uid", "gene_name", "locid", "genbank_record",
    "scaffold", "slot", "strand", "gene_start", "gene_end", "genomic_size",
    "nested_inside", "CDS_Nest", "UTR_Nest", "target_nesting_gene", "notes",
]

# A file must have these to be merged. The rest of HEADER is filled with
# blanks if absent.
REQUIRED = [
    "schema_version", "target_gene", "species", "feature_uid", "parent_uid",
    "gene_name", "slot", "strand", "gene_start", "gene_end",
]

DMEL = "D. melanogaster"
MERGED_SUFFIX = "_merged.tsv"
SAMPLE_SUFFIX = "_sample.tsv"

Row = Dict[str, str]


class MergeError(Exception):
    """The merge cannot go ahead. The message says why and what to do."""


# ----------------------------------------------------------------------
# Reading and writing (same conventions as the Table Builder)
# ----------------------------------------------------------------------

def _norm_header(h: str) -> str:
    h = unicodedata.normalize("NFKC", h or "")
    h = h.strip().lstrip("\ufeff")
    return "".join(ch for ch in h if ch.isprintable())


def _read_table(path: str) -> Tuple[List[str], List[Row]]:
    """Return (column names, rows). Every value is text; blanks stay blank."""
    with open(path, "rb") as fh:
        raw = fh.read()
    text = None
    for enc in ("utf-8-sig", "utf-8", "windows-1252", "macroman"):
        try:
            text = raw.decode(enc)
            break
        except UnicodeDecodeError:
            continue
    if text is None:
        text = raw.decode("utf-8", errors="replace")
    text = text.replace("\r\n", "\n").replace("\r", "\n")

    first_line = next((ln for ln in text.splitlines() if ln.strip()), "")
    delim = "\t" if ("\t" in first_line or "," not in first_line) else ","

    reader = csv.DictReader(io.StringIO(text, newline=""), delimiter=delim)
    if not reader.fieldnames:
        return [], []
    names = [_norm_header(h) for h in reader.fieldnames]
    field_map = dict(zip(reader.fieldnames, names))

    rows: List[Row] = []
    for rec in reader:
        row: Row = {}
        for k, v in rec.items():
            if k is None:          # more cells than columns; ignore the extras
                continue
            row[field_map.get(k, _norm_header(k))] = (v or "").strip()
        if any(row.values()):      # skip completely blank lines
            rows.append(row)
    return names, rows


def _write_table(path: str, rows: List[Row]) -> None:
    """Atomic write, in the Table Builder's own format."""
    tmp = path + ".tmp"
    with open(tmp, "w", encoding="utf-8", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=HEADER, delimiter="\t",
                           extrasaction="ignore")
        w.writeheader()
        for r in rows:
            w.writerow({h: r.get(h, "") for h in HEADER})
    os.replace(tmp, path)


# ----------------------------------------------------------------------
# Small helpers
# ----------------------------------------------------------------------

def _to_int(s: str) -> Optional[int]:
    try:
        return int(str(s).strip())
    except (TypeError, ValueError):
        try:
            f = float(str(s).strip())
        except (TypeError, ValueError):
            return None
        return int(f) if f == int(f) else None


def _safe_name(s: str) -> str:
    out = "".join(ch if (ch.isalnum() or ch in "_-") else "_" for ch in s.strip())
    return out or "merged"


def _tables_in(folder: str) -> List[str]:
    return sorted(
        (fn for fn in os.listdir(folder)
         if fn.lower().endswith(".tsv")
         and not fn.lower().endswith(MERGED_SUFFIX)
         and not fn.startswith(".")
         and os.path.isfile(os.path.join(folder, fn))),
        key=str.lower,
    )


def _find_sample_files(folder: str) -> List[str]:
    """Sample tables shipped with the suite. Listed, never merged."""
    return [fn for fn in _tables_in(folder)
            if fn.lower().endswith(SAMPLE_SUFFIX)]


def _find_input_files(folder: str) -> List[str]:
    return [fn for fn in _tables_in(folder)
            if not fn.lower().endswith(SAMPLE_SUFFIX)]


def default_folder() -> str:
    """
    The suite's shared `output` folder, when this file is in its usual place
    (<suite>/merge/ beside <suite>/output/). Anywhere else, the folder the
    notebook is running in.
    """
    try:
        here = os.path.dirname(os.path.abspath(__file__))
    except NameError:
        return os.getcwd()
    shared = os.path.join(os.path.dirname(here), "output")
    return shared if os.path.isdir(shared) else os.getcwd()


def _container_name(row: Row, by_uid: Dict[str, Row]) -> str:
    """Name of the gene this one is nested inside, or '' if top level."""
    parent = by_uid.get(row.get("parent_uid", ""))
    if parent is not None:
        return parent.get("gene_name", "")
    return row.get("nested_inside", "")


def _gene_key(row: Row, by_uid: Dict[str, Row]) -> Tuple[str, ...]:
    """What two copies of a species are compared on."""
    return (
        row.get("gene_name", ""),
        str(_to_int(row.get("slot", ""))),
        row.get("strand", ""),
        row.get("scaffold", ""),
        str(_to_int(row.get("gene_start", ""))),
        str(_to_int(row.get("gene_end", ""))),
        _container_name(row, by_uid),
    )


def _describe(key: Tuple[str, ...]) -> str:
    name, slot, strand, scaffold, start, end, container = key
    text = f"{name}  slot {slot}  strand {strand}  {scaffold}:{start}-{end}"
    if container:
        text += f"  (inside {container})"
    return text


def _file_list(files: List[str], show: int = 4) -> str:
    """'3 files: a, b, c' - shortened when a whole class agrees."""
    n = len(files)
    label = "1 file" if n == 1 else f"{n} files"
    if n <= show:
        return f"{label}: {', '.join(files)}"
    return f"{label}: {', '.join(files[:show])} and {n - show} more"


def _species_sort_key(species: str, rows: List[Row]):
    n = _to_int(rows[0].get("species_number", "")) if rows else None
    return (n if n and n > 0 else 9999, species)


# ----------------------------------------------------------------------
# The merge
# ----------------------------------------------------------------------

def merge(folder: str = ".", output: Optional[str] = None,
          quiet: bool = False) -> str:
    """
    Merge every Table Builder file in `folder` and write the merged table and
    its report there. Returns the path of the merged table.

    `output` overrides the merged table's file name. Raises MergeError, with a
    message meant to be read by the person running it, if the merge cannot go
    ahead.
    """
    folder = os.path.abspath(os.path.expanduser(folder))
    if not os.path.isdir(folder):
        raise MergeError(f"Not a folder: {folder}")

    names = _find_input_files(folder)
    samples = _find_sample_files(folder)
    if not names:
        raise MergeError(
            f"No tables to merge in {folder}\n"
            f"(files ending in {MERGED_SUFFIX} or _Sample.tsv are not read "
            "as input).")

    report: List[str] = []
    check: List[str] = []          # things a person should look at

    def say(line: str = "") -> None:
        report.append(line)

    # ---- Read ---------------------------------------------------------
    tables: "OrderedDict[str, List[Row]]" = OrderedDict()
    skipped: List[str] = [f"{fn}: a sample table, never merged" for fn in samples]
    old_schema: List[str] = []
    for fn in names:
        cols, rows = _read_table(os.path.join(folder, fn))
        if not rows:
            skipped.append(f"{fn}: no rows (an empty project)")
            continue
        missing = [c for c in REQUIRED if c not in cols]
        versions = {r.get("schema_version", "") for r in rows} - {""}
        if missing or versions != {SCHEMA_VERSION}:
            why = (f"missing column(s) {', '.join(missing)}" if missing
                   else f"schema {', '.join(sorted(versions)) or 'not recorded'}")
            old_schema.append(f"  {fn}: {why}")
            continue
        tables[fn] = rows

    if old_schema:
        raise MergeError(
            "These files are not Table Builder schema v7 files, so they "
            "cannot be merged:\n" + "\n".join(old_schema) + "\n\n"
            "If they are not Table Builder files, move them out of the "
            "folder. If they were made with an older Table Builder, they "
            "need converting to v7 first.")
    if not tables:
        raise MergeError("Every .tsv file in the folder was empty:\n  "
                         + "\n  ".join(skipped))

    # ---- One target gene ----------------------------------------------
    file_targets: Dict[str, List[str]] = {}
    for fn, rows in tables.items():
        file_targets[fn] = sorted({r["target_gene"] for r in rows
                                   if r["target_gene"]})
    no_target = [fn for fn, t in file_targets.items() if not t]
    if no_target:
        raise MergeError("These files have no target gene recorded:\n  "
                         + "\n  ".join(no_target))

    spellings = Counter(t for ts in file_targets.values() for t in ts)
    if len({t.lower() for t in spellings}) > 1:
        lines = [f"  {fn}: {', '.join(ts)}" for fn, ts in file_targets.items()]
        raise MergeError(
            "The files in this folder are not all for the same target gene, "
            "so nothing was merged.\n" + "\n".join(lines) + "\n\n"
            "Put the files for one target gene in a folder of their own and "
            "run the merge there.")

    # Most common spelling wins; ties go to the first alphabetically.
    target = sorted(spellings, key=lambda t: (-spellings[t], t))[0]
    if len(spellings) > 1:
        others = ", ".join(f"{t} ({n} file{'s' if n != 1 else ''})"
                           for t, n in sorted(spellings.items()))
        check.append(
            f"The target gene is spelled more than one way: {others}. "
            f"Using '{target}'. Check the spelling on FlyBase.")

    # Bring the target's own spelling into line everywhere, so the visualizer
    # recognises it in every species. No other gene name is touched.
    for rows in tables.values():
        for r in rows:
            r["target_gene"] = target
            if r["gene_name"].lower() == target.lower():
                r["gene_name"] = target
            if r.get("nested_inside", "").lower() == target.lower():
                r["nested_inside"] = target

    # ---- Header of the report -----------------------------------------
    say("MicroSynteny merge report")
    say("=========================")
    say(f"Merge version:  {MERGE_VERSION}")
    say(f"Date:           {datetime.datetime.now().strftime('%Y-%m-%d %H:%M')}")
    say(f"Target gene:    {target}")
    say(f"Folder:         {os.path.basename(folder)}")
    say()
    say(f"Files read ({len(tables)})")
    say("-" * 40)
    for fn, rows in tables.items():
        sp = list(OrderedDict.fromkeys(r["species"] for r in rows))
        say(f"  {fn}: {len(rows)} rows; {', '.join(sp)}")
    if skipped:
        say()
        say("Files skipped")
        say("-" * 40)
        for s in skipped:
            say(f"  {s}")

    # ---- Drop unusable rows -------------------------------------------
    bad_rows: List[str] = []
    for fn, rows in tables.items():
        keep = []
        for r in rows:
            problems = []
            if not r["species"]:
                problems.append("no species")
            if not r["gene_name"]:
                problems.append("no gene name")
            for c in ("slot", "gene_start", "gene_end"):
                if _to_int(r[c]) is None:
                    problems.append(f"{c} is '{r[c]}'")
            if problems:
                bad_rows.append(
                    f"  {fn}: {r['species'] or '?'} / {r['gene_name'] or '?'}"
                    f" left out ({'; '.join(problems)})")
            else:
                keep.append(r)
        tables[fn] = keep

    # ---- One copy of each species -------------------------------------
    # copies[species][file] = that file's rows for the species
    copies: "OrderedDict[str, OrderedDict[str, List[Row]]]" = OrderedDict()
    for fn, rows in tables.items():
        for r in rows:
            copies.setdefault(r["species"], OrderedDict()).setdefault(fn, []).append(r)

    kept: "OrderedDict[str, List[Row]]" = OrderedDict()
    agree_lines: List[str] = []
    disagree_lines: List[str] = []

    for species, by_file in copies.items():
        if len(by_file) == 1:
            kept[species] = next(iter(by_file.values()))
            continue

        # Group the copies into versions by content.
        versions: "OrderedDict[Tuple, List[str]]" = OrderedDict()
        keysets: Dict[str, Counter] = {}
        for fn, rows in by_file.items():
            by_uid = {r["feature_uid"]: r for r in rows if r["feature_uid"]}
            ks = Counter(_gene_key(r, by_uid) for r in rows)
            keysets[fn] = ks
            versions.setdefault(tuple(sorted(ks.items())), []).append(fn)

        ranked = sorted(
            versions.values(),
            key=lambda fns: (-len(fns), -len(by_file[fns[0]]), fns[0].lower()))
        winner_files = ranked[0]
        winner = winner_files[0]
        kept[species] = by_file[winner]

        if len(versions) == 1:
            agree_lines.append(
                f"  {species}: {len(by_file)} copies, all in agreement. "
                f"Kept the copy in {winner}.")
            continue

        disagree_lines.append(
            f"  {species}: {len(by_file)} copies in {len(versions)} versions.")
        tie = len(ranked) > 1 and len(ranked[1]) == len(winner_files)
        disagree_lines.append(
            f"    KEPT the version in {_file_list(winner_files)}")
        if tie:
            disagree_lines.append(
                "    TIE: no version is in more files than the other. This one "
                "was kept only for having more genes or the earlier file "
                "name. Please check which is right.")
        for other_files in ranked[1:]:
            other = other_files[0]
            only_other = keysets[other] - keysets[winner]
            only_kept = keysets[winner] - keysets[other]
            disagree_lines.append(
                f"    NOT USED, the version in {_file_list(other_files)}")
            for key in sorted(only_other.elements()):
                disagree_lines.append(f"        this version has:  {_describe(key)}")
            for key in sorted(only_kept.elements()):
                disagree_lines.append(f"        kept version has:  {_describe(key)}")
        disagree_lines.append("")

    say()
    say("Species found in more than one file")
    say("-" * 40)
    if not agree_lines and not disagree_lines:
        say("  None. Every species came from exactly one file.")
    for line in agree_lines:
        say(line)
    if disagree_lines:
        say()
        say("  DISAGREEMENTS (each is an error or a finding; please decide which)")
        for line in disagree_lines:
            say(line)
        n = sum(1 for line in disagree_lines if " copies in " in line)
        check.append(f"{n} species had copies that disagree. See "
                     "'Species found in more than one file'.")

    # ---- Assemble, in phylogenetic order ------------------------------
    merged: List[Row] = []
    for species in sorted(kept, key=lambda s: _species_sort_key(s, kept[s])):
        merged.extend(kept[species])

    # ---- Checks on the merged table -----------------------------------
    notes: List[str] = list(bad_rows)

    # A nested gene must point at a gene in its own species' kept copy.
    for species, rows in kept.items():
        uids = {r["feature_uid"] for r in rows}
        for r in rows:
            if r["parent_uid"] and r["parent_uid"] not in uids:
                notes.append(
                    f"  {species} / {r['gene_name']}: marked as nested, but its "
                    "container gene is not in the table. Treated as top level.")
                r["parent_uid"] = ""

    dup_uids = [u for u, n in Counter(r["feature_uid"] for r in merged).items()
                if u and n > 1]
    if dup_uids:
        notes.append(
            f"  {len(dup_uids)} gene ID(s) occur more than once in the merged "
            "table. This usually means one species' rows were copied by hand "
            "into another species. The figure may draw these genes wrongly.")

    # The visualizer needs the target at slot 0 in every species.
    for species, rows in kept.items():
        at_zero = [r["gene_name"] for r in rows if _to_int(r["slot"]) == 0]
        if target not in at_zero:
            found = ", ".join(sorted(set(at_zero))) or "nothing"
            notes.append(
                f"  {species}: no gene named {target} at the target position "
                f"(found there: {found}). The figure cannot orient this "
                "species until that is corrected in the Table Builder.")

    if notes:
        say()
        say("Problems in the tables")
        say("-" * 40)
        for line in notes:
            say(line)
        check.append(f"{len(notes)} problem(s) in the tables. See 'Problems in "
                     "the tables'.")

    # ---- Gene names to check ------------------------------------------
    species_of: Dict[str, set] = {}
    for r in merged:
        species_of.setdefault(r["gene_name"], set()).add(r["species"])

    by_lower: Dict[str, List[str]] = {}
    for name in species_of:
        by_lower.setdefault(name.lower(), []).append(name)
    case_lines = []
    in_case_group = set()
    for group in sorted(by_lower.values(), key=lambda g: g[0].lower()):
        if len(group) > 1:
            in_case_group.update(group)
            parts = [f"{n} ({len(species_of[n])} species: "
                     f"{', '.join(sorted(species_of[n]))})"
                     for n in sorted(group, key=lambda n: -len(species_of[n]))]
            case_lines.append("  " + "  vs  ".join(parts))

    single_lines = []
    if len(kept) >= 3:
        for name in sorted(species_of, key=str.lower):
            if len(species_of[name]) == 1 and name not in in_case_group:
                single_lines.append(f"  {name}  ({next(iter(species_of[name]))})")

    say()
    say("Gene names to check on FlyBase")
    say("-" * 40)
    if not case_lines and not single_lines:
        say("  Nothing to check.")
    if case_lines:
        say("  Same name with different capital letters. The figure treats")
        say("  these as different genes:")
        for line in case_lines:
            say("  " + line)
    if single_lines:
        if case_lines:
            say()
        say("  Names found in only one species. Each is either a gene that")
        say("  really is unique to that species here, or a misspelling:")
        for line in single_lines:
            say("  " + line)
    if case_lines or single_lines:
        say()
        say("  Nothing was changed. Confirm each spelling on FlyBase")
        say("  (flybase.org) and correct it in the Table Builder if needed,")
        say("  then run the merge again.")
        check.append(f"{len(case_lines) + len(single_lines)} gene name(s) to "
                     "check on FlyBase.")

    # ---- Write --------------------------------------------------------
    out_name = output or f"{_safe_name(target)}{MERGED_SUFFIX}"
    out_path = out_name if os.path.isabs(out_name) else os.path.join(folder, out_name)
    report_path = os.path.join(
        os.path.dirname(out_path), f"{_safe_name(target)}_merge_report.txt")

    summary = [
        "Result",
        "-" * 40,
        f"  Merged table:  {os.path.basename(out_path)}",
        f"  Rows:          {len(merged)}",
        f"  Species ({len(kept)}):  "
        + ", ".join(sorted(kept, key=lambda s: _species_sort_key(s, kept[s]))),
        "",
        "To check" if check else "Nothing needs checking.",
    ]
    summary += [f"  - {c}" for c in check]

    # The result and the to-do list go at the top of the report, where the
    # reader starts; the detail follows.
    head, body = report[:7], report[7:]
    full = head + summary + [""] + body

    _write_table(out_path, merged)
    with open(report_path, "w", encoding="utf-8", newline="\n") as fh:
        fh.write("\n".join(full).rstrip() + "\n")

    if not quiet:
        print("\n".join(full))
        print()
        print(f"Saved: {os.path.basename(out_path)}")
        print(f"Saved: {os.path.basename(report_path)}")
        print(f"In:    {folder}")
    return out_path


# ----------------------------------------------------------------------
# Running it
# ----------------------------------------------------------------------

def _ask(prompt: str) -> Optional[str]:
    try:
        return input(prompt)
    except EOFError:
        return None
    except Exception as e:
        if type(e).__name__ == "StdinNotImplementedError":
            return None
        raise


def run(folder: Optional[str] = None) -> Optional[str]:
    """
    Run the merge with questions, for a notebook cell. Returns the merged
    table's path, or None if nothing was merged. Never raises on an expected
    problem and never exits the kernel.
    """
    print(f"MicroSynteny Merge {MERGE_VERSION}")
    try:
        if folder is None:
            folder = default_folder()
            while True:
                names = _find_input_files(folder) if os.path.isdir(folder) else []
                print(f"\nFolder: {folder}")
                if names:
                    print(f"Files to merge ({len(names)}):")
                    for fn in names:
                        print(f"  {fn}")
                else:
                    print("No .tsv files to merge here.")
                ans = _ask("Merge these files? (y = yes, n = choose another "
                           "folder, q = quit) ")
                if ans is None:
                    # No input box available: go ahead with this folder.
                    break
                ans = ans.strip().lower()
                if ans in ("y", "yes") and names:
                    break
                if ans in ("q", ":q", "quit", "exit"):
                    print("Nothing merged.")
                    return None
                if ans in ("n", "no") or not names:
                    raw = _ask("Type the full path to the folder: ")
                    if raw is None or not raw.strip():
                        print("Nothing merged.")
                        return None
                    folder = os.path.abspath(os.path.expanduser(
                        raw.strip().strip('"').strip("'")))
        print()
        return merge(folder)
    except MergeError as e:
        print(f"\nMERGE STOPPED\n{e}")
        return None
    except KeyboardInterrupt:
        print("\nStopped. Nothing merged.")
        return None


if __name__ == "__main__":
    _args = [a for a in sys.argv[1:] if not a.startswith("-")]
    # Under %run in a notebook sys.argv holds only the script name, so this
    # asks its questions; from a terminal with a folder argument it does not.
    if _args:
        try:
            merge(_args[0])
        except MergeError as _e:
            print(f"\nMERGE STOPPED\n{_e}")
            if "ipykernel" not in sys.modules:
                sys.exit(1)
    else:
        run()
