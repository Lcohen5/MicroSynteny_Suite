#!/usr/bin/env python3
r"""
MicroSynteny Table Builder V12.J - JupyterLab edition (Anaconda Navigator)

Same program as MicroSynteny_Table_Builder_V12.py. Only the input/output layer
differs; the parser, the entry flows and the correction editor are identical
code, so the two editions cannot drift apart.

Schema version: v7 (UNCHANGED). V9, V10 and V11 project files open directly.


HOW TO RUN IT IN ANACONDA NAVIGATOR
-----------------------------------
1. Open Anaconda Navigator and launch JupyterLab.
2. In JupyterLab's file list, open the suite's main folder (the one that
   holds table_builder, merge, visualization and output).
3. Make a new notebook there (Python 3) and put this in the first cell:

       %run table_builder/MicroSynteny_Table_Builder.py

   Run the cell. Answer the questions in the input box that appears under it.
   Your tables are saved in the suite's output folder, next to table_builder.
   (If this file has been moved away from the suite, tables are saved in the
   notebook's own folder instead. Either way the program shows the folder and
   asks before using it.)

   If you would rather keep the notebook somewhere else:

       import sys; sys.path.append("C:/path/to/folder")
       from MicroSynteny_Table_Builder import run
       run(folder="C:/path/to/folder")

   Forward slashes work in Windows paths in Python, and save you
   worrying about backslashes.

4. Type :q at any prompt to save and quit. Your work is saved after every
   change regardless, so an interrupted kernel loses nothing.


ABOUT PASTING IN JUPYTERLAB
---------------------------
JupyterLab's input box is a single-line field. When you paste several lines
into it the browser strips the line breaks, so the gene information block
arrives as one run-on string:

    Genbank Record: NC_004354.4Position: chr2L:10,242,512-10,247,078Strand: +

That is exactly what made pasting painful in earlier versions. V11 parses
labelled fields wherever they appear rather than line by line, so the run-on
form is read correctly and the paste ends by itself - no blank line needed.

If a paste ever does come out wrong, two escapes are always available at a
paste prompt:

    clip     read the block straight from your system clipboard, so the text
             never travels through the input box at all. (This reads the
             clipboard of the machine running the kernel, which for a normal
             Anaconda install is your own computer.)
    cancel   back out of this gene

And when a paste is missing a field, the program asks only for the field it
could not find instead of making you paste everything again.

V12.1 fix: NCBI blocks now end with "CDS Start: complete" and "CDS End:
complete" after the Gene Symbol line. In a run-on paste nothing terminated the
symbol value, so those lines were swallowed into the LOCID. The parser now
recognises them as labels and, independently, keeps only the first token of
the symbol. Project files built before this fix may carry the extra text in
their locid column; see clean_locid_column() at the bottom of this file.


NOTEBOOK-SPECIFIC NOTES
-----------------------
- Set CLEAR_BETWEEN_STEPS = True below if you would rather the cell redrew
  itself each time instead of growing a long scroll. It is off by default
  because the scroll is a record of what you did.
- The program never calls sys.exit(); quitting returns control to the cell so
  the kernel and your variables stay alive.
- Interrupting the kernel (the stop button) saves before it stops.
"""

import csv
import io
import os
import re
import shutil
import subprocess
import sys
import unicodedata
from dataclasses import dataclass, asdict
from typing import Callable, Dict, List, Optional, Tuple
from uuid import uuid4

SCHEMA_VERSION = "v7"
PROGRAM_VERSION = "V12"

SPECIES_FULL = [
    "D. melanogaster",
    "D. mauritiana",
    "D. sechellia",
    "D. simulans",
    "D. yakuba",
    "D. santomea",
    "D. teissieri",
    "D. erecta",
    "D. ficusphila",
    "D. suzukii",
    "D. subpulchrella",
    "D. biarmipes",
    "D. takahashii",
    "D. eugracilis",
    "D. rhopaloa",
    "D. elegans",
    "D. kikkawai",
    "D. serrata",
    "D. bipectinata",
    "D. ananassae",
    "D. pseudoobscura",
    "D. persimilis",
    "D. miranda",
    "D. guanche",
    "D. subobscura",
    "D. obscura",
    "D. willistoni",
    "D. arizonae",
    "D. mojavensis",
    "D. navojoa",
    "D. hydei",
    "D. virilis",
    "D. novamexicana",
    "D. albomicans",
    "D. grimshawi",
    "D. busckii",
]

DMEL = "D. melanogaster"


def phylogenetic_number(species: str) -> int:
    """1-based index of `species` in SPECIES_FULL (phylogenetic order); 0 if unknown."""
    try:
        return SPECIES_FULL.index(species) + 1
    except ValueError:
        return 0


HEADER = [
    "schema_version",
    "target_gene",
    "species",
    "species_number",
    "feature_uid",
    "parent_uid",
    "gene_name",
    "locid",
    "genbank_record",
    "scaffold",
    "slot",
    "strand",
    "gene_start",
    "gene_end",
    "genomic_size",
    "nested_inside",
    "CDS_Nest",
    "UTR_Nest",
    "target_nesting_gene",
    "notes",
]

QUIT_TOKENS = {":q", ":wq", "save and quit", "quit", "exit"}
BACK_TOKENS = {"back", ":b", "go back", "previous"}
END_TOKENS = {"end", ":e", ".", ";;", "done"}
MANUAL_TOKENS = {"manual"}
CLIP_TOKENS = {"clip", "clipboard", ":c"}
CANCEL_TOKENS = {"cancel", "skip", ":x"}


@dataclass
class Feature:
    target_gene: str
    species: str
    species_number: int
    feature_uid: str
    parent_uid: str
    gene_name: str
    locid: str
    genbank_record: str
    scaffold: str
    slot: int
    strand: str
    gene_start: int
    gene_end: int
    genomic_size: int
    nested_inside: str = ""
    CDS_Nest: str = ""
    UTR_Nest: str = ""
    target_nesting_gene: str = "n"
    notes: str = ""
    schema_version: str = SCHEMA_VERSION


class QuitSignal(Exception):
    """Raised by any prompt when the student types a quit token."""


class ParseError(Exception):
    """Raised when pasted text cannot be turned into a gene, or entry is aborted."""


class ProjectError(Exception):
    """A project file cannot be opened. Reported to the student, not raised as SystemExit."""


# Sentinel returned by back-aware prompts when the student wants to step back.
BACK = object()


# ==== ENV BEGIN ====
# Everything between the ENV markers is the only environment-specific code in
# this program. The terminal edition (V11) is identical except for this block,
# the module docstring and the footer.

ENV_NAME = "jupyterlab"
ENV_IS_JUPYTER = True

# Redraw the cell instead of letting it scroll. Off by default: the scroll is a
# useful record of what was entered.
CLEAR_BETWEEN_STEPS = False


def _in_ipython() -> bool:
    try:
        from IPython import get_ipython

        return get_ipython() is not None
    except Exception:
        return False


# Notebook output is UTF-8, so the arrows and box-drawing characters are safe.
UNICODE_OK = True

GLYPH = {
    "right": "\u2192",
    "left": "\u2190",
    "unknown": "?",
    "branch": "\u2514\u2500 ",
    "bar": "\u2502",
}


def clear_screen() -> None:
    if not CLEAR_BETWEEN_STEPS:
        return
    try:
        from IPython.display import clear_output

        clear_output(wait=True)
    except Exception:
        pass


def read_line(prompt: str = "") -> Optional[str]:
    """
    One line of input from the notebook's stdin box.

    Returns None on EOF, which the callers treat as a quit. A kernel with no
    stdin at all - papermill, nbconvert, "Run All" in some setups - raises
    StdinNotImplementedError, which is worth an explicit message because the
    program is interactive by nature and cannot run unattended.
    """
    try:
        return input(prompt)
    except EOFError:
        return None
    except Exception as e:
        if type(e).__name__ == "StdinNotImplementedError":
            print(
                "\nThis notebook cannot accept typed input, so the builder cannot "
                "run here.\nRun the cell interactively (click into the notebook and "
                "press Shift+Enter)\nrather than through nbconvert, papermill or a "
                "headless 'Run All'."
            )
            return None
        raise


def clipboard_text() -> Optional[str]:
    """
    System clipboard contents, or None if it cannot be read.

    This reads the clipboard of the machine running the KERNEL. For a normal
    local Anaconda install that is the student's own computer. On a remote or
    hosted kernel it is not, so paste normally there.

    tkinter is deliberately not used as a fallback here: creating a Tk root
    inside a kernel can block it, and a hung kernel is a worse outcome than a
    clipboard read that politely fails.
    """
    cmds: List[List[str]] = []
    if sys.platform == "darwin":
        cmds.append(["pbpaste"])
    elif sys.platform.startswith("win"):
        cmds.append([
            "powershell", "-NoProfile", "-Command",
            "[Console]::OutputEncoding=[System.Text.Encoding]::UTF8; Get-Clipboard",
        ])
    else:
        cmds.append(["xclip", "-selection", "clipboard", "-o"])
        cmds.append(["xsel", "--clipboard", "--output"])
        cmds.append(["wl-paste", "--no-newline"])

    for cmd in cmds:
        if shutil.which(cmd[0]) is None:
            continue
        try:
            out = subprocess.run(cmd, capture_output=True, timeout=10)
        except Exception:
            continue
        if out.returncode == 0:
            text = out.stdout.decode("utf-8", errors="replace")
            if text.strip():
                return text
    return None


PASTE_HELP = (
    "  Paste the block into the input box and press Enter. If the lines run "
    "together\n  into one long string, that is fine - it is read correctly "
    "either way.\n"
    "  Shortcuts: 'clip' reads it straight from your clipboard, 'cancel' backs "
    "out,\n  ':q' saves and quits."
)

BANNER = (
    f"MicroSynteny Table Builder {PROGRAM_VERSION}.J (JupyterLab edition)\n"
    "Answer the prompts in the input box below. Type :q at any prompt to save "
    "and quit."
)
# ==== ENV END ====


# ----------------------------------------------------------------------
# Prompt helpers
# ----------------------------------------------------------------------

def is_quit(s: str) -> bool:
    return s.strip().lower() in QUIT_TOKENS


def is_back(s: str) -> bool:
    return s.strip().lower() in BACK_TOKENS


def is_manual(s: str) -> bool:
    return s.strip().lower() in MANUAL_TOKENS


def is_cancel(s: str) -> bool:
    return s.strip().lower() in CANCEL_TOKENS


def warn(msg: str) -> None:
    print(f"  ! {msg}")


def note(msg: str) -> None:
    print(f"  - {msg}")


def prompt_raw(prompt: str, allow_empty: bool = False) -> str:
    while True:
        s = read_line(prompt)
        if s is None:
            raise QuitSignal()
        s = s.strip()
        if is_quit(s):
            raise QuitSignal()
        if s == "" and not allow_empty:
            note("This one needs a value. (Type :q to save and quit.)")
            continue
        return s


def prompt_int(prompt: str, minv: Optional[int] = None, maxv: Optional[int] = None) -> int:
    while True:
        s = prompt_raw(prompt)
        try:
            v = int(s.replace(",", ""))
        except ValueError:
            warn("Please enter a whole number.")
            continue
        if minv is not None and v < minv:
            warn(f"Must be at least {minv}.")
            continue
        if maxv is not None and v > maxv:
            warn(f"Must be at most {maxv}.")
            continue
        return v


def prompt_yesno(prompt: str) -> bool:
    while True:
        s = prompt_raw(prompt).lower()
        if s in {"y", "yes"}:
            return True
        if s in {"n", "no"}:
            return False
        warn("Please enter y or n.")


def prompt_yesno_cancel(prompt: str) -> Optional[bool]:
    """
    Three-way choice: True (yes), False (no), None (cancel).

    Use this at gates where a mistaken y or n commits the student to a
    chain of downstream data entry - e.g. "is X nested inside another
    gene?" answering yes triggers container-gene collection, and answering
    no locks the gene in as top-level.

    The caller decides what "cancel" means in context (usually: return to
    an earlier prompt without committing).
    """
    while True:
        s = prompt_raw(prompt).lower()
        if s in {"y", "yes"}:
            return True
        if s in {"n", "no"}:
            return False
        if s in {"c", "cancel", "back", "wait"}:
            return None
        warn("Please enter y, n, or 'cancel'.")


def prompt_choice(prompt: str, options: Dict[str, str]) -> str:
    """
    Choose one of a small set of options. Accepts the number as shown, the
    full keyword, or an unambiguous prefix. Returns the chosen keyword.
    """
    keys = list(options.keys())
    print("\n" + prompt)
    for i, (k, desc) in enumerate(options.items(), start=1):
        print(f"   {i}) {k:<10} {desc}")
    number_hint = "/".join(str(i) for i in range(1, len(keys) + 1))
    while True:
        s = prompt_raw(f"Choose [{number_hint}]: ").lower()
        if s.isdigit():
            n = int(s)
            if 1 <= n <= len(keys):
                return keys[n - 1]
            warn(f"Enter a number between 1 and {len(keys)}.")
            continue
        if s in options:
            return s
        matches = [k for k in keys if k.startswith(s)]
        if len(matches) == 1:
            return matches[0]
        warn(f"Enter a number ({number_hint}) or one of: {', '.join(keys)}.")


def prompt_raw_back(prompt: str, allow_empty: bool = False):
    """Like prompt_raw, but 'back' returns the BACK sentinel."""
    while True:
        s = read_line(prompt)
        if s is None:
            raise QuitSignal()
        s = s.strip()
        if is_quit(s):
            raise QuitSignal()
        if is_back(s):
            return BACK
        if s == "" and not allow_empty:
            note("This one needs a value. (Type 'back' to step back, or :q to "
                 "save and quit.)")
            continue
        return s


def prompt_field_value(prompt: str, allow_manual: bool = True):
    """Single free-text field. Returns BACK, the literal 'Manual', or the text."""
    s = prompt_raw_back(prompt, allow_empty=False)
    if s is BACK:
        return BACK
    if allow_manual and is_manual(s):
        return "Manual"
    return s


def prompt_int_field(prompt: str, minv: Optional[int] = None, allow_manual: bool = True):
    """Integer field that may be 'Manual' or 'back'. Returns BACK, 'Manual', or int."""
    while True:
        s = prompt_raw_back(prompt, allow_empty=False)
        if s is BACK:
            return BACK
        if allow_manual and is_manual(s):
            return "Manual"
        try:
            v = int(s.replace(",", "").replace(" ", ""))
        except ValueError:
            warn("Please enter a whole number (or 'Manual', or 'back').")
            continue
        if minv is not None and v < minv:
            warn(f"Must be at least {minv}.")
            continue
        return v


Step = Tuple[str, Callable]


def run_steps(steps: List[Step], state: Dict, back_out: bool = False):
    """
    Run a linear sequence of steps with 'go back a step' support. Each step fn
    returns True (advance), False (abort the sequence) or BACK (step back).

    'back' at the first step normally has nowhere to go. With back_out=True it
    leaves the sequence instead, and run_steps returns BACK so the caller can
    return the student to whatever came before.
    """
    i = 0
    while i < len(steps):
        _name, fn = steps[i]
        result = fn(state)
        if result is BACK:
            if i == 0:
                if back_out:
                    return BACK
                note("You're at the first question; there's nothing to go back to.")
                continue
            i -= 1
            continue
        if result is False:
            return False
        i += 1
    return True


def fmt_int(n: int) -> str:
    try:
        return f"{int(n):,}"
    except (TypeError, ValueError):
        return str(n)


# ----------------------------------------------------------------------
# TSV read / write
# ----------------------------------------------------------------------

def sanitize_filename(s: str) -> str:
    s = re.sub(r"\s+", "_", s.strip())
    s = re.sub(r"[^A-Za-z0-9_\-]+", "", s)
    return s or "project"


def project_filename(project_name: str) -> str:
    return sanitize_filename(project_name) + ".tsv"


def list_project_files(folder: str) -> List[str]:
    # The sample table that ships with the suite (ImpL2_Sample.tsv) shares the
    # output folder with everyone's projects. It is left out of the list so it
    # cannot be opened and changed by accident.
    return sorted(
        fn
        for fn in os.listdir(folder)
        if fn.lower().endswith(".tsv")
        and not fn.lower().endswith("_sample.tsv")
        and os.path.isfile(os.path.join(folder, fn))
    )


def _norm_header(h: str) -> str:
    h = unicodedata.normalize("NFKC", h)
    h = h.strip().lstrip("\ufeff")
    return "".join(ch for ch in h if ch.isprintable())


def _detect_delimiter(sample: str) -> str:
    if "\t" in sample:
        return "\t"
    if "," in sample:
        return ","
    return "\t"


def repair_to_utf8(path: str) -> None:
    b = open(path, "rb").read()
    for enc in ("utf-8-sig", "utf-8", "windows-1252", "macroman"):
        try:
            s = b.decode(enc)
        except UnicodeDecodeError:
            continue
        with open(path, "w", encoding="utf-8", newline="") as fh:
            fh.write(s)
        return
    s = b.decode("utf-8", errors="replace")
    with open(path, "w", encoding="utf-8", newline="") as fh:
        fh.write(s)


def _new_uid() -> str:
    return uuid4().hex


def write_tsv(path: str, feats: List[Feature]) -> None:
    """Atomic write: build a .tmp beside the target, then replace."""
    tmp = path + ".tmp"
    with open(tmp, "w", encoding="utf-8", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=HEADER, delimiter="\t")
        w.writeheader()
        for f in feats:
            d = asdict(f)
            d["slot"] = int(d["slot"])
            d["gene_start"] = int(d["gene_start"])
            d["gene_end"] = int(d["gene_end"])
            d["genomic_size"] = int(d["genomic_size"])
            w.writerow(d)
    os.replace(tmp, path)


OLD_PROJECT_MESSAGE = (
    "This project was created with an older schema version. V11 uses schema "
    "v7 (biological slot semantics: upstream is always negative, downstream "
    "always positive, regardless of target strand). Older schemas (v4/v5/v6) "
    "stored slot in genomic-coordinate space and are NOT auto-migrated on "
    "open, because a silent migration could corrupt data.\n\n"
    "To migrate old files:\n"
    "    python migrate_v6_to_v7.py <folder_containing_tsvs>\n\n"
    "The migration script flips slot signs for species where the target is on "
    "the - strand and bumps schema_version to v7. Old files are backed up "
    "before rewriting."
)


def _repair_container_links(feats: List[Feature]) -> int:
    """
    V10 recorded a container gene by setting the target's `nested_inside` name
    but never setting `parent_uid`, so both sat at the top level in the same
    slot and only one was ever drawn. Re-link any such orphan to its named
    container. Returns the number of rows repaired.
    """
    fixed = 0
    by_species: Dict[str, List[Feature]] = {}
    for f in feats:
        by_species.setdefault(f.species, []).append(f)

    for _species, group in by_species.items():
        for f in group:
            if f.parent_uid.strip() or not f.nested_inside.strip():
                continue
            candidates = [
                g
                for g in group
                if g.gene_name == f.nested_inside
                and g.feature_uid != f.feature_uid
                and g.slot == f.slot
            ]
            if len(candidates) == 1:
                f.parent_uid = candidates[0].feature_uid
                fixed += 1
    return fixed


def read_tsv(path: str) -> List[Feature]:
    """Read a v7 TSV. Aborts with the old-project message on any other schema."""
    b = open(path, "rb").read()

    def decode_bytes() -> str:
        for enc in ("utf-8-sig", "utf-8", "windows-1252", "macroman"):
            try:
                return b.decode(enc)
            except UnicodeDecodeError:
                continue
        return b.decode("utf-8", errors="replace")

    text = decode_bytes().replace("\r\n", "\n").replace("\r", "\n")
    first_line = next((ln for ln in text.splitlines() if ln.strip()), "")
    delim = _detect_delimiter(first_line)

    fh = io.StringIO(text, newline="")
    r = csv.DictReader(fh, delimiter=delim)
    if not r.fieldnames:
        raise ProjectError(f"Could not read a header row from {os.path.basename(path)}.")

    norm_fields = [_norm_header(h) for h in r.fieldnames]
    if "schema_version" not in norm_fields:
        raise ProjectError(OLD_PROJECT_MESSAGE)
    field_map = dict(zip(r.fieldnames, norm_fields))

    rows: List[Dict[str, str]] = []
    for row in r:
        nrow: Dict[str, str] = {}
        for k, v in row.items():
            nrow[field_map.get(k, _norm_header(k))] = v
        rows.append(nrow)

    if not rows:
        return []

    versions = {(row.get("schema_version") or "").strip() for row in rows}
    versions.discard("")
    if versions and SCHEMA_VERSION not in versions:
        raise ProjectError(OLD_PROJECT_MESSAGE)

    def s(nrow: Dict[str, str], key: str, default: str = "") -> str:
        return (nrow.get(key, default) or default).strip()

    feats: List[Feature] = []
    for nrow in rows:
        species_name = s(nrow, "species")
        stored_number = int(s(nrow, "species_number", "0") or 0)
        if stored_number == 0 and species_name:
            stored_number = phylogenetic_number(species_name)
        try:
            feats.append(
                Feature(
                    target_gene=s(nrow, "target_gene"),
                    species=species_name,
                    species_number=stored_number,
                    feature_uid=s(nrow, "feature_uid") or _new_uid(),
                    parent_uid=s(nrow, "parent_uid"),
                    gene_name=s(nrow, "gene_name"),
                    locid=s(nrow, "locid"),
                    genbank_record=s(nrow, "genbank_record"),
                    scaffold=s(nrow, "scaffold"),
                    slot=int(s(nrow, "slot", "0") or 0),
                    strand=s(nrow, "strand"),
                    gene_start=int(s(nrow, "gene_start", "0") or 0),
                    gene_end=int(s(nrow, "gene_end", "0") or 0),
                    genomic_size=int(s(nrow, "genomic_size", "0") or 0),
                    nested_inside=s(nrow, "nested_inside"),
                    CDS_Nest=s(nrow, "CDS_Nest"),
                    UTR_Nest=s(nrow, "UTR_Nest"),
                    target_nesting_gene=s(nrow, "target_nesting_gene", "n") or "n",
                    notes=s(nrow, "notes"),
                    schema_version=SCHEMA_VERSION,
                )
            )
        except ValueError:
            warn(f"Skipping a malformed row: {nrow}")
            continue

    repaired = _repair_container_links(feats)
    if repaired:
        note(
            f"Repaired {repaired} container link(s) written by an older version "
            "so that container genes appear in the layout."
        )
    return feats


# ----------------------------------------------------------------------
# Folder + project selection
# ----------------------------------------------------------------------

def choose_from_list(prompt: str, options: List[str], allow_back: bool = False):
    """
    Pick one of `options` by number. With allow_back=True the student may type
    'back' instead, and BACK is returned so the caller can show the previous
    question again.
    """
    if not options:
        raise ProjectError("There is nothing to choose from.")
    print("\n" + prompt)
    for i, opt in enumerate(options, start=1):
        print(f"  {i:2d}. {opt}")
    ask = "Enter number (or 'back'): " if allow_back else "Enter number: "
    while True:
        raw = prompt_raw(ask)
        if allow_back and is_back(raw):
            return BACK
        try:
            idx = int(raw.replace(",", ""))
        except ValueError:
            warn("Please enter a whole number.")
            continue
        if 1 <= idx <= len(options):
            return options[idx - 1]
        warn("That number isn't on the list.")


def choose_species() -> str:
    return choose_from_list("Select species:", SPECIES_FULL)


def species_in_table(feats: List[Feature]) -> List[str]:
    present = {f.species for f in feats if (f.species or "").strip()}
    return sorted(present, key=lambda s: (phylogenetic_number(s) or 999, s))


def default_project_folder() -> str:
    """
    Where tables are kept by default: the suite's shared `output` folder, when
    this file is in its usual place (<suite>/table_builder/ beside
    <suite>/output/). Anywhere else, the folder the notebook is running in.
    """
    try:
        here = os.path.dirname(os.path.abspath(__file__))
    except NameError:
        return os.getcwd()
    shared = os.path.join(os.path.dirname(here), "output")
    return shared if os.path.isdir(shared) else os.getcwd()


def choose_folder() -> str:
    folder = default_project_folder()
    print(f"\nProject folder: {folder}")
    existing = list_project_files(folder) if os.path.isdir(folder) else []
    if existing:
        note(f"{len(existing)} table(s) here.")
    else:
        note("No tables here yet.")
    if prompt_yesno("Use this folder? (y/n) "):
        return folder
    while True:
        raw = prompt_raw("Type the full path to the folder you want: ")
        cand = os.path.abspath(os.path.expanduser(raw.strip().strip('"').strip("'")))
        if os.path.isdir(cand):
            return cand
        warn(f"Not a folder: {cand}")


def choose_project(folder: str) -> Tuple[str, str, List[Feature]]:
    files = list_project_files(folder)
    print("\nWhat do you want to work on?")
    for i, fn in enumerate(files, start=1):
        print(f"  {i:2d}. {fn}")
    print(f"  {len(files) + 1:2d}. New project")

    choice = prompt_int("Enter number: ", minv=1, maxv=len(files) + 1)

    if choice == len(files) + 1:
        print(
            "\nName the project file. A good convention is genename_lastname "
            "(for example chico_Smith)."
        )
        project = prompt_raw("Project name: ")
        path = os.path.join(folder, project_filename(project))
        if os.path.exists(path):
            warn(f"{os.path.basename(path)} already exists; opening it instead.")
            feats = read_tsv(path)
            target = feats[0].target_gene if feats else prompt_raw("What is your target gene? ")
            return target, path, feats
        write_tsv(path, [])
        print(f"  Created {os.path.basename(path)}")
        target = prompt_raw("What is your target gene? ")
        return target, path, []

    path = os.path.join(folder, files[choice - 1])
    try:
        repair_to_utf8(path)
    except Exception:
        pass
    feats = read_tsv(path)
    target = feats[0].target_gene if feats else os.path.splitext(files[choice - 1])[0]
    print(f"  Loaded {files[choice - 1]} ({len(feats)} rows)")
    return target, path, feats


# ----------------------------------------------------------------------
# Slot / position labels
# ----------------------------------------------------------------------
# Slot is stored BIOLOGICALLY: upstream negative, downstream positive, target
# zero, regardless of which strand the target sits on. The value the student
# sees at the prompt is the value written to the TSV.

def parse_position_label(s: str, up: int, down: int) -> int:
    t = re.sub(r"\s+", " ", s.strip().lower())
    if t in {"t", "target", "0"}:
        return 0
    m = re.match(r"^(u|upstream)\s*([0-9]+)$", t)
    if m:
        k = int(m.group(2))
        if 1 <= k <= up:
            return -k
        raise ValueError(f"Upstream index must be between 1 and {up}.")
    m = re.match(r"^(d|downstream)\s*([0-9]+)$", t)
    if m:
        k = int(m.group(2))
        if 1 <= k <= down:
            return k
        raise ValueError(f"Downstream index must be between 1 and {down}.")
    raise ValueError("Expected Target, Upstream1, Downstream1 (and so on).")


def prompt_position_label(up: int, down: int, prompt: Optional[str] = None) -> int:
    up_range = f"Upstream1..Upstream{up}" if up > 0 else "(no upstream positions)"
    down_range = f"Downstream1..Downstream{down}" if down > 0 else "(no downstream positions)"
    hint = prompt or f"Enter Target, {up_range}, or {down_range}: "
    while True:
        s = prompt_raw(hint)
        try:
            return parse_position_label(s, up, down)
        except ValueError as e:
            warn(str(e))


def position_label_from_slot(slot: int) -> str:
    if slot == 0:
        return "Target"
    if slot < 0:
        return f"Upstream{abs(slot)}"
    return f"Downstream{slot}"


# ----------------------------------------------------------------------
# Gene name entry (always manual, with paste assistance)
# ----------------------------------------------------------------------

def _extract_outermost_parens(text: str) -> Optional[str]:
    """
    Contents of the LAST top-level parenthetical group in `text`, or None.
    'Top-level' means the group opening at depth 0, so a gene name that itself
    contains parentheses - Su(var)205, E(spl)m8 - is returned whole rather than
    truncated at the first inner ')'.
    """
    groups: List[str] = []
    depth = 0
    start = -1
    for i, ch in enumerate(text):
        if ch == "(":
            if depth == 0:
                start = i + 1
            depth += 1
        elif ch == ")" and depth > 0:
            depth -= 1
            if depth == 0 and start != -1:
                groups.append(text[start:i])
                start = -1
    for g in reversed(groups):
        if g.strip():
            return g.strip()
    return None


# -PA / -PB (polypeptide) and -RA / -RB (transcript) isoform suffixes.
_ISOFORM_RE = re.compile(r"-(?:P|R)[A-Z]{1,2}$")


def looks_like_header_line(text: str) -> bool:
    """
    True when the text looks like a pasted SPALN / TBLASTN / FASTA header
    rather than a bare gene name. A bare name with parentheses, such as
    Su(var)205, is a single token and is NOT treated as a header.
    """
    t = text.strip()
    if not t:
        return False
    if t.startswith(">"):
        return True
    if re.search(r"(?i)\bquery\s*=", t):
        return True
    if re.search(r"(?i)\b(alignment|proteins?|blast|spaln|subject|score|expect)\b", t):
        return True
    return len(t.split()) >= 3


def extract_gene_name_candidate(text: str) -> str:
    """
    Best guess at the gene name inside a pasted header line.

    Only two signals are trusted enough to pull a name out of a longer line:
    a 'Query=' field (TBLASTN) and a top-level parenthetical (SPALN). Anything
    else is merely tidied - leading '>', a trailing [Species] tag, stray
    punctuation, an isoform suffix - and handed back for the student to accept
    or overrule. Guessing harder than that produces confident nonsense like
    proposing 'protein' as a gene name.
    """
    s = text.strip()
    if not s:
        return ""

    if s.startswith(">"):
        s = s.lstrip(">").strip()
    s = re.sub(r"\s*\[[^\]]*\]\s*$", "", s)          # trailing [Drosophila x]

    m = re.search(r"(?i)\bquery\s*=\s*(\S+)", s)
    if m:
        s = m.group(1)
    elif looks_like_header_line(s):
        inner = _extract_outermost_parens(s)
        if inner:
            s = inner

    s = s.strip().strip(",;:|").strip()
    s = _ISOFORM_RE.sub("", s).strip()
    return s


GENE_NAME_GUIDANCE = (
    "\n  Gene symbol. This is the short official name (chico, Sik3, Su(var)205),\n"
    "  not a descriptive full name. Copy it out of the SPALN or TBLASTN header\n"
    "  and paste it here - that is the surest way to get the spelling and\n"
    "  capitalisation right.\n"
    "  Pasting the whole header line is fine; the symbol will be pulled out of\n"
    "  it for you to confirm.\n"
    "      SpaIn Alignment of D. melanogaster Proteins (chico-PB)   ->  chico\n"
    "      Query= Su(var)205-PA                                     ->  Su(var)205"
)


def prompt_gene_name(role: str, prefill: str = "") -> str:
    """
    Ask for a gene symbol (the short official name, e.g. chico, Sik3).
    Always manual. `role` is a short description used in the prompt
    ("the target gene", "the gene at Upstream2", ...). If `prefill` is
    given, pressing Enter keeps it.
    """
    print(GENE_NAME_GUIDANCE)
    while True:
        if prefill:
            print(f"\n  Current symbol for {role}: {prefill}")
            raw = prompt_raw(
                "  Press Enter to keep it, or paste/type a different symbol: ",
                allow_empty=True,
            )
            if raw == "":
                return prefill
        else:
            raw = prompt_raw(f"\n  Gene symbol for {role}: ")

        if looks_like_info_block(raw):
            warn(
                "That is a gene information block, not a gene symbol. Type or "
                "paste just the symbol here; the block is pasted at the next step."
            )
            continue

        candidate = extract_gene_name_candidate(raw)

        if candidate and candidate != raw:
            print(f"\n  You pasted:  {raw}")
            print(f"  Gene symbol: {candidate}")
            pick = prompt_choice(
                "Which should be recorded as the gene symbol?",
                {
                    "name": f"{candidate}",
                    "aspasted": f"{raw}",
                    "retype": "Neither - let me type it again",
                },
            )
            if pick == "name":
                return candidate
            if pick == "aspasted":
                return raw
            continue

        chosen = candidate or raw
        if prompt_yesno(f'  Record "{chosen}" as the gene symbol? (y/n) '):
            return chosen


# ----------------------------------------------------------------------
# Info-block parsing
# ----------------------------------------------------------------------
# The parser is deliberately newline-agnostic. It locates each labelled field
# anywhere in the text and takes everything up to the next label as its value.
# That means a block pasted as five clean lines and the same block mashed into
# one run-on line parse identically - which is what makes the program usable in
# JupyterLab, whose input box strips newlines out of a paste.

_LABEL_RE = re.compile(
    r"(?i)("
    r"genbank\s+record|flybase\s+record|ncbi\s+record|refseq\s+record|record"
    r"|position|location|coordinates"
    r"|genomic\s+size|size"
    r"|strand|orientation"
    r"|gene\s+symbol|symbol"
    # Trailing lines the NCBI block carries after Gene Symbol. They are not
    # used, but they must be recognised as labels: in a run-on Jupyter paste
    # there is no newline to stop the symbol value, so without these the
    # LOCID came out as "LOC108110920 CDS Start: complete CDS End: complete".
    # Do NOT generalise this to "any word followed by a colon" - the Position
    # value (NC_057930:14571081-14576048) contains one and would be truncated.
    r"|cds\s+start|cds\s+end"
    r")\s*:"
)

_LABEL_CANON = {
    "record": "record",
    "genbank record": "record",
    "flybase record": "record",
    "ncbi record": "record",
    "refseq record": "record",
    "position": "position",
    "location": "position",
    "coordinates": "position",
    "genomic size": "size",
    "size": "size",
    "strand": "strand",
    "orientation": "strand",
    "gene symbol": "symbol",
    "symbol": "symbol",
}

_POSITION_RE = re.compile(
    r"([A-Za-z0-9_.|:\-]*[A-Za-z0-9_.\-])\s*:\s*"
    r"([\d,]+)\s*(?:\.\.|--|-|\u2013|\u2014|to)\s*([\d,]+)"
)


def _canon_label(raw: str) -> str:
    key = re.sub(r"\s+", " ", raw.strip().lower())
    return _LABEL_CANON.get(key, key)


def split_labeled_fields(text: str) -> Dict[str, str]:
    """Map canonical label -> value, tolerant of missing newlines."""
    out: Dict[str, str] = {}
    matches = list(_LABEL_RE.finditer(text))
    for i, m in enumerate(matches):
        canon = _canon_label(m.group(1))
        end = matches[i + 1].start() if i + 1 < len(matches) else len(text)
        value = text[m.end():end]
        # Every one of these values lives on a single line in the browser, so
        # anything after the first newline belongs to an unlabelled line.
        value = value.split("\n")[0].strip()
        if canon not in out or not out[canon]:
            out[canon] = value
    return out


def normalize_strand(value: str) -> str:
    t = value.strip().lower()
    if not t:
        return ""
    if t.startswith("+") or t.split()[0] in {"plus", "forward", "fwd", "sense", "1"}:
        return "+"
    if t.startswith("-") or t.split()[0] in {"minus", "reverse", "rev", "antisense", "complement"}:
        return "-"
    return ""


def parse_position_value(value: str) -> Optional[Tuple[str, int, int]]:
    """
    (scaffold, start, end) from a Position value, or None if unparseable.
    Accepts NC_057930:14571081-14576048, chr2L:10,242,512-10,247,078,
    NC_057930.1 : 14571081..14576048, and trailing junk after the range.
    """
    m = _POSITION_RE.search(value.strip())
    if not m:
        return None
    scaffold = m.group(1).strip().rstrip(":").strip()
    if not scaffold:
        return None
    start = int(m.group(2).replace(",", ""))
    end = int(m.group(3).replace(",", ""))
    if start > end:
        start, end = end, start
    return scaffold, start, end


@dataclass
class ParsedBlock:
    genbank_record: str = ""
    scaffold: str = ""
    gene_start: int = 0
    gene_end: int = 0
    genomic_size: int = 0
    strand: str = ""
    gene_symbol: str = ""


def block_missing(p: ParsedBlock) -> List[str]:
    """
    Which required fields are still absent. Genomic Size is not required - it
    is computed from the coordinates when the browser did not supply it.
    """
    missing: List[str] = []
    if not p.genbank_record:
        missing.append("record")
    if not p.scaffold or not p.gene_start or not p.gene_end:
        missing.append("position")
    if p.strand not in {"+", "-"}:
        missing.append("strand")
    return missing


def parse_gene_block(text: str) -> Tuple[ParsedBlock, List[str]]:
    """
    Parse whatever can be found in `text`. Returns (parsed, missing).
    Never raises: a partial parse is useful, because the caller can ask for
    only the fields that are missing instead of demanding a whole new paste.
    """
    fields = split_labeled_fields(text)
    p = ParsedBlock()

    p.genbank_record = fields.get("record", "").strip()

    pos = parse_position_value(fields.get("position", ""))
    if pos is None and "position" not in fields:
        # Some browsers emit a bare scaffold:range line with no label at all.
        pos = parse_position_value(text)
    if pos is not None:
        p.scaffold, p.gene_start, p.gene_end = pos

    p.strand = normalize_strand(fields.get("strand", ""))

    size_raw = fields.get("size", "")
    m = re.search(r"[\d,]+", size_raw)
    if m:
        p.genomic_size = int(m.group(0).replace(",", ""))
    elif p.gene_start and p.gene_end:
        p.genomic_size = p.gene_end - p.gene_start + 1

    # A gene symbol never contains whitespace, so the first token is the whole
    # symbol. This is the belt to the label list's braces: if NCBI adds a line
    # after Gene Symbol that the label list does not know yet, the symbol still
    # cannot swallow it.
    symbol_raw = fields.get("symbol", "").strip()
    p.gene_symbol = symbol_raw.split()[0] if symbol_raw else ""
    return p, block_missing(p)


def looks_like_info_block(text: str) -> bool:
    """
    True when `text` is a pasted gene information block (FlyBase or NCBI)
    rather than a gene symbol or a SPALN / TBLASTN header line.
    """
    fields = split_labeled_fields(text)
    return "record" in fields or "position" in fields


def block_text_is_complete(text: str) -> bool:
    _p, missing = parse_gene_block(text)
    return not missing


def looks_like_dmel_scaffold(scaffold: str) -> bool:
    return scaffold.lower().startswith("chr")


def looks_like_refseq_scaffold(scaffold: str) -> bool:
    # Modern RefSeq molecules: NC_, NG_, NT_, NW_, NZ_ and X-prefixed records.
    if re.match(r"^[NX][A-Z]_\d", scaffold):
        return True
    # Older GenBank assembly accessions (CP012523, CM000123) and WGS-style
    # accessions (AABU01000001), with an optional .version suffix.
    if re.match(r"^[A-Z]{1,2}\d{5,}(\.\d+)?$", scaffold):
        return True
    if re.match(r"^[A-Z]{4}\d{6,}(\.\d+)?$", scaffold):
        return True
    return False


def scaffold_format_ok(species: str, scaffold: str) -> bool:
    if not scaffold or scaffold == "Manual":
        return True
    if species == DMEL:
        return looks_like_dmel_scaffold(scaffold)
    return looks_like_refseq_scaffold(scaffold)


# Transcript-style accessions. Older NCBI RefSeq blocks from pre-RefSeq
# assemblies report one of these and carry no Gene Symbol line at all.
_TRANSCRIPT_RECORD_RE = re.compile(r"^[NX][MR]_\d")


def looks_like_legacy_refseq(parsed: ParsedBlock) -> bool:
    if parsed.gene_symbol.strip():
        return False
    return bool(_TRANSCRIPT_RECORD_RE.match(parsed.genbank_record.strip()))


# ----------------------------------------------------------------------
# Nesting model
# ----------------------------------------------------------------------

def children_of(feats: List[Feature], parent_uid: str) -> List[Feature]:
    if not parent_uid:
        return []
    kids = [f for f in feats if f.parent_uid == parent_uid]
    kids.sort(key=lambda x: (x.gene_start, x.gene_end, x.gene_name))
    return kids


def parent_of(feats: List[Feature], f: Feature) -> Optional[Feature]:
    if not f.parent_uid.strip():
        return None
    for g in feats:
        if g.feature_uid == f.parent_uid:
            return g
    return None


def descendants_of(feats: List[Feature], f: Feature) -> List[Feature]:
    out: List[Feature] = []
    stack = list(children_of(feats, f.feature_uid))
    while stack:
        cur = stack.pop()
        out.append(cur)
        stack.extend(children_of(feats, cur.feature_uid))
    return out


def compute_utr_nest(child: Feature, parent: Feature) -> str:
    """'y' when the child's coordinate range sits inside the parent's."""
    return "y" if (child.gene_start >= parent.gene_start
                   and child.gene_end <= parent.gene_end) else "n"


def top_level_genes_for_species(feats: List[Feature], species: str) -> List[Feature]:
    genes = [f for f in feats if f.species == species and not f.parent_uid.strip()]
    genes.sort(key=lambda x: (x.slot, x.gene_start, x.gene_end, x.gene_name))
    return genes


def genes_for_species(feats: List[Feature], species: str) -> List[Feature]:
    return [f for f in feats if f.species == species]


def erase_slot(feats: List[Feature], species: str, slot: int) -> None:
    feats[:] = [f for f in feats if not (f.species == species and f.slot == slot)]


def erase_all_upstream(feats: List[Feature], species: str) -> None:
    feats[:] = [f for f in feats if not (f.species == species and f.slot < 0)]


def erase_all_downstream(feats: List[Feature], species: str) -> None:
    feats[:] = [f for f in feats if not (f.species == species and f.slot > 0)]


def get_target_strand(feats: List[Feature], species: str, target: str) -> str:
    """
    The target's genomic strand, used only to decide whether the preview needs
    flipping. The target may be nested inside a container gene, so look it up
    by name first and fall back to whatever occupies slot 0.
    """
    for f in feats:
        if f.species == species and f.slot == 0 and f.gene_name == target \
                and f.strand in {"+", "-"}:
            return f.strand
    for f in feats:
        if f.species == species and f.slot == 0 and not f.parent_uid.strip() \
                and f.strand in {"+", "-"}:
            return f.strand
    return "+"


def arrow(strand: str) -> str:
    if strand == "+":
        return GLYPH["right"]
    if strand == "-":
        return GLYPH["left"]
    return GLYPH["unknown"]


def flip_strand(s: str) -> str:
    return "-" if s == "+" else ("+" if s == "-" else s)


def set_slot_recursive(feats: List[Feature], f: Feature, new_slot: int) -> None:
    """A nested gene lives at its parent's locus, so descendants follow along."""
    f.slot = new_slot
    for d in descendants_of(feats, f):
        d.slot = new_slot


def refresh_nesting_fields(feats: List[Feature], f: Feature) -> None:
    """Recompute nested_inside and UTR_Nest for `f` and everything under it."""
    p = parent_of(feats, f)
    if p is None:
        f.nested_inside = ""
        f.UTR_Nest = ""
        if f.CDS_Nest:
            f.CDS_Nest = ""
    else:
        f.nested_inside = p.gene_name
        f.UTR_Nest = compute_utr_nest(f, p)
    for k in children_of(feats, f.feature_uid):
        refresh_nesting_fields(feats, k)


# ----------------------------------------------------------------------
# Neighborhood index + rendering
# ----------------------------------------------------------------------
# Every gene in the layout, nested genes included, gets a reference label.
# The label is what the student types to correct one specific gene.
#
# Labels are derived from the POSITION, not from entry order:
#
#     U4  U3  U2  U1  T  D1  D2  D3  D4
#
# with nested genes taking a letter after their container (U1a, U1b), and
# deeper nesting another (U1aa). So the label for a gene never changes as
# other genes are entered around it. Numbering by list position, as earlier
# versions did, meant every reference shifted whenever a gene was added
# further out -- which made the layout hard to follow when entering the
# neighborhood in any order other than strictly outside-in.

def _letter(i: int) -> str:
    return chr(ord("a") + i) if i < 26 else f"({i + 1})"


def slot_ref(slot: int) -> str:
    """Short position label: U4..U1, T, D1..D4."""
    if slot == 0:
        return "T"
    return f"U{abs(slot)}" if slot < 0 else f"D{slot}"


def neighborhood_index(
    feats: List[Feature], species: str, up: int, down: int
) -> List[Tuple[str, Feature, int]]:
    """[(label, feature, depth)] in display order, nested genes included."""
    tops = [f for f in top_level_genes_for_species(feats, species)
            if -up <= f.slot <= down]
    entries: List[Tuple[str, Feature, int]] = []

    def walk(f: Feature, label: str, depth: int) -> None:
        entries.append((label, f, depth))
        for i, k in enumerate(children_of(feats, f.feature_uid)):
            walk(k, label + _letter(i), depth + 1)

    # Two top-level genes should never share a slot, but if data from an older
    # file ever does, disambiguate rather than emitting a duplicate reference
    # the student cannot use.
    seen: Dict[str, int] = {}
    for f in tops:
        base = slot_ref(f.slot)
        n = seen.get(base, 0)
        seen[base] = n + 1
        walk(f, base if n == 0 else f"{base}#{n + 1}", 0)
    return entries


def render_neighborhood(
    feats: List[Feature],
    species: str,
    target: str,
    up: int,
    down: int,
    detail: bool = True,
) -> str:
    """
    The layout in the student frame: the target always points right, and every
    other arrow is flipped with it when the target is on the - strand. Flipping
    is purely visual; the stored `strand` is never touched.
    """
    tstrand = get_target_strand(feats, species, target)

    def disp_strand(f: Feature) -> str:
        if f.slot == 0 and f.gene_name == target:
            return "+"
        return f.strand if tstrand == "+" else flip_strand(f.strand)

    entries = neighborhood_index(feats, species, up, down)
    by_slot: Dict[int, List[Tuple[str, Feature, int]]] = {}
    for lab, f, depth in entries:
        by_slot.setdefault(f.slot, []).append((lab, f, depth))

    lines: List[str] = []
    lines.append(f"Layout for {species}   (target: {target})")
    lines.append("")

    for slot in list(range(-up, 0)) + [0] + list(range(1, down + 1)):
        label = position_label_from_slot(slot)
        occupants = by_slot.get(slot, [])
        if not occupants:
            lines.append(f"  {label:<12} --- empty ---")
            continue
        first = True
        for lab, f, depth in occupants:
            indent = "" if depth == 0 else "   " * depth + GLYPH["branch"]
            head = f"  {label:<12}" if first else f"  {'':<12}"
            lines.append(f"{head}[{lab}] {indent}{f.gene_name}{arrow(disp_strand(f))}")
            first = False

    if tstrand == "-":
        lines.append("")
        lines.append("  (the target is on the - strand in the genome; the arrows above are")
        lines.append("   flipped so the picture reads left to right in transcription order)")

    if detail and entries:
        lines.append("")
        lines.append(
            f"  {'ref':<5}{'gene':<18}{'str':<5}{'scaffold':<16}"
            f"{'start':>13}{'end':>13}{'size':>11}  nested in"
        )
        for lab, f, _d in entries:
            p = parent_of(feats, f)
            lines.append(
                f"  {lab:<5}{f.gene_name[:17]:<18}{(f.strand or '?'):<5}"
                f"{f.scaffold[:15]:<16}{fmt_int(f.gene_start):>13}"
                f"{fmt_int(f.gene_end):>13}{fmt_int(f.genomic_size):>11}"
                f"  {p.gene_name if p else ''}"
            )
    return "\n".join(lines)


def show_neighborhood(
    feats: List[Feature], species: str, target: str, up: int, down: int,
    detail: bool = True,
) -> None:
    clear_screen()
    print("\n" + render_neighborhood(feats, species, target, up, down, detail))


# ----------------------------------------------------------------------
# Reading a pasted block
# ----------------------------------------------------------------------

def read_paste_block(
    prompt_text: str,
    complete_check: Optional[Callable[[str], bool]] = None,
) -> str:
    """
    Read pasted text. Returns "" if the student cancels.

    Termination, in order of preference:
      1. The very first delivery already contains a complete block. This is the
         normal case in JupyterLab, where a multi-line paste arrives as one
         run-on line, and it means no extra keypress is needed. Auto-ending is
         restricted to the first delivery on purpose: in a terminal a paste
         arrives one line at a time, and stopping early would leave the unread
         lines to be swallowed by the next prompt.
      2. A blank line.
      3. 'end' / '.' on a line of its own.
      4. 'clip', which reads the clipboard and finishes immediately.
    """
    print(prompt_text)
    print(PASTE_HELP)

    lines: List[str] = []
    iteration = 0
    while True:
        raw = read_line()
        if raw is None:
            raise QuitSignal()
        if is_quit(raw):
            raise QuitSignal()

        stripped = raw.strip()
        low = stripped.lower()

        if low in CANCEL_TOKENS:
            return ""

        if low in CLIP_TOKENS:
            text = clipboard_text()
            if text:
                cleaned = text.replace("\r\n", "\n").replace("\r", "\n")
                lines.extend(cleaned.split("\n"))
                print(f"  Read {len(cleaned)} characters from the clipboard.")
                break
            warn("Could not read the clipboard from here. Paste the text instead.")
            continue

        if low in END_TOKENS and lines:
            break

        if stripped == "":
            if lines:
                break
            note("Waiting for the block. If it is already on your clipboard, "
                 "type 'clip' to read it straight from there - or 'cancel' to "
                 "back out of this gene.")
            continue

        if "\n" in raw or "\r" in raw:
            cleaned = raw.replace("\r\n", "\n").replace("\r", "\n")
            lines.extend(cleaned.split("\n"))
        else:
            lines.append(raw)

        if iteration == 0 and complete_check and complete_check("\n".join(lines)):
            break
        iteration += 1

    return "\n".join(lines).strip("\n")


def fill_missing_block_fields(parsed: ParsedBlock, missing: List[str], species: str) -> bool:
    """
    Ask for only the fields the parser could not find. Returns False if the
    student would rather start the paste again.
    """
    print("\n  The paste was missing: " + ", ".join(missing) + ".")
    what = prompt_choice(
        "How would you like to handle that?",
        {
            "fill": "Type in just the missing field(s)",
            "repaste": "Paste the whole block again",
            "cancel": "Abandon this gene",
        },
    )
    if what == "repaste":
        return False
    if what == "cancel":
        raise ParseError("Entry cancelled.")

    for field in missing:
        if field == "record":
            parsed.genbank_record = prompt_raw(
                "  Record / accession (for example NC_057930.1): "
            )
        elif field == "position":
            while True:
                raw = prompt_raw(
                    "  Position as scaffold:start-end "
                    "(for example NC_057930:14571081-14576048): "
                )
                pos = parse_position_value(raw)
                if pos:
                    parsed.scaffold, parsed.gene_start, parsed.gene_end = pos
                    break
                warn("Could not read a scaffold and a coordinate range from that.")
        elif field == "strand":
            while True:
                s = normalize_strand(prompt_raw("  Strand (+ or -): "))
                if s:
                    parsed.strand = s
                    break
                warn("Enter + or - (plus/minus/forward/reverse also work).")

    if parsed.gene_start and parsed.gene_end and not parsed.genomic_size:
        parsed.genomic_size = parsed.gene_end - parsed.gene_start + 1
    return True


def _validate_parsed(
    parsed: ParsedBlock, species: str, gene_name: str
) -> Optional[Tuple[str, str]]:
    """
    Returns the LOCID for non-Dmel species, or None if the student wants to
    start over. Raises ParseError only on an explicit abort.
    """
    # Size versus range.
    implied = parsed.gene_end - parsed.gene_start + 1
    if parsed.genomic_size != implied:
        warn(
            f"Genomic Size ({fmt_int(parsed.genomic_size)}) does not match the range "
            f"{fmt_int(parsed.gene_start)}-{fmt_int(parsed.gene_end)}, which spans "
            f"{fmt_int(implied)} bases."
        )
        pick = prompt_choice(
            "Which do you want to keep?",
            {
                "computed": f"Use the computed size, {fmt_int(implied)}",
                "pasted": f"Keep the pasted size, {fmt_int(parsed.genomic_size)}",
                "repaste": "Neither - start this paste again",
            },
        )
        if pick == "computed":
            parsed.genomic_size = implied
        elif pick == "repaste":
            return None

    # Scaffold style.
    if not scaffold_format_ok(species, parsed.scaffold):
        expected = ("a FlyBase scaffold such as chr2L or chr3R" if species == DMEL
                    else "a RefSeq or GenBank accession such as NC_057930 or CP012523")
        warn(
            f"Species is {species} but the scaffold ({parsed.scaffold}) does not look "
            f"like {expected}. Did the paste come from the wrong browser?"
        )
        if not prompt_yesno("  Use it anyway? (y/n) "):
            return None

    # Gene symbol handling.
    if species == DMEL:
        symbol = parsed.gene_symbol.strip()
        if not gene_name:
            if symbol:
                return "", symbol
            warn("No Gene Symbol line was found in the block.")
            return "", prompt_gene_name("this gene")
        if symbol and symbol.lower() != gene_name.lower():
            warn(
                f"The block's Gene Symbol is {symbol}, but you typed {gene_name}."
            )
            pick = prompt_choice(
                "Which name should be recorded?",
                {
                    "typed": f"{gene_name} (what you typed)",
                    "symbol": f"{symbol} (from the block)",
                },
            )
            if pick == "symbol":
                # The student chose the block's spelling. Return it so it
                # reaches the row that gets written, not just the summary.
                parsed.gene_symbol = symbol
                return "", symbol
        return "", gene_name

    locid = parsed.gene_symbol.strip()
    if not locid:
        if looks_like_legacy_refseq(parsed):
            note(
                "This looks like an older NCBI RefSeq record (no Gene Symbol / LOC "
                "id). LOCID will be left blank, which is expected here."
            )
        else:
            warn(
                "No Gene Symbol / LOC id was found. Older NCBI assemblies legitimately "
                "omit it, but a truncated paste looks the same."
            )
            if not prompt_yesno("  Leave LOCID blank and carry on? (y/n) "):
                return None
    elif not locid.upper().startswith("LOC"):
        warn(
            f"The Gene Symbol ({locid}) does not start with 'LOC'. For species other "
            "than D. melanogaster it normally does. Was this a FlyBase block?"
        )
        if not prompt_yesno("  Use it anyway? (y/n) "):
            return None
    return locid, gene_name


def summarize_parsed(parsed: ParsedBlock, gene_name: str, locid: str, species: str) -> None:
    print("\n  Parsed:")
    print(f"    Gene name   {gene_name}")
    print(f"    Scaffold    {parsed.scaffold}")
    print(f"    Range       {fmt_int(parsed.gene_start)} - {fmt_int(parsed.gene_end)}"
          f"   (size {fmt_int(parsed.genomic_size)})")
    print(f"    Strand      {parsed.strand}")
    print(f"    Record      {parsed.genbank_record}")
    if species != DMEL:
        print(f"    LOCID       {locid or '(none)'}")


def collect_gene_block(
    species: str, gene_name: str, role: str = ""
) -> Tuple[ParsedBlock, str, str]:
    """
    Paste-and-parse one gene info block. Returns (parsed, gene_name, locid).

    `gene_name` comes back because the student may overrule what they typed in
    favour of the Gene Symbol in the pasted block; the caller needs the name
    that was actually chosen.
    Raises ParseError if the student abandons the gene.

    `role` is a short description of where in the neighborhood this gene sits
    (e.g. "Upstream2"). If provided, it's shown in the paste prompt so the
    student can confirm they're pasting for the right gene.
    """
    # Build a context line for the paste prompt. The source-view names differ
    # by species: FlyBase for D. melanogaster, NCBI Gene for everything else.
    source_view = "FlyBase" if species == DMEL else "NCBI Gene"
    if not gene_name:
        context = f"for {role}" if role else "for this gene"
    elif role:
        context = f"for {gene_name} at {role}"
    else:
        context = f"for {gene_name}"

    while True:
        text = read_paste_block(
            f"\n  Paste the {source_view} gene information block {context} - all\n"
            f"  of it, from the Record line through Gene Symbol. Record IDs\n"
            f"  starting with NM_, XM_, NC_, NR_, or NW_ are all valid.",
            complete_check=block_text_is_complete,
        )
        if not text:
            raise ParseError("Entry cancelled.")

        parsed, missing = parse_gene_block(text)
        if missing:
            if not fill_missing_block_fields(parsed, missing, species):
                continue

        validated = _validate_parsed(parsed, species, gene_name)
        if validated is None:
            continue
        locid, gene_name = validated

        summarize_parsed(parsed, gene_name, locid, species)
        if prompt_yesno("\n  Is this correct? (y/n) "):
            parsed.gene_symbol = gene_name if species == DMEL else locid
            return parsed, gene_name, locid

        pick = prompt_choice(
            "What would you like to do?",
            {
                "fix": "Correct one field here and now",
                "repaste": "Paste the whole block again",
                "cancel": "Abandon this gene",
            },
        )
        if pick == "repaste":
            continue
        if pick == "cancel":
            raise ParseError("Entry cancelled.")
        parsed, locid = edit_parsed_fields(parsed, locid, species, gene_name)
        summarize_parsed(parsed, gene_name, locid, species)
        if prompt_yesno("\n  Is this correct now? (y/n) "):
            parsed.gene_symbol = gene_name if species == DMEL else locid
            return parsed, gene_name, locid


def edit_parsed_fields(
    parsed: ParsedBlock, locid: str, species: str, gene_name: str
) -> Tuple[ParsedBlock, str]:
    """Field-by-field correction of a freshly parsed block, before it is saved."""
    while True:
        print("\n  Which field is wrong?")
        print(f"    1) Scaffold   {parsed.scaffold}")
        print(f"    2) Start      {fmt_int(parsed.gene_start)}")
        print(f"    3) End        {fmt_int(parsed.gene_end)}")
        print(f"    4) Size       {fmt_int(parsed.genomic_size)}")
        print(f"    5) Strand     {parsed.strand}")
        print(f"    6) Record     {parsed.genbank_record}")
        if species != DMEL:
            print(f"    7) LOCID      {locid or '(none)'}")
        print("    d) done - nothing else is wrong")

        choice = prompt_raw("  Enter a number, or d: ").lower()
        if choice == "d":
            return parsed, locid

        if choice == "1":
            parsed.scaffold = prompt_raw(f"  Scaffold [{parsed.scaffold}]: ")
            if not scaffold_format_ok(species, parsed.scaffold):
                warn("That does not look like a scaffold name for this species.")
        elif choice in {"2", "3"}:
            field = "gene_start" if choice == "2" else "gene_end"
            cur = getattr(parsed, field)
            val = prompt_int(f"  {'Start' if choice == '2' else 'End'} [{fmt_int(cur)}]: ", minv=0)
            setattr(parsed, field, val)
            if parsed.gene_start > parsed.gene_end > 0:
                warn("Start is now greater than End; swapping them.")
                parsed.gene_start, parsed.gene_end = parsed.gene_end, parsed.gene_start
            implied = parsed.gene_end - parsed.gene_start + 1
            if implied != parsed.genomic_size:
                note(f"The range now spans {fmt_int(implied)} bases.")
                if prompt_yesno(f"  Set Size to {fmt_int(implied)}? (y/n) "):
                    parsed.genomic_size = implied
        elif choice == "4":
            parsed.genomic_size = prompt_int(f"  Size [{fmt_int(parsed.genomic_size)}]: ", minv=0)
            implied = parsed.gene_end - parsed.gene_start + 1
            if parsed.genomic_size != implied:
                warn(f"Size no longer matches the range, which spans {fmt_int(implied)} bases.")
        elif choice == "5":
            while True:
                s = normalize_strand(prompt_raw(f"  Strand [{parsed.strand}] (+ or -): "))
                if s:
                    parsed.strand = s
                    break
                warn("Enter + or -.")
        elif choice == "6":
            parsed.genbank_record = prompt_raw(f"  Record [{parsed.genbank_record}]: ")
        elif choice == "7" and species != DMEL:
            locid = prompt_raw(f"  LOCID [{locid or 'none'}]: ", allow_empty=True)
        else:
            warn("Not one of the options.")


def collect_gene_block_manual(
    species: str, gene_name: str
) -> Optional[Tuple[ParsedBlock, str, str]]:
    """
    Guided field-by-field entry, for when the browser or Spaln data is missing
    or wrong. Returns (parsed, locid, why_note), None if abandoned, or BACK if
    the student typed 'back' at the first question to leave manual entry.
    'back' steps to the previous field at any point.
    """
    print(
        "\n  Manual entry. Each field is asked for in turn, with an example of the\n"
        "  format. Type 'back' to return to the previous field, or 'Manual' for a\n"
        "  text field whose value is genuinely unavailable.\n"
        "  'back' at the first question leaves manual entry."
    )

    state: Dict = {
        "why": "", "record": "", "scaffold": "",
        "start": None, "end": None, "strand": "", "locid": "",
    }

    def step_why(st):
        v = prompt_raw_back(
            "  Why is this gene being entered by hand? "
            "(no Spaln hit, RefSeq wrong, ...): "
        )
        if v is BACK:
            return BACK
        st["why"] = v
        return True

    def step_record(st):
        v = prompt_field_value("  Record / accession (for example XM_033294685): ")
        if v is BACK:
            return BACK
        st["record"] = v
        return True

    def step_scaffold(st):
        example = "chr2L" if species == DMEL else "NW_024545513"
        v = prompt_field_value(f"  Scaffold (for example {example}): ")
        if v is BACK:
            return BACK
        st["scaffold"] = v
        return True

    def step_start(st):
        v = prompt_int_field("  Gene start coordinate (for example 4720727): ",
                             minv=0, allow_manual=False)
        if v is BACK:
            return BACK
        st["start"] = v
        return True

    def step_end(st):
        while True:
            v = prompt_int_field("  Gene end coordinate (for example 4762202): ",
                                 minv=0, allow_manual=False)
            if v is BACK:
                return BACK
            if st["start"] is not None and v < st["start"]:
                warn(
                    f"End ({fmt_int(v)}) is smaller than Start ({fmt_int(st['start'])}). "
                    "That can be a typo, or a - strand gene read off the browser "
                    "in transcription order."
                )
                if prompt_choice(
                    "Proceed or re-enter?",
                    {"proceed": "Keep it - they will be stored low to high",
                     "reenter": "Type the end coordinate again"},
                ) == "reenter":
                    continue
            st["end"] = v
            return True

    def step_strand(st):
        while True:
            v = prompt_raw_back("  Strand, + or - (for example -): ")
            if v is BACK:
                return BACK
            s = normalize_strand(v)
            if s:
                st["strand"] = s
                return True
            warn("Strand must be + or - (or 'back').")

    def step_locid(st):
        v = prompt_field_value("  LOCID (for example LOC6534321), or 'Manual' if none: ")
        if v is BACK:
            return BACK
        st["locid"] = v
        return True

    steps: List[Step] = [
        ("why", step_why), ("record", step_record), ("scaffold", step_scaffold),
        ("start", step_start), ("end", step_end), ("strand", step_strand),
    ]
    if species != DMEL:
        steps.append(("locid", step_locid))

    outcome = run_steps(steps, state, back_out=True)
    if outcome is BACK:
        return BACK
    if not outcome:
        return None

    lo = min(state["start"], state["end"])
    hi = max(state["start"], state["end"])
    parsed = ParsedBlock(
        genbank_record=state["record"],
        scaffold=state["scaffold"],
        gene_start=lo,
        gene_end=hi,
        genomic_size=hi - lo + 1,
        strand=state["strand"],
        gene_symbol=state["locid"] if species != DMEL else gene_name,
    )
    locid = "" if species == DMEL else state["locid"]

    if not scaffold_format_ok(species, parsed.scaffold):
        warn(f"The scaffold ({parsed.scaffold}) does not match the usual format for {species}.")
        if not prompt_yesno("  Save it as entered? (y = save, n = start these fields again) "):
            return collect_gene_block_manual(species, gene_name)

    print("\n  Manual entry summary:")
    print(f"    Gene name   {gene_name}")
    print(f"    Scaffold    {parsed.scaffold}")
    print(f"    Range       {fmt_int(parsed.gene_start)} - {fmt_int(parsed.gene_end)}"
          f"   (size {fmt_int(parsed.genomic_size)})")
    print(f"    Strand      {parsed.strand}")
    print(f"    Record      {parsed.genbank_record}")
    if species != DMEL:
        print(f"    LOCID       {locid or '(none)'}")
    print(f"    Reason      {state['why']}")

    if not prompt_yesno("\n  Is this correct? (y/n) "):
        return collect_gene_block_manual(species, gene_name)
    return parsed, locid, state["why"]


def collect_gene_data(
    species: str, role: str, prefill_name: str = ""
) -> Tuple[ParsedBlock, str, str, str]:
    """
    Collect one gene: name first (always typed), then the information block.
    Returns (parsed, gene_name, locid, why_note).
    """
    if species == DMEL:
        # FlyBase's block carries the official symbol on its Gene Symbol line,
        # so for D. melanogaster the name is read from the block instead of
        # being asked for first. (For every other species that line holds a
        # LOC id, and the name has to come from the SPALN / TBLASTN header.)
        while True:
            source = prompt_choice(
                f"How would you like to enter {role} for {species}?",
                {
                    "paste": "Paste the FlyBase gene info block (recommended; the "
                             "gene symbol is read from it)",
                    "manual": "Type the fields in by hand",
                },
            )
            if source == "paste":
                parsed, gene_name, locid = collect_gene_block(species, "", role=role)
                return parsed, gene_name, locid, ""
            gene_name = prompt_gene_name(role, prefill=prefill_name)
            result = collect_gene_block_manual(species, gene_name)
            if result is BACK:
                note("Left manual entry.")
                continue
            if result is None:
                raise ParseError("Manual entry was abandoned.")
            parsed, locid, why = result
            return parsed, gene_name, locid, why

    gene_name = prompt_gene_name(role, prefill=prefill_name)

    while True:
        source = prompt_choice(
            f"How would you like to enter the information for {gene_name}?",
            {
                "paste": "Paste the gene info block from the browser (recommended)",
                "manual": "Type the fields in by hand",
            },
        )

        if source == "paste":
            parsed, gene_name, locid = collect_gene_block(species, gene_name, role=role)
            return parsed, gene_name, locid, ""

        result = collect_gene_block_manual(species, gene_name)
        if result is BACK:
            note("Left manual entry.")
            continue
        if result is None:
            raise ParseError("Manual entry was abandoned.")
        parsed, locid, why = result
        return parsed, gene_name, locid, why


# ----------------------------------------------------------------------
# Building features and the four entry flows
# ----------------------------------------------------------------------

def build_feature(
    target_gene: str,
    species: str,
    slot: int,
    parsed: ParsedBlock,
    gene_name: str,
    locid: str,
    notes_label: str,
    why_note: str,
    parent: Optional[Feature],
    target_nesting_gene: str = "n",
) -> Feature:
    if why_note:
        notes = f"{notes_label} | manual: {why_note}" if notes_label else f"manual: {why_note}"
    else:
        notes = notes_label

    f = Feature(
        target_gene=target_gene,
        species=species,
        species_number=phylogenetic_number(species),
        feature_uid=_new_uid(),
        parent_uid=parent.feature_uid if parent else "",
        gene_name=gene_name,
        locid=locid,
        genbank_record=parsed.genbank_record,
        scaffold=parsed.scaffold,
        slot=slot,
        strand=parsed.strand,
        gene_start=parsed.gene_start,
        gene_end=parsed.gene_end,
        genomic_size=parsed.genomic_size,
        nested_inside=parent.gene_name if parent else "",
        CDS_Nest="",
        UTR_Nest="",
        target_nesting_gene=target_nesting_gene,
        notes=notes,
    )
    if parent is not None:
        f.UTR_Nest = compute_utr_nest(f, parent)
    return f


def append_end_of_gene_notes(f: Feature) -> None:
    n = prompt_raw(
        f"\n  Any notes for {f.gene_name}? (press Enter to skip): ",
        allow_empty=True,
    )
    if n:
        f.notes = f"{f.notes} | {n}" if f.notes else n


def prompt_children_loop(
    parent_feature: Feature,
    target_gene: str,
    species: str,
    feats: List[Feature],
    already: int = 0,
) -> None:
    """
    Collect any genes nested inside `parent_feature`, one at a time. `already`
    is the number of genes known to be inside it before this loop starts, so a
    container that was just entered around the target asks about "any other
    genes" rather than implying the target is not there.
    """
    added = already
    while True:
        word = "any other genes" if added else "any genes"
        if not prompt_yesno(f"\n  Are there {word} nested in {parent_feature.gene_name}? (y/n) "):
            break

        # The CDS judgement is captured while the student is still looking at
        # the browser, before their attention moves to the paste.
        print(
            f"\n  Look at the browser now: does the nested gene sit inside an "
            f"intron of {parent_feature.gene_name}'s CDS?"
        )
        print(
            "  (Overlap between the two CDSs is a 'no' - CDS-nesting means the "
            "smaller gene is entirely within an intron of the larger one.)"
        )
        cds_answer = prompt_yesno_cancel(
            "  CDS nested? (y/n, or 'cancel' to back out) "
        )
        if cds_answer is None:
            note("Cancelled. Returning to the previous prompt.")
            continue
        cds_nest = "y" if cds_answer else "n"

        try:
            child = collect_nested_gene(target_gene, species, parent_feature, feats)
            child.CDS_Nest = cds_nest
            added += 1
        except ParseError as e:
            warn(f"{e} Returning to {parent_feature.gene_name}.")
            continue


def collect_target_gene(
    target_gene: str, species: str, slot: int, feats: List[Feature]
) -> Feature:
    """
    The target gene: slot 0, no parent. Its name is prefilled with the project
    target so a student normally just presses Enter, but it stays editable.
    """
    print(f"\n=== Target gene for {species} ===")
    parsed, gene_name, locid, why = collect_gene_data(
        species, "the target gene", prefill_name=target_gene
    )

    if gene_name != target_gene:
        warn(
            f"You named this gene {gene_name}, but the project target is "
            f"{target_gene}. Recording it as {target_gene} keeps the table "
            "consistent across species."
        )
        if not prompt_yesno(f"  Record it as {target_gene}? (y = yes, n = keep {gene_name}) "):
            note(f"Keeping {gene_name}. The project target is still {target_gene}.")
        else:
            gene_name = target_gene

    f = build_feature(
        target_gene, species, slot, parsed, gene_name, locid,
        notes_label="target", why_note=why, parent=None,
    )
    feats.append(f)

    if prompt_yesno(f"\n  Is {f.gene_name} nested inside another gene? (y/n) "):
        print(
            f"\n  Look at the browser now: does {f.gene_name} sit inside an "
            "intron of the larger gene's CDS?"
        )
        print(
            "  (Overlap between the two CDSs is a 'no' - CDS-nesting means "
            f"{f.gene_name} is entirely within an intron of the larger gene.)"
        )
        cds_answer = prompt_yesno_cancel(
            "  CDS nested? (y/n, or 'cancel' to back out) "
        )
        if cds_answer is None:
            note(
                f"Cancelled. {f.gene_name} stays at the top level; you can "
                "add a container from the correction menu later."
            )
        else:
            cds_nest = "y" if cds_answer else "n"
            try:
                collect_container_gene(f, target_gene, species, feats, cds_nest)
            except ParseError as e:
                warn(f"{e} The container gene was not recorded; {f.gene_name} stays "
                     "at the top level.")

    prompt_children_loop(f, target_gene, species, feats)
    append_end_of_gene_notes(f)
    return f


def collect_synteny_gene(
    target_gene: str, species: str, slot: int, feats: List[Feature]
) -> Feature:
    """An upstream or downstream neighbour. Slot is biological."""
    label = position_label_from_slot(slot)
    print(f"\n=== {label} gene for {species} ===")
    parsed, gene_name, locid, why = collect_gene_data(species, f"the gene at {label}")

    f = build_feature(
        target_gene, species, slot, parsed, gene_name, locid,
        notes_label="", why_note=why, parent=None,
    )
    feats.append(f)

    prompt_children_loop(f, target_gene, species, feats)
    append_end_of_gene_notes(f)
    return f


def collect_nested_gene(
    target_gene: str, species: str, parent: Feature, feats: List[Feature]
) -> Feature:
    """
    A gene nested inside `parent`. It inherits the parent's slot, because a
    nested gene sits at the same locus. The caller sets CDS_Nest afterwards.
    """
    print(f"\n=== Gene nested in {parent.gene_name} ({species}) ===")
    parsed, gene_name, locid, why = collect_gene_data(
        species, f"the gene nested in {parent.gene_name}"
    )

    if parsed.scaffold != parent.scaffold and not why:
        warn(
            f"This gene's scaffold ({parsed.scaffold}) does not match "
            f"{parent.gene_name}'s ({parent.scaffold}). A nested gene must sit on "
            "the same scaffold as the gene containing it."
        )
        if not prompt_yesno("  Save it anyway? (y = save, n = discard this gene) "):
            raise ParseError("Scaffold mismatch between the nested gene and its parent.")

    f = build_feature(
        target_gene, species, parent.slot, parsed, gene_name, locid,
        notes_label="nested", why_note=why, parent=parent,
    )
    feats.append(f)

    prompt_children_loop(f, target_gene, species, feats)
    append_end_of_gene_notes(f)
    return f


def collect_container_gene(
    target_feature: Feature,
    target_gene: str,
    species: str,
    feats: List[Feature],
    target_cds_nest_value: str,
) -> Feature:
    """
    The gene that CONTAINS the target. It is stored as a top-level row at the
    target's slot with its own name, and the target is re-parented onto it so
    that both appear in the layout. (In V10 the link was recorded by name only,
    which left the container invisible in the preview.)
    """
    print(f"\n=== Gene containing {target_feature.gene_name} ({species}) ===")
    parsed, gene_name, locid, why = collect_gene_data(
        species, f"the gene that contains {target_feature.gene_name}"
    )

    if parsed.scaffold != target_feature.scaffold and not why:
        warn(
            f"The container's scaffold ({parsed.scaffold}) does not match the "
            f"target's ({target_feature.scaffold}). A gene and the gene containing "
            "it must be on the same scaffold."
        )
        if not prompt_yesno("  Save it anyway? (y = save, n = discard) "):
            raise ParseError("Scaffold mismatch between the container and the target.")

    f = build_feature(
        target_gene, species, target_feature.slot, parsed, gene_name, locid,
        notes_label=f"contains target gene {target_feature.gene_name}",
        why_note=why, parent=None, target_nesting_gene="y",
    )
    feats.append(f)

    target_feature.parent_uid = f.feature_uid
    target_feature.nested_inside = f.gene_name
    target_feature.CDS_Nest = target_cds_nest_value
    target_feature.UTR_Nest = compute_utr_nest(target_feature, f)

    prompt_children_loop(f, target_gene, species, feats, already=1)
    append_end_of_gene_notes(f)
    return f


def collect_gene_at_slot(
    target_gene: str, species: str, slot: int, feats: List[Feature]
) -> Feature:
    """Dispatch to the right entry flow for a top-level position."""
    if slot == 0:
        return collect_target_gene(target_gene, species, slot, feats)
    return collect_synteny_gene(target_gene, species, slot, feats)


# ----------------------------------------------------------------------
# Granular corrections
# ----------------------------------------------------------------------
# Nothing here ever discards more than the student asked it to. Every action
# works on one gene, or one field of one gene, re-validates what it touched,
# recomputes the values that depend on it, and saves straight away.

def pick_gene(
    feats: List[Feature],
    species: str,
    target: str,
    up: int,
    down: int,
    question: str = "Which gene?",
    exclude: Optional[set] = None,
) -> Optional[Feature]:
    """Choose one gene by its reference label. None means the student backed out."""
    exclude = exclude or set()
    entries = [(lab, f, d) for lab, f, d in neighborhood_index(feats, species, up, down)
               if f.feature_uid not in exclude]
    if not entries:
        note("There is no gene available to choose here.")
        return None
    while True:
        raw = prompt_raw(
            f"  {question} Enter its reference (for example "
            f"{entries[0][0]}), or 'cancel': "
        )
        if is_cancel(raw):
            return None
        want = raw.strip().lower()
        for lab, f, _d in entries:
            if lab.lower() == want:
                return f
        warn(f"No gene here has reference '{raw}'. The references are shown in "
             "square brackets in the layout above.")


def rename_gene(feats: List[Feature], f: Feature, target: str, path: str) -> str:
    """
    Rename one gene, keeping the names stored on its children in step. Returns
    the (possibly new) project target name.
    """
    old = f.gene_name
    new = prompt_gene_name(f"this gene (currently {old})", prefill=old)
    if new == old:
        note("Name unchanged.")
        return target

    f.gene_name = new
    for k in children_of(feats, f.feature_uid):
        k.nested_inside = new
    for g in feats:
        if g.notes and old in g.notes and "contains target gene" in g.notes:
            g.notes = g.notes.replace(f"contains target gene {old}",
                                      f"contains target gene {new}")
    note(f"Renamed {old} to {new}.")

    if old == target:
        warn(
            f"{old} is this project's target gene, recorded in the target_gene "
            "column of every row, for every species."
        )
        if prompt_yesno(f"  Rename the project target to {new} everywhere too? (y/n) "):
            for g in feats:
                g.target_gene = new
            note(f"Project target is now {new}.")
            return new
        note(f"The project target is still {old}.")
    return target


def edit_gene_field(
    feats: List[Feature],
    f: Feature,
    target: str,
    species: str,
    up: int,
    down: int,
    path: str,
) -> str:
    """
    Correct individual fields of one gene. Loops so several fields can be fixed
    in one visit. Returns the (possibly new) project target name.
    """
    while True:
        p = parent_of(feats, f)
        print(f"\n  Editing {f.gene_name}  ({position_label_from_slot(f.slot)}"
              f"{', nested in ' + p.gene_name if p else ''})")
        print(f"     1) Gene name       {f.gene_name}")
        print(f"     2) LOCID           {f.locid or '(none)'}")
        print(f"     3) Record          {f.genbank_record or '(none)'}")
        print(f"     4) Scaffold        {f.scaffold or '(none)'}")
        print(f"     5) Start           {fmt_int(f.gene_start)}")
        print(f"     6) End             {fmt_int(f.gene_end)}")
        print(f"     7) Genomic size    {fmt_int(f.genomic_size)}")
        print(f"     8) Strand          {f.strand or '(none)'}")
        print(f"     9) Position        {position_label_from_slot(f.slot)}"
              f"{'  (follows its parent)' if p else ''}")
        print(f"    10) CDS nested?     {f.CDS_Nest or '(n/a)' if p else '(n/a - not nested)'}")
        print(f"    11) Notes           {f.notes or '(none)'}")
        print("     d) done with this gene")

        choice = prompt_raw("  Enter a number, or d: ").strip().lower()

        if choice == "d":
            write_tsv(path, feats)
            return target

        if choice == "1":
            target = rename_gene(feats, f, target, path)

        elif choice == "2":
            if species == DMEL:
                note("D. melanogaster rows do not carry a LOCID; the gene symbol "
                     "is the name.")
            else:
                f.locid = prompt_raw(f"  LOCID [{f.locid or 'none'}]: ", allow_empty=True)

        elif choice == "3":
            f.genbank_record = prompt_raw(f"  Record [{f.genbank_record or 'none'}]: ")

        elif choice == "4":
            new = prompt_raw(f"  Scaffold [{f.scaffold or 'none'}]: ")
            if not scaffold_format_ok(species, new):
                warn(f"That does not look like a scaffold name for {species}.")
                if not prompt_yesno("  Use it anyway? (y/n) "):
                    continue
            f.scaffold = new
            kids = children_of(feats, f.feature_uid)
            mismatched = [k.gene_name for k in kids if k.scaffold != new]
            if mismatched:
                warn("These nested genes are still on a different scaffold: "
                     + ", ".join(mismatched))
            if p and p.scaffold != new:
                warn(f"{f.gene_name} is nested in {p.gene_name}, which is on "
                     f"{p.scaffold}.")

        elif choice in {"5", "6"}:
            field = "gene_start" if choice == "5" else "gene_end"
            word = "Start" if choice == "5" else "End"
            val = prompt_int(f"  {word} [{fmt_int(getattr(f, field))}]: ", minv=0)
            setattr(f, field, val)
            if f.gene_start > f.gene_end:
                warn("Start is now greater than End. Coordinates are stored low to "
                     "high, so they will be swapped; strand is unaffected.")
                f.gene_start, f.gene_end = f.gene_end, f.gene_start
            _reconcile_size(f)
            refresh_nesting_fields(feats, f)

        elif choice == "7":
            f.genomic_size = prompt_int(f"  Genomic size [{fmt_int(f.genomic_size)}]: ", minv=0)
            implied = f.gene_end - f.gene_start + 1
            if f.genomic_size != implied:
                warn(f"The coordinate range spans {fmt_int(implied)} bases, so the "
                     "size no longer matches it.")

        elif choice == "8":
            while True:
                s = normalize_strand(prompt_raw(f"  Strand [{f.strand or 'none'}] (+ or -): "))
                if s:
                    if f.slot == 0 and f.gene_name == target and s != f.strand:
                        warn("This is the target gene. Changing its strand flips the "
                             "direction the whole layout is drawn in.")
                    f.strand = s
                    break
                warn("Enter + or -.")

        elif choice == "9":
            if p:
                note(f"{f.gene_name} is nested in {p.gene_name}, so it sits at "
                     "the same position. Move the parent, or un-nest this gene "
                     "first (the 'nest' option).")
            else:
                move_gene(feats, f, species, up, down)

        elif choice == "10":
            if not p:
                note("This gene is not nested inside anything, so CDS_Nest does "
                     "not apply.")
            else:
                print(f"\n  Does {f.gene_name} sit inside an intron of "
                      f"{p.gene_name}'s CDS?")
                print(
                    "  (Overlap between the two CDSs is a 'no' - CDS-nesting "
                    f"means {f.gene_name} is entirely within an intron of "
                    f"{p.gene_name}.)"
                )
                f.CDS_Nest = "y" if prompt_yesno("  CDS nested? (y/n) ") else "n"

        elif choice == "11":
            print(f"  Current notes: {f.notes or '(none)'}")
            pick = prompt_choice(
                "What would you like to do with the notes?",
                {"append": "Add to them", "replace": "Replace them",
                 "clear": "Delete them", "keep": "Leave them alone"},
            )
            if pick == "append":
                extra = prompt_raw("  Text to add: ")
                f.notes = f"{f.notes} | {extra}" if f.notes else extra
            elif pick == "replace":
                f.notes = prompt_raw("  New notes: ", allow_empty=True)
            elif pick == "clear":
                f.notes = ""
        else:
            warn("Not one of the options.")

        write_tsv(path, feats)


def _reconcile_size(f: Feature) -> None:
    implied = f.gene_end - f.gene_start + 1
    if f.genomic_size == implied:
        return
    note(f"The range now spans {fmt_int(implied)} bases; the stored size is "
         f"{fmt_int(f.genomic_size)}.")
    if prompt_yesno(f"  Update the size to {fmt_int(implied)}? (y/n) "):
        f.genomic_size = implied


def move_gene(feats: List[Feature], f: Feature, species: str, up: int, down: int) -> None:
    """Move a top-level gene to a different position, taking its nested genes along."""
    cur = position_label_from_slot(f.slot)
    new_slot = prompt_position_label(
        up, down,
        prompt=f"  Move {f.gene_name} from {cur} to which position? "
               f"(Target, Upstream1, Downstream1, ...): ",
    )
    if new_slot == f.slot:
        note("That is where it already is.")
        return

    occupants = [g for g in top_level_genes_for_species(feats, species)
                 if g.slot == new_slot and g.feature_uid != f.feature_uid]
    if occupants:
        names = ", ".join(g.gene_name for g in occupants)
        warn(f"{position_label_from_slot(new_slot)} already holds {names}.")
        pick = prompt_choice(
            "What should happen?",
            {
                "swap": f"Swap them - {names} moves to {cur}",
                "both": "Put both at that position",
                "cancel": "Leave everything where it is",
            },
        )
        if pick == "cancel":
            return
        if pick == "swap":
            old_slot = f.slot
            for g in occupants:
                set_slot_recursive(feats, g, old_slot)

    set_slot_recursive(feats, f, new_slot)
    note(f"{f.gene_name} is now at {position_label_from_slot(new_slot)}.")


def renest_gene(
    feats: List[Feature], f: Feature, species: str, target: str, up: int, down: int
) -> None:
    """Change what a gene is nested inside, or lift it back to the top level."""
    p = parent_of(feats, f)
    print(f"\n  {f.gene_name} is currently "
          + (f"nested in {p.gene_name}." if p else "at the top level (not nested)."))

    pick = prompt_choice(
        f"Where does {f.gene_name} belong?",
        {
            "inside": "Nested inside another gene",
            "top": "Not nested inside anything - a top-level gene",
            "cancel": "Leave it as it is",
        },
    )
    if pick == "cancel":
        return

    if pick == "top":
        if p is None:
            note(f"{f.gene_name} is already a top-level gene.")
            return
        f.parent_uid = ""
        f.CDS_Nest = ""
        refresh_nesting_fields(feats, f)
        slot = prompt_position_label(
            up, down, prompt=f"  Which position does {f.gene_name} occupy now? ",
        )
        set_slot_recursive(feats, f, slot)
        note(f"{f.gene_name} is now a top-level gene at "
             f"{position_label_from_slot(slot)}.")
        return

    # A gene cannot be nested inside itself or inside one of its own
    # descendants; that would make a cycle the tree walkers never escape.
    forbidden = {f.feature_uid} | {d.feature_uid for d in descendants_of(feats, f)}
    new_parent = pick_gene(
        feats, species, target, up, down,
        question=f"Which gene contains {f.gene_name}?",
        exclude=forbidden,
    )
    if new_parent is None:
        return

    f.parent_uid = new_parent.feature_uid
    set_slot_recursive(feats, f, new_parent.slot)
    if f.scaffold and new_parent.scaffold and f.scaffold != new_parent.scaffold:
        warn(f"{f.gene_name} is on {f.scaffold} but {new_parent.gene_name} is on "
             f"{new_parent.scaffold}. A nested gene should share its parent's scaffold.")
    print(f"\n  Does {f.gene_name} sit inside an intron of {new_parent.gene_name}'s CDS?")
    print(
        "  (Overlap between the two CDSs is a 'no' - CDS-nesting means "
        f"{f.gene_name} is entirely within an intron of {new_parent.gene_name}.)"
    )
    f.CDS_Nest = "y" if prompt_yesno("  CDS nested? (y/n) ") else "n"
    refresh_nesting_fields(feats, f)
    note(f"{f.gene_name} is now nested in {new_parent.gene_name}, and its UTR_Nest "
         f"has been recomputed as '{f.UTR_Nest}'.")


def delete_gene(feats: List[Feature], f: Feature, species: str) -> None:
    """Delete one gene, keeping or promoting whatever was nested inside it."""
    kids = children_of(feats, f.feature_uid)
    parent = parent_of(feats, f)

    if kids:
        names = ", ".join(k.gene_name for k in kids)
        print(f"\n  {f.gene_name} has {len(kids)} gene(s) nested inside it: {names}")
        pick = prompt_choice(
            "What should happen to them?",
            {
                "promote": (f"Keep them, attached to {parent.gene_name}"
                            if parent else "Keep them as top-level genes"),
                "delete": "Delete them too, and anything nested in them",
                "cancel": "Do not delete anything",
            },
        )
        if pick == "cancel":
            return
        if pick == "delete":
            doomed = {f.feature_uid} | {d.feature_uid for d in descendants_of(feats, f)}
            n = len(doomed)
            feats[:] = [g for g in feats if g.feature_uid not in doomed]
            note(f"Deleted {f.gene_name} and {n - 1} gene(s) nested within it.")
            return
        for k in kids:
            k.parent_uid = parent.feature_uid if parent else ""
            refresh_nesting_fields(feats, k)
    else:
        if not prompt_yesno(f"  Delete {f.gene_name}? (y/n) "):
            return

    feats[:] = [g for g in feats if g.feature_uid != f.feature_uid]
    note(f"Deleted {f.gene_name}.")


def mirror_neighborhood(feats: List[Feature], species: str, target: str) -> None:
    """
    Swap upstream and downstream for one species. This is the fix for a
    neighborhood entered the wrong way round: only the slot signs change, so
    every gene keeps its own coordinates and strand.
    """
    affected = [f for f in feats if f.species == species and f.slot != 0]
    if not affected:
        note("There is nothing either side of the target to mirror.")
        return
    print(f"\n  This flips every gene for {species}: Upstream1 becomes "
          "Downstream1, and so on. Coordinates and strands are untouched.")
    if not prompt_yesno(f"  Mirror the {len(affected)} gene(s) around {target}? (y/n) "):
        return
    for f in affected:
        f.slot = -f.slot
    note("Mirrored.")


def bulk_erase_menu(
    feats: List[Feature], species: str, target: str, path: str
) -> str:
    """The coarse operations, kept behind their own menu. Returns a status string."""
    pick = prompt_choice(
        "Larger changes",
        {
            "upstream": "Erase every upstream gene for this species",
            "downstream": "Erase every downstream gene for this species",
            "mirror": "Mirror the neighborhood (upstream and downstream swap)",
            "species": f"Erase everything for {species}",
            "cancel": "Do nothing",
        },
    )
    if pick == "cancel":
        return "ok"
    if pick == "upstream":
        if prompt_yesno("  Erase every upstream gene? (y/n) "):
            erase_all_upstream(feats, species)
            write_tsv(path, feats)
            note("Upstream positions erased.")
            return "restart_upstream"
        return "ok"
    if pick == "downstream":
        if prompt_yesno("  Erase every downstream gene? (y/n) "):
            erase_all_downstream(feats, species)
            write_tsv(path, feats)
            note("Downstream positions erased.")
            return "restart_downstream"
        return "ok"
    if pick == "mirror":
        mirror_neighborhood(feats, species, target)
        write_tsv(path, feats)
        return "ok"
    warn(f"This deletes every row for {species}, the target gene included.")
    if prompt_yesno(f"  Really erase all of {species}? (y/n) "):
        feats[:] = [f for f in feats if f.species != species]
        write_tsv(path, feats)
        note(f"All {species} rows erased.")
    return "ok"


def correction_menu(
    target: str,
    species: str,
    feats: List[Feature],
    up: int,
    down: int,
    path: str,
) -> Tuple[str, str]:
    """
    The correction loop shown whenever something in the layout is wrong.
    Returns (status, target) where status is 'ok', 'restart_upstream' or
    'restart_downstream'.
    """
    while True:
        show_neighborhood(feats, species, target, up, down)

        pick = prompt_choice(
            "What needs correcting?",
            {
                "done": "Nothing - the layout is right",
                "field": "One field of one gene is wrong",
                "regene": "Re-enter one whole gene",
                "delete": "Delete a gene",
                "move": "Move a gene to a different position",
                "nest": "Change what a gene is nested inside",
                "add": "Add a gene at an empty position",
                "more": "Larger changes (erase a side, mirror, wipe species)",
                "quit": "Save and quit",
            },
        )

        if pick == "done":
            write_tsv(path, feats)
            return "ok", target

        if pick == "quit":
            raise QuitSignal()

        if pick == "add":
            slot = prompt_position_label(up, down)
            existing = [g for g in top_level_genes_for_species(feats, species)
                        if g.slot == slot]
            if existing:
                warn(f"{position_label_from_slot(slot)} already holds "
                     + ", ".join(g.gene_name for g in existing) + ".")
                if not prompt_yesno("  Add another gene there anyway? (y/n) "):
                    continue
            try:
                collect_gene_at_slot(target, species, slot, feats)
            except ParseError as e:
                warn(f"{e} Nothing was added.")
            write_tsv(path, feats)
            continue

        if pick == "more":
            status = bulk_erase_menu(feats, species, target, path)
            if status != "ok":
                return status, target
            continue

        f = pick_gene(feats, species, target, up, down)
        if f is None:
            continue

        if pick == "field":
            target = edit_gene_field(feats, f, target, species, up, down, path)
        elif pick == "delete":
            delete_gene(feats, f, species)
        elif pick == "move":
            if parent_of(feats, f):
                note(f"{f.gene_name} is nested, so it follows its parent's "
                     "position. Use 'nest' to detach it first.")
            else:
                move_gene(feats, f, species, up, down)
        elif pick == "nest":
            renest_gene(feats, f, species, target, up, down)
        elif pick == "regene":
            parent = parent_of(feats, f)
            slot = f.slot
            name = f.gene_name
            kids = children_of(feats, f.feature_uid)
            if kids:
                warn(f"{name} has {len(kids)} nested gene(s). Re-entering it "
                     "keeps them attached.")
            try:
                if parent is None:
                    new = collect_gene_at_slot(target, species, slot, feats)
                else:
                    new = collect_nested_gene(target, species, parent, feats)
            except ParseError as e:
                warn(f"{e} {name} was left as it was.")
                continue
            for k in kids:
                k.parent_uid = new.feature_uid
                refresh_nesting_fields(feats, k)
            feats[:] = [g for g in feats if g.feature_uid != f.feature_uid]
            note(f"{name} was replaced by {new.gene_name}.")

        write_tsv(path, feats)


# ----------------------------------------------------------------------
# Position and species workflow
# ----------------------------------------------------------------------

def collect_position_with_checkpoint(
    target: str,
    species: str,
    slot: int,
    feats: List[Feature],
    up: int,
    down: int,
    path: str,
    replace: bool = True,
) -> Tuple[str, str]:
    """
    Enter the gene at one position, then show the layout and offer corrections.
    Returns (status, target).
    """
    if replace:
        existing = [f for f in feats if f.species == species and f.slot == slot]
        if existing:
            names = ", ".join(f.gene_name for f in existing)
            warn(f"{position_label_from_slot(slot)} already holds {names}.")
            pick = prompt_choice(
                "What would you like to do?",
                {
                    "replace": "Replace what is there",
                    "keep": "Keep it and correct it instead",
                    "cancel": "Do nothing",
                },
            )
            if pick == "cancel":
                return "ok", target
            if pick == "keep":
                return correction_menu(target, species, feats, up, down, path)
            erase_slot(feats, species, slot)

    try:
        collect_gene_at_slot(target, species, slot, feats)
    except ParseError as e:
        warn(f"{e} This position was not saved.")
        if prompt_yesno("  Try this position again? (y/n) "):
            return collect_position_with_checkpoint(
                target, species, slot, feats, up, down, path, replace=False
            )
        return "ok", target

    write_tsv(path, feats)
    show_neighborhood(feats, species, target, up, down)

    if prompt_yesno("\n  Does this all look right? (y/n) "):
        return "ok", target
    return correction_menu(target, species, feats, up, down, path)


def edit_species_session(
    target: str,
    species: str,
    feats: List[Feature],
    window: Dict[str, int],
    path: str,
) -> str:
    """Work on one species until the student is done with it. Returns the target name."""
    for f in feats:
        if f.species == species:
            f.target_gene = target
            if not f.species_number:
                f.species_number = phylogenetic_number(species)
    write_tsv(path, feats)

    while True:
        up, down = window["up"], window["down"]
        show_neighborhood(feats, species, target, up, down)

        pick = prompt_choice(
            f"What would you like to do for {species}?",
            {
                "enter": "Enter the gene at a position",
                "fix": "Correct something in the layout",
                "window": f"Change the window (currently {up} up, {down} down)",
                "done": f"Finished with {species}",
                "quit": "Save and quit",
            },
        )

        if pick == "done":
            return target
        if pick == "quit":
            raise QuitSignal()

        if pick == "window":
            window["up"] = prompt_int("  How many positions upstream? ", minv=0)
            window["down"] = prompt_int("  How many positions downstream? ", minv=0)
            continue

        if pick == "fix":
            status, target = correction_menu(
                target, species, feats, window["up"], window["down"], path
            )
            continue

        slot = prompt_position_label(window["up"], window["down"])
        status, target = collect_position_with_checkpoint(
            target, species, slot, feats, window["up"], window["down"], path
        )


# ----------------------------------------------------------------------
# Main
# ----------------------------------------------------------------------

def choose_species_to_work_on(feats: List[Feature]) -> str:
    present = species_in_table(feats)
    if not present:
        return choose_species()

    print("\nSpecies already in the table:")
    for s in present:
        n = len([f for f in feats if f.species == s])
        print(f"  - {s}  ({n} row{'s' if n != 1 else ''})")

    # 'back' at either species list returns to this first question.
    while True:
        pick = prompt_choice(
            "Which species do you want to work on?",
            {"existing": "One already in the table", "new": "Add a species"},
        )
        missing = [s for s in SPECIES_FULL if s not in set(present)]
        if pick == "existing":
            chosen = choose_from_list(
                "Select species to work on:", present, allow_back=True)
        elif not missing:
            note("Every species is already in the table.")
            chosen = choose_from_list(
                "Select species to work on:", present, allow_back=True)
        else:
            chosen = choose_from_list(
                "Which species would you like to add?", missing, allow_back=True)
        if chosen is not BACK:
            return chosen


def main(folder: Optional[str] = None) -> int:
    print(BANNER)

    path: Optional[str] = None
    feats: Optional[List[Feature]] = None

    try:
        if folder is None:
            folder = choose_folder()
        elif not os.path.isdir(folder):
            print(f"Error: not a folder: {folder}")
            return 1

        target, path, feats = choose_project(folder)

        window = {
            "up": prompt_int(
                "\nHow many positions upstream (not counting nested genes)? ", minv=0
            ),
            "down": prompt_int(
                "How many positions downstream (not counting nested genes)? ", minv=0
            ),
        }

        for f in feats:
            f.target_gene = target
            if not f.feature_uid:
                f.feature_uid = _new_uid()
        write_tsv(path, feats)

        while True:
            species = choose_species_to_work_on(feats)
            target = edit_species_session(target, species, feats, window, path)
            write_tsv(path, feats)
            if not prompt_yesno("\nWork on another species? (y/n) "):
                break

        print("\nDone.")
        print(f"  File:   {os.path.basename(path)}")
        print(f"  Folder: {os.path.dirname(os.path.abspath(path))}")
        print(f"  Rows:   {len(feats)}")
        return 0

    except ProjectError as e:
        print(f"\n{e}")
        return 1

    except (QuitSignal, KeyboardInterrupt):
        if path is not None and feats is not None:
            try:
                write_tsv(path, feats)
                print(f"\nSaved and quit: {os.path.basename(path)} ({len(feats)} rows)")
            except Exception as e:
                print(f"\nWarning: could not save on quit: {e}")
        else:
            print("\nQuit before any project was opened; nothing to save.")
        return 0


# ==== FOOTER BEGIN ====
def clean_locid_column(path: str, dry_run: bool = True) -> int:
    """
    Repair a project file built before the V12.1 LOCID fix.

    Rows whose locid carries trailing text ("LOC108110920 CDS Start: complete
    CDS End: complete") are trimmed to their first token. Nothing else is
    touched. Run with dry_run=True first to see what would change; with
    dry_run=False a copy of the original is kept beside it as
    <name>.pre_locid_fix before the file is rewritten.

        from MicroSynteny_Table_Builder import clean_locid_column
        clean_locid_column("Sik3_Stanga.tsv")                 # report only
        clean_locid_column("Sik3_Stanga.tsv", dry_run=False)  # fix in place

    Returns the number of rows that were (or would be) changed.
    """
    feats = read_tsv(path)
    changed = 0
    for f in feats:
        raw = (f.locid or "").strip()
        fixed = raw.split()[0] if raw else ""
        if fixed != raw:
            changed += 1
            print(f"  {f.species} / {f.gene_name}: {raw!r} -> {fixed!r}")
            if not dry_run:
                f.locid = fixed
    if changed == 0:
        print(f"{os.path.basename(path)}: no locid values need cleaning.")
        return 0
    if dry_run:
        print(f"{os.path.basename(path)}: {changed} row(s) would be changed. "
              f"Re-run with dry_run=False to apply.")
        return changed
    shutil.copy2(path, path + ".pre_locid_fix")
    write_tsv(path, feats)
    print(f"{os.path.basename(path)}: {changed} row(s) cleaned; original kept "
          f"as {os.path.basename(path)}.pre_locid_fix")
    return changed


def run(folder: Optional[str] = None) -> int:
    """
    Start the builder from a notebook cell.

        from MicroSynteny_Table_Builder import run
        run()                      # use the notebook's own folder
        run(folder="C:/data")      # use a specific folder

    Returns 0 on a normal finish, 1 if a project file could not be opened.
    """
    if not _in_ipython():
        print(
            "Note: this is the JupyterLab edition. In a plain terminal, "
            "MicroSynteny_Table_Builder_V11.py is the one to use - it handles "
            "multi-line pastes and clipboard reads better outside a notebook."
        )
    return main(folder)


if __name__ == "__main__":
    run()
# ==== FOOTER END ====
