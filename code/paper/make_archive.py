# make_archive.py -- assembles the reproducibility archive that Section 4 of EDGE paper 2 promises and
# that Statistics in Medicine requires as a Data File ("any supporting computer code or simulations that
# allow readers to institute any new methodology proposed").
#
#   python make_archive.py [--out <dir>]
#
# Read-only with respect to the study: it copies, hashes and lists. Nothing is regenerated here, so the
# archive holds exactly the files the paper's tables were built from. Every declaration's stored sha256
# is re-checked on the way in, and a file whose hash has moved stops the build.
import argparse
import hashlib
import json
import os
import shutil
import sys

# The author's working tree. This script lives in <tree>/paper_EDGE/paper2/, so the tree is two levels up;
# EDGE_SOURCE_ROOT overrides it. (The copy deposited in code/paper/ assembles the archive; it is not needed to
# use the archive.)
ROOT = os.environ.get("EDGE_SOURCE_ROOT") or os.path.abspath(
    os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", ".."))
SIM = os.path.join(ROOT, "simulations")
PAPER = os.path.join(ROOT, "paper_EDGE", "paper2")     # this script and the archive's static documents
MS = os.path.join(ROOT, "paper_EDGE", "paper3")        # the manuscript that is submitted (paper 3, 2026-09-30)
THEORY = os.path.join(ROOT, "paper_EDGE", "theory")

# the blocks the paper reports; the main power battery (1a, 1b, 2-7) is summarised in battery/analysis
BLOCKS = ["8", "8L", "8R", "8S", "9", "9R", "9R2", "9R2n", "9R3", "9b", "9bL", "9c", "9d", "9dX",
          "C1b", "R", "R2", "T",
          "R_p3", "R2_p3",   # paper 3: the cohort blocks rerun on one admission per patient, split by patient
          "RQ",              # paper 3: robust quasi-deviance versions of Stukel's test (declaration f7354147)
          "EXT",             # paper 3: external validation against the usual tests (declaration 754cf309)
          "STREAM",          # paper 3: edge.stream() under drift and repeated looks (declaration 6af98e73)
          "TW", "TW_HH",     # paper 3: what protects a test (declaration 80b5a83f and addenda 8bf292c3, de5ebb60)
          "TW_ext"]          # paper 3: the same on frozen predictions (addendum 3, 22e3db97)
MAIN_BLOCKS = ["0", "1a", "1b", "2", "3", "4", "5", "6", "7"]
SKIP_DIR = {"_cache", "_review", "_test", "_test2", "dryrun", "dryrun_first_attempt",
            "_archive"}   # _archive holds exploratory work that predates the declared study

# Inputs that the scripts read from beside themselves in the author's tree, and where the archive keeps them.
# The path rewriting below sends every read of these names to the archive location.
SIM_INPUTS = {"uis_data.rds": "code/simulations",                # UIS data (Hosmer-Lemeshow), real-data cells
              "_census_power_paper2.csv": "results/analysis"}    # the census cache of _fig13_scenarios.R


# ------------------------------------------------------------------------------------------------------------
# Making the deposited R code run from the archive (referee request, 2026-10-01)
#
# The study's scripts were written against absolute paths on the author's machine. As they are copied, every
# such path is rewritten, mechanically, to a path under the archive root. The root is resolved at run time:
#   1. the environment variable EDGE_ARCHIVE_ROOT, if set;
#   2. else, under Rscript, two folders above the running script (code/simulations/x.R -> the archive root);
#   3. else the working directory.
# The map, author's tree -> archive (the same table is in README.md, "How to run"):
#   simulations/battery/analysis/...      -> results/analysis/...
#   simulations/battery/<anything else>   -> results/blocks/<...>
#   simulations/runI_p3_*.csv             -> results/cohort/
#   simulations/<SIM_INPUTS name>         -> the folder SIM_INPUTS gives
#   simulations/...                       -> code/simulations/...
#   paper_EDGE/paper3/...                 -> manuscript/...
#   paper_EDGE/theory/...                 -> declarations/...
#   the tree's root                       -> code/   (one legacy figure script)
#   other output folders of the author's  -> output/<name>/, created on demand, not part of the deposit
#   session scratch folders               -> tempdir()
#   ebrahim.gof development tree          -> the installed package: load_all(PKG) becomes library(ebrahim.gof);
#                                            PKG itself becomes Sys.getenv("EBRAHIM_GOF_SRC") (a source checkout,
#                                            needed only by the package-development checks pkg280_*, map_pkg_*)
#   the author's pre-package function files (an older projects folder)
#                                         -> code/legacy_not_deposited/  (not deposited; see code/INDEX.md)
# A path that matches none of these stops the build, so nothing machine-specific can slip through.
# ------------------------------------------------------------------------------------------------------------
import re

R_PRELUDE = """\
## ---- archive paths (inserted by make_archive.py; the convention is in README.md, "How to run") -------------
EDGE_ARCHIVE_ROOT <- local({
  r <- Sys.getenv("EDGE_ARCHIVE_ROOT")
  if (!nzchar(r)) {
    f <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))
    r <- if (length(f)) file.path(dirname(normalizePath(f[1], winslash = "/")), "..", "..") else getwd()
  }
  r <- normalizePath(r, winslash = "/", mustWork = FALSE)
  Sys.setenv(EDGE_ARCHIVE_ROOT = r)        # so worker processes started from here resolve the same root
  r
})
edge_path <- function(...) file.path(EDGE_ARCHIVE_ROOT, ...)
edge_battery <- function(...) {           # the author's battery/ folder: analysis/ -> results/analysis, rest -> results/blocks
  if (...length() == 0L) return(edge_path("results", "blocks"))
  p <- file.path(...)
  ifelse(p == "analysis" | startsWith(p, "analysis/"), edge_path("results", p), edge_path("results", "blocks", p))
}
edge_out <- function(...) { d <- edge_path("output", ...); dir.create(d, showWarnings = FALSE, recursive = TRUE); d }
## ---------------------------------------------------------------------------------------------------------------
"""

def _norm(path):
    return path.replace("\\", "/").rstrip("/").lower()


# The machine-specific prefixes are derived from this machine, so no path is typed here.
_HOME = _norm(os.path.expanduser("~"))                       # the author's home folder
_TREE = _norm(ROOT)                                          # the working tree (see ROOT)
_SEP = r"(?:/|\\\\|\\)"                                      # /, an escaped \\ in R source, or a single \
_HOME_RE = r"(?:[A-Za-z]:)?" + "".join(_SEP + re.escape(part) for part in _HOME.split("/")[1:])
_ABS_LITERAL = re.compile(r"""(["'])(""" + _HOME_RE + _SEP + r"""[^"'\n]*)\1""", re.I)
_MACHINE = re.compile(_HOME_RE + r"(?:" + _SEP + r"|\b)", re.I)   # used by the final check
_DEV = _HOME + "/.cursor-tutor/projects"                     # the author's projects folder outside the tree
_PKG_TREES = (_DEV + "/ebrahim frangton", _DEV + "/ebrahim.gof-proj-bagoft")   # ebrahim.gof working copies
PKG_COMMENT = ("  # archive: was load_all() of the author's development tree; "
               "install ebrahim.gof (>= 2.9.0) from CRAN instead")


def _rstr(s):
    return '"' + s + '"'


def _tail_args(tail):
    """'/a/b' -> ', "a/b"' (empty when there is no tail)."""
    tail = tail.strip("/")
    return (", " + _rstr(tail)) if tail else ""


def map_abs_path(raw):
    """Map one absolute path of the author's machine to an R expression under the archive root."""
    p = raw.replace("\\\\", "/").replace("\\", "/")
    while "//" in p:
        p = p.replace("//", "/")
    low = p.lower().rstrip("/")
    p = p.rstrip("/")

    def under(prefix):
        return low == prefix or low.startswith(prefix + "/")

    def tail(prefix):
        return p[len(prefix):]

    sim = _TREE + "/simulations"
    if under(sim + "/battery/analysis"):
        return "edge_path(%s)" % _rstr("results/analysis" + tail(sim + "/battery/analysis"))
    if under(sim + "/battery"):
        t = tail(sim + "/battery").strip("/")
        return "edge_battery(%s)" % (_rstr(t) if t else "")
    if under(sim):
        t = tail(sim).strip("/")
        if re.fullmatch(r"runI_p3_\w+\.csv", t):
            return "edge_path(%s)" % _rstr("results/cohort/" + t)
        if t in SIM_INPUTS:
            return "edge_path(%s)" % _rstr(SIM_INPUTS[t] + "/" + t)
        return "edge_path(%s)" % _rstr("code/simulations" + ("/" + t if t else ""))
    for sub, dest in (("/paper_edge/paper3", "manuscript"), ("/paper_edge/theory", "declarations")):
        if under(_TREE + sub):
            return "edge_path(%s)" % _rstr(dest + tail(_TREE + sub))
    for sub, dest in (("/paper_edge/paper2/figures", "paper2_figures"),
                      ("/paper_edge/submission_statistics_in_medicine", "submission_statistics_in_medicine"),
                      ("/figures", "figures")):
        if under(_TREE + sub):
            return "edge_out(%s%s)" % (_rstr(dest), _tail_args(tail(_TREE + sub)))
    if low == _TREE:
        return 'edge_path("code")'
    m = re.match(re.escape(_HOME) + r"/appdata/local/temp/[^/]+/[^/]+/[^/]+/scratchpad", low)
    if m:
        return "file.path(tempdir()%s)" % _tail_args(p[m.end():])
    if under(_HOME + "/desktop"):
        return "file.path(tempdir()%s)" % _tail_args(tail(_HOME + "/desktop"))
    for dev in _PKG_TREES:
        if under(dev):
            return 'Sys.getenv("EBRAHIM_GOF_SRC")'
    if under(_DEV):
        return 'edge_path("code", "legacy_not_deposited"%s)' % _tail_args(tail(_DEV))
    raise ValueError("no archive mapping for the path " + raw)


_MAKE_CLUSTER = re.compile(
    r'(?m)^([ \t]*)(?:[^#\n]*?;[ \t]*)?([A-Za-z_.][\w.]*)[ \t]*(?:<-|=)[ \t]*(?:tryCatch\([ \t]*)?'
    r'(?:parallel::)?make(?:PSOCK)?[Cc]luster\(')
CLUSTER_EXPORT = ('parallel::clusterExport(%s, c("EDGE_ARCHIVE_ROOT", "edge_path", "edge_battery", "edge_out"), '
                  'envir = environment())  # archive: the path helpers, for the workers')


def _statement_end(text, start):
    """Index just past the newline that ends the R statement beginning at `start` (strings and comments skipped)."""
    depth, i, n = 0, start, len(text)
    while i < n:
        c = text[i]
        if c in "\"'":
            q, i = c, i + 1
            while i < n and text[i] != q:
                i += 2 if text[i] == "\\" else 1
        elif c == "#":
            while i < n and text[i] != "\n":
                i += 1
            continue
        elif c in "([{":
            depth += 1
        elif c in ")]}":
            depth -= 1
        elif c == "\n" and depth <= 0:
            return i + 1
        i += 1
    return n


def _export_helpers_to_clusters(text):
    out, pos, k = [], 0, 0
    for m in _MAKE_CLUSTER.finditer(text):
        if m.start() < pos:
            continue
        end = _statement_end(text, m.start())
        out.append(text[pos:end])
        out.append(m.group(1) + CLUSTER_EXPORT % m.group(2) + "\n")
        pos, k = end, k + 1
    out.append(text[pos:])
    return "".join(out), k


def portable_r(text, name):
    """Rewrite one R script so that it runs from the archive. Returns (text, number of rewrites)."""
    n = 0
    pkg_was_abs = False

    # A. every absolute path literal of the author's machine
    def sub_abs(m):
        nonlocal n, pkg_was_abs
        expr = map_abs_path(m.group(2))
        if expr.startswith("Sys.getenv(\"EBRAHIM_GOF_SRC"):
            pkg_was_abs = True
        n += 1
        return expr
    text = _ABS_LITERAL.sub(sub_abs, text)

    # the variables that now hold code/simulations (SIMDIR, SIM, BT_SIMDIR, sim_dir, here, ...)
    simvars = set(re.findall(r'(?:^|[;{])\s*([A-Za-z_.][\w.]*)\s*(?:<-|=)\s*edge_path\("code/simulations"\)',
                             text, flags=re.M))

    # B. file.path(<simulations dir>, "battery", ...) -> edge_battery(...)
    def sub_bat(m):
        nonlocal n
        if m.group(1) not in simvars:
            return m.group(0)
        n += 1
        return "edge_battery()" if m.group(2) == ")" else "edge_battery("
    text = re.sub(r'file\.path\(\s*([A-Za-z_.][\w.]*)\s*,\s*"battery/?"\s*(\)|,\s*)', sub_bat, text)

    # C. a variable bound to the battery folder, then extended with file.path(VAR, ...)
    batvars = re.findall(r'(?:^|[;{])\s*([A-Za-z_.][\w.]*)\s*(?:<-|=)\s*edge_battery\(\)', text, flags=re.M)
    for v in sorted(set(batvars)):
        rebinds = re.findall(r'(?:^|[;{(,])\s*' + re.escape(v) + r'\s*(?:<-|=)\s*([^\n;]*)', text, flags=re.M)
        if any(not r.strip().startswith("edge_battery()") for r in rebinds):
            continue                                   # rebound elsewhere in the file: leave it alone
        pat = r'file\.path\(\s*' + re.escape(v) + r'\s*(\)|,\s*)'
        text, k = re.subn(pat, lambda m: "edge_battery()" if m.group(1) == ")" else "edge_battery(", text)
        n += k

    # D. inputs and outputs that the archive keeps elsewhere than beside the scripts
    def input_expr(fn):
        dest = "results/cohort" if fn.startswith("runI_p3_") else SIM_INPUTS[fn]
        return "edge_path(%s)" % _rstr(dest + "/" + fn)

    def sub_input_fp(m):
        nonlocal n
        if m.group(1) not in simvars:
            return m.group(0)
        n += 1
        return input_expr(m.group(2))

    def sub_input_bare(m):                 # a bare name, read or written in the working directory
        nonlocal n
        n += 1
        return input_expr(m.group(1))
    names = r"(runI_p3_\w+\.csv|" + "|".join(re.escape(k) for k in SIM_INPUTS if k.endswith(".csv")) + ")"
    text = re.sub(r'file\.path\(\s*([A-Za-z_.][\w.]*)\s*,\s*"' + names + r'"\s*\)', sub_input_fp, text)
    text = re.sub(r'"' + names + r'"', sub_input_bare, text)

    # E. the package: the development tree is replaced by the CRAN release
    if pkg_was_abs:
        lines = text.split("\n")
        for i, line in enumerate(lines):
            new = re.sub(r'(?:devtools|pkgload)::load_all\(\s*PKG\b[^()]*\)', "library(ebrahim.gof)", line)
            if new != line:
                lines[i] = new + PKG_COMMENT
                n += 1
        text = "\n".join(lines)

    if n == 0:
        return text, 0

    # F. worker processes: a function shipped to a PSOCK worker looks its free variables up in the worker's
    #    global environment, so the path helpers are exported right after each cluster is made
    text, k = _export_helpers_to_clusters(text)
    n += k

    # the prelude goes after the script's header comment
    lines = text.split("\n")
    k = 0
    while k < len(lines) and lines[k].startswith("#"):
        k += 1
    text = "\n".join(lines[:k] + R_PRELUDE.rstrip("\n").split("\n") + lines[k:])
    return text, n


# Text corrections made in the deposited copy only (the source file needs the same fix by hand).
TEXT_FIXES = {
    # bt_hlw() in _battery_tests.R is the Hosmer-Lemeshow statistic on G equal-width risk intervals,
    # not a weighted decile test; the comment named the wrong test.
    "analyse_census_extended.R": [("the Hosmer-Hjort weighted decile test (HL_w)",
                                   "the Hosmer-Lemeshow test on equal-width risk groups (HL_w)")],
    "analyse_corruption_partitions.R": [("Hosmer-Hjort weighted deciles HL_w",
                                         "Hosmer-Lemeshow on equal-width risk groups HL_w")],
    # so that a reader can write the tables to a folder of their own and compare them with manuscript/tables
    "make_tables_paper3.R": [('OUT <- edge_path("manuscript/tables")',
                              'OUT <- Sys.getenv("EDGE_TABLES_OUT", edge_path("manuscript/tables"))'
                              '   # archive: set EDGE_TABLES_OUT to write elsewhere')],
}

REWRITE_LOG = []      # (archive path, rewrites) for the build report


def transform_code(src, dst):
    """Copy hook for code files: rewrite R paths, apply text fixes, make the assembler portable."""
    fn = os.path.basename(src)
    if not fn.lower().endswith(".r"):
        return False
    with open(src, "r", encoding="utf-8", errors="surrogateescape", newline="") as f:
        text = f.read()
    text, n = portable_r(text, fn)
    for old, new in TEXT_FIXES.get(fn, []):
        if old in text:
            text = text.replace(old, new)
            n += 1
    os.makedirs(os.path.dirname(dst), exist_ok=True)
    with open(dst, "w", encoding="utf-8", errors="surrogateescape", newline="") as f:
        f.write(text)
    shutil.copystat(src, dst)
    if n:
        REWRITE_LOG.append((dst, n))
    return True


def rmtree(path):
    """Windows keeps a handle on a directory for a moment after it is listed; retry instead of dying."""
    import stat
    import time

    def onexc(func, p, exc):
        try:
            os.chmod(p, stat.S_IWRITE)
        except OSError:
            pass
        time.sleep(0.2)
        func(p)

    shutil.rmtree(path, onexc=onexc)


def sha256(path):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def copy_into(src, dst, pred=None, transform=None, recurse=True):
    """Copy a file or a directory tree, returning the number of files and their total bytes.

    transform(src, dst) may write dst itself and return True; otherwise the file is copied unchanged."""
    n = b = 0
    if os.path.isfile(src):
        os.makedirs(os.path.dirname(dst), exist_ok=True)
        if not (transform and transform(src, dst)):
            shutil.copy2(src, dst)
        return 1, os.path.getsize(dst)
    for dirpath, dirnames, filenames in os.walk(src):
        dirnames[:] = [d for d in dirnames if d not in SKIP_DIR] if recurse else []
        for fn in filenames:
            s = os.path.join(dirpath, fn)
            if pred and not pred(s):
                continue
            d = os.path.join(dst, os.path.relpath(s, src))
            os.makedirs(os.path.dirname(d), exist_ok=True)
            if not (transform and transform(s, d)):
                shutil.copy2(s, d)
            n += 1
            b += os.path.getsize(d)
    return n, b


def check_declarations():
    """Every PREDECLARATION with a stored hash must still hash to it.

    Two sidecar conventions are in use and both are accepted. The plain one holds a single hash. The
    older one is a LOG: the document was revised more than once before its block ran, and each revision
    was recorded with its hash, its time and its reason, newest last. The live document must match the
    LAST hash in such a log; matching an earlier one would mean it had been reverted, which is reported
    separately rather than passed.
    """
    import re
    ok, moved, reverted, unhashed = [], [], [], []
    for fn in sorted(os.listdir(THEORY)):
        if not fn.startswith("PREDECLARATION") or not fn.endswith(".md"):
            continue
        p = os.path.join(THEORY, fn)
        side = next((c for c in (p + ".sha256", p[:-3] + ".sha256") if os.path.exists(c)), None)
        if side is None:
            unhashed.append(fn)
            continue
        logged = re.findall(r"\b[0-9a-f]{64}\b", open(side, encoding="utf-8", errors="replace").read())
        h = sha256(p)
        if not logged:
            unhashed.append(fn)
        elif h == logged[-1]:
            ok.append(fn)
        elif h in logged:
            reverted.append(fn)
        else:
            moved.append(fn)
    return ok, moved, reverted, unhashed


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", default=os.path.join(ROOT, "paper_EDGE", "paper2_archive"))
    a = ap.parse_args()
    out = a.out

    ok, moved, reverted, unhashed = check_declarations()
    print("declarations: %d hashes verified, %d moved, %d reverted, %d without a stored hash"
          % (len(ok), len(moved), len(reverted), len(unhashed)))
    for fn in moved:
        print("  MOVED:", fn)
    for fn in reverted:
        print("  REVERTED to an earlier logged hash:", fn)
    for fn in unhashed:
        print("  no stored hash:", fn)
    if moved or reverted:
        sys.exit("a declaration's content no longer matches its frozen hash; archive not built")

    # Rebuild the tree in place, but never delete .git: the archive is a git repository whose history
    # is the record of what was deposited when, and an earlier version of this script removed it.
    if os.path.exists(out):
        for name in os.listdir(out):
            if name == ".git":
                continue
            p = os.path.join(out, name)
            if os.path.isdir(p):
                rmtree(p)
            else:
                os.remove(p)
    else:
        os.makedirs(out)

    tally = {}

    # 1. the declarations and their hashes
    tally["declarations"] = copy_into(
        THEORY, os.path.join(out, "declarations"),
        pred=lambda p: os.path.basename(p).startswith(("PREDECLARATION", "PILOT_B_PREDECLARATION")))

    # 2. the code: every R file of the study, plus the paper's table and figure builders
    #    (R files are rewritten on the way in so that they run from the archive root; see portable_r)
    tally["code/simulations"] = copy_into(
        SIM, os.path.join(out, "code", "simulations"),
        pred=lambda p: p.lower().endswith(".r") and os.sep + "battery" + os.sep not in p
        and not re.match(r"map_c[a-z]+_", os.path.basename(p)),   # exploratory cost maps, not used by the paper
        transform=transform_code)
    tally["code/paper"] = copy_into(
        MS, os.path.join(out, "code", "paper"),
        pred=lambda p: p.lower().endswith((".r", ".py", ".sh")) and os.sep + "tables" + os.sep not in p,
        transform=transform_code)
    tally["code/paper/make_archive.py"] = copy_into(
        os.path.join(PAPER, "make_archive.py"), os.path.join(out, "code", "paper", "make_archive.py"))
    # the data inputs the scripts read from beside themselves
    for fn, dest in SIM_INPUTS.items():
        tally[dest + "/" + fn] = copy_into(os.path.join(SIM, fn), os.path.join(out, *dest.split("/"), fn))

    # 3. the results: the declared analysis files, and the per-replicate deposits of the reported blocks
    tally["results/analysis"] = copy_into(
        os.path.join(SIM, "battery", "analysis"), os.path.join(out, "results", "analysis"))
    for blk in BLOCKS:
        src = os.path.join(SIM, "battery", blk)
        if os.path.isdir(src):
            tally["results/blocks/" + blk] = copy_into(src, os.path.join(out, "results", "blocks", blk))
    # the main power battery (blocks 0-7): its per-replicate files (~500 MB) stay out, but its block summaries,
    # progress logs and identity checks go in, so the analyses that read the summaries run from the archive
    for blk in MAIN_BLOCKS:
        src = os.path.join(SIM, "battery", blk)
        if os.path.isdir(src):
            tally["results/blocks/" + blk + " (summaries)"] = copy_into(
                src, os.path.join(out, "results", "blocks", blk), recurse=False,
                pred=lambda p: "_pvalues" not in os.path.basename(p))
    # the battery's design files (cell table, sample-size plan, launch plan, seed checks, weighting rule),
    # which the runners and analyses read from the top of the battery folder
    tally["results/blocks (design files)"] = copy_into(
        os.path.join(SIM, "battery"), os.path.join(out, "results", "blocks"), recurse=False,
        pred=lambda p: p.endswith((".csv", ".txt")))

    # 3b. the cohort analysis of paper 3's Section 7: its outputs sit beside the scripts, not in a block
    tally["results/cohort"] = copy_into(
        SIM, os.path.join(out, "results", "cohort"),
        pred=lambda p: os.path.basename(p).startswith("runI_p3_") and p.endswith(".csv")
        and os.path.dirname(p) == SIM)

    # 4. the manuscript, so the archive can rebuild every table and figure it reports
    for sub in ("sections", "tables", "figures"):
        tally["manuscript/" + sub] = copy_into(
            os.path.join(MS, sub), os.path.join(out, "manuscript", sub))
    for fn in ("edge3.tex", "refs.bib", "edge3.bbl", "edge3.pdf", "supplement.tex", "supplement.pdf"):
        s = os.path.join(MS, fn)
        if os.path.exists(s):
            tally["manuscript/" + fn] = copy_into(s, os.path.join(out, "manuscript", fn))

    # 5. the archive's own documents, kept under version control beside the paper
    static = os.path.join(PAPER, "archive_static")
    for fn in sorted(os.listdir(static)):
        tally["(root) " + fn] = copy_into(os.path.join(static, fn), os.path.join(out, fn))

    # 5b. no deposited code may still point into the author's machine
    left = []
    machine = _MACHINE
    for dirpath, dirnames, filenames in os.walk(os.path.join(out, "code")):
        for fn in filenames:
            if fn.lower().endswith((".r", ".py", ".sh")):
                p = os.path.join(dirpath, fn)
                for i, line in enumerate(open(p, encoding="utf-8", errors="replace"), 1):
                    if machine.search(line):
                        left.append("%s:%d" % (os.path.relpath(p, out), i))
    print("code: %d R files rewritten (%d path rewrites); %d lines still name the author's machine"
          % (len(REWRITE_LOG), sum(k for _, k in REWRITE_LOG), len(left)))
    for x in left:
        print("  STILL ABSOLUTE:", x)
    if left:
        sys.exit("absolute paths remain in the deposited code; archive not finished")

    # 6. the manifest: sha256 of everything deposited
    lines, total, count = [], 0, 0
    for dirpath, dirnames, filenames in os.walk(out):
        dirnames[:] = [d for d in dirnames if d != ".git"]   # the repository is not part of the deposit
        for fn in sorted(filenames):
            p = os.path.join(dirpath, fn)
            rel = os.path.relpath(p, out).replace(os.sep, "/")
            lines.append("%s  %s" % (sha256(p), rel))
            total += os.path.getsize(p)
            count += 1
    with open(os.path.join(out, "MANIFEST.sha256"), "w", encoding="utf-8", newline="\n") as f:
        f.write("\n".join(sorted(lines)) + "\n")

    print("\narchive: %s" % out)
    for k in sorted(tally):
        n, b = tally[k]
        print("  %-28s %5d files  %8.1f MB" % (k, n, b / 1e6))
    print("  %-28s %5d files  %8.1f MB" % ("TOTAL (before manifest)", count, total / 1e6))


if __name__ == "__main__":
    main()
