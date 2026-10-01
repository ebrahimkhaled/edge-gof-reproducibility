# fix_bbl.py -- repair the spurious year suffixes wileyNJD-Chicago.bst writes into the .bbl.
#
# The style appends \natexlab{a}, {b}, ... to entries it believes share an author-and-year label.
# On this bibliography it attaches one to entries that share nothing: the compiled bibliography
# printed "Hampel ... 1986j", "Carroll ... 2006k", "Rousseeuw and Leroy 1987w". The counter is a
# single character, so once the bibliography passes twenty-six such entries it runs off the end of
# the alphabet and emits "{" itself:
#
#     2014{\natexlab{{}}.
#
# which is an unbalanced brace and stops LaTeX with "Missing } inserted" at \end{thebibliography}.
# That is what adding five references to this paper triggered.
#
# This script removes the suffixes. They are wrong when they print and fatal when they overflow, and
# no entry in this bibliography needs disambiguating: no two share an author and a year.
#
#   python fix_bbl.py edge2.bbl
#
# Run it between bibtex and the second pdflatex. It is idempotent and it reports what it changed.
import re
import sys

# two shapes: the well-formed suffix {\natexlab{j}}, and the overflow {\natexlab{{}} , which carries
# one closing brace fewer because the style emitted "{" as the counter's character
# past the 26th suffix the style emits the ASCII characters after "z": "{", "|", "}" and on.
# The braces print as {\natexlab{{}} and {\natexlab{}}} and unbalance the file; "|" prints as a
# stray bar. Any single-character suffix is removed, the two brace forms included.
SUFFIX = re.compile(r"\{\\natexlab\{(?:[^{}]\}\}|\{\}\}|\}\}\})")


def main(path):
    src = open(path, encoding="utf-8", errors="replace").read()
    found = SUFFIX.findall(src)
    out = SUFFIX.sub("", src)
    if out.count("{") != out.count("}"):
        sys.exit("fix_bbl: braces still unbalanced after the repair; not written")
    if out != src:
        open(path, "w", encoding="utf-8", newline="\n").write(out)
    print("fix_bbl: removed %d year suffix%s from %s" % (len(found), "" if len(found) == 1 else "es", path))


if __name__ == "__main__":
    main(sys.argv[1] if len(sys.argv) > 1 else "edge2.bbl")
