#!/bin/sh
# build.sh -- the full build of EDGE paper 3, in the order the pieces depend on each other.
#
#   sh build.sh          rebuild the tables from the deposits, then compile
#   sh build.sh --quick  compile only
#
# The fix_bbl.py step is not optional. wileyNJD-Chicago.bst attaches a one-character year suffix to
# entries it believes need disambiguating, gets it wrong on this bibliography (it printed "1986j"),
# and once past twenty-six such entries emits "{" as the character, which is an unbalanced brace and
# stops LaTeX at \end{thebibliography}. fix_bbl.py removes the suffixes and refuses to write if the
# braces do not balance afterwards.
set -e
cd "$(dirname "$0")"

if [ "$1" != "--quick" ]; then
  echo "== tables from the deposits =="
  Rscript make_tables_paper3.R
fi

echo "== pass 1 =="
pdflatex -interaction=nonstopmode edge3.tex > /dev/null || true   # a stale .bbl from a failed run may error here; bibtex rebuilds it
echo "== bibtex =="
bibtex edge3 > /dev/null || true
python fix_bbl.py edge3.bbl
echo "== pass 2 =="
pdflatex -interaction=nonstopmode edge3.tex > /dev/null
echo "== pass 3 =="
pdflatex -interaction=nonstopmode edge3.tex > /dev/null

echo
echo "undefined references or citations: $(grep -ciE 'undefined' edge3.log)"
echo "LaTeX errors:                      $(grep -cE '^! ' edge3.log)"
pdfinfo edge3.pdf | grep -i '^pages'
