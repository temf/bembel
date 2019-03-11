#!/bin/bash
cp LaTeX/*.tex build/test/
cd build/test/
pdflatex LaplaceSingle.tex
pdflatex HelmholtzSingle.tex
pdflatex MaxwellSingle.tex
printf "\n  *   *   *   *   *   *   *   *   *   *   *   *   *   *   *   *\nHave a look at the .pdf's generated in the build/test/ directory!\n  *   *   *   *   *   *   *   *   *   *   *   *   *   *   *   *\n"