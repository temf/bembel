#!/bin/bash
cd build
cp ../geo/sphere.dat .
./spline_testsuite.out
./bem_testsuite.out
mkdir test
echo "
-> Running the LaplaceSingle example..." 
./LaplaceSingle.out
echo "
-> Running the HelmholtzSingle example..." 
./HelmholtzSingle.out
echo "
-> Running the MaxwellSingle example..." 
./MaxwellSingle.out
echo "
*   *   *   *   *   *   *   *   *   *   *   *
You may run run_latex.sh to visualize output.
*   *   *   *   *   *   *   *   *   *   *   *
"