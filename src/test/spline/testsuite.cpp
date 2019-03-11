// This file is part of Bembel, the higher order C++ boundary element library.
// It was written as part of a cooperation of J. Doelz, H. Harbrecht, S. Kurz,
// M. Multerer, S. Schoeps, and F. Wolf at Technische Universtaet Darmstadt,
// Universitaet Basel, and Universita della Svizzera italiana, Lugano. This
// source code is subject to the GNU General Public License version 3 and
// provided WITHOUT ANY WARRANTY, see <http://www.bembel.eu> for further
// information.
#include "spltest.h"

int main() {
  const int testnum = 6;
  int done = 0;

  std::srand(time(NULL));
  std::cout

      << "\n==================================================================="
         "===="
         "\n";

  std::cout << "This is the \033[1;4mC++ "
               "Spline-Framework\033[0m, written by Felix "
               "Wolf.\n\033[1;4mCopyright (c), "
               "2018, Technische Universitaet Darmstadt.\033[0m\n\n";
  std::cout
      << "This software was written as a basis for numerical codes of "
         "higher\norder. It relies heavily on template recursion, and thus, "
         "compiler\noptimisation. Please, only use this if compiled with -O2 "
         "or -O3. It\nincludes higher order polynomial bases, as well as "
         "routines to handle\nNURBS geometries, and some routines for "
         "visualisation. Many functions\nhave been tailored for usage in "
         "BEMBEL, see <http://www.bembel.eu/> for\nmore information. This "
         "program is part of free software, software which\nyou can "
         "redistribute and/or modify under the terms of the GNU "
         "General\nPublic License as published by the Free Software "
         "Foundation, either\nversion 3 of the License, or (at your option) "
         "any later version. The\nsoftware is distributed in the hope that it "
         "will be useful, but WITHOUT\nANY WARRANTY; without even the implied "
         "warranty of MERCHANTABILITY or\nFITNESS FOR A PARTICULAR PURPOSE. "
         "See the GNU General Public License\nfor more details.  You should "
         "have received a copy of the GNU General\nPublic License along with "
         "this software.\nIf not, see <http://www.gnu.org/licenses/>.\n";

  std::cout
      << "\n\033[1;4mI am running automated tests with random input.\033[0m\n";
  int failcount = 0;

  failcount += run_test("bernstein.h", ++done, testnum, test_bernstein());
  failcount += run_test("deBoor.h", ++done, testnum, test_deBoor());
  failcount += run_test("knots.h\t", ++done, testnum, test_knots());
  failcount += run_test("localise.h", ++done, testnum, test_localise());
  failcount += run_test("patch.h\t", ++done, testnum, test_patch());
  failcount += run_test("shape.cpp", ++done, testnum, test_import());

  if (failcount == 0) {
    std::cout << "\033[1;32mAll tests of Spline-Framework\t\tOK!\033[0m\n";
    std::cout << "My job is done !\t\t\tExiting...\n";
    std::cout << "============================================================="
                 "=========="
                 "\n";
    return 0;
  }

  std::cout << "\033[1;31m" << failcount
            << " tests of Spline-Framework FAILED!\033[0m\n";
  std::cout << "You should check whats wrong\t\tExiting...\n\n";

  return failcount;
}
