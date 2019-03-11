// This file is part of Bembel, the higher order C++ boundary element library.
// It was written as part of a cooperation of J. Doelz, H. Harbrecht, S. Kurz,
// M. Multerer, S. Schoeps, and F. Wolf at Technische Universtaet Darmstadt,
// Universitaet Basel, and Universita della Svizzera italiana, Lugano. This
// source code is subject to the GNU General Public License version 3 and
// provided WITHOUT ANY WARRANTY, see <http://www.bembel.eu> for further
// information.
#include "bemtest.h"

using namespace Bembel;
using namespace Bembel::Test;

/**
 *  \brief 				 This is an automated testsuite for the
 * bem library.
 */
int main() {
  int done = 0;
  int testnum = 7;

  std::cout

      << "\n==================================================================="
         "===="
         "\n"
         "This is the testsuite for \033[1;4mBembel\033[0m, a boundary element "
         "library written by\nJuergen Doelz, Helmut Harbrecht, Michael "
         "Multerer "
         "and Felix Wolf.\n"
         "\033[1;4mCopyright (c), 2018, Technische Universitaet "
         "Darmstadt.\033[0m\n\n";

  std::cout
      << "This software is a basis for boundary element methods of "
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

  std::cout << "\n\033[1;4mI am running automtated tests.\033[0m\n";
  int failcount = 0;
  failcount += run_test("glue.cpp", ++done, testnum, test_glue());
  failcount += run_test("projector.cpp", ++done, testnum, test_projector());
  failcount += run_test("phi()\t", ++done, testnum, test_phi());
  failcount += run_test("phiphi()", ++done, testnum, test_phiphi());
  failcount += run_test("divdiv()", ++done, testnum, test_divdiv());
  failcount += run_test("div-conforming", ++done, testnum, test_divconf());
  failcount +=
      run_test("c++ wrappers", ++done, testnum, test_HMatrixWrappert());

  assert(failcount == 0 && "A test failed, see console output!");

  if (failcount == 0) {
    std::cout << "\033[1;32mAll tests of Bembel\t\t\tOK!\033[0m\n";
    std::cout << "My job is done !\t\t\tExiting...\n";
    std::cout << "============================================================="
                 "=========="
              << std::endl;
    return 0;
  }
}
