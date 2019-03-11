// This file is part of Bembel, the higher order C++ boundary element library.
// It was written as part of a cooperation of J. Doelz, H. Harbrecht, S. Kurz,
// M. Multerer, S. Schoeps, and F. Wolf at Technische Universtaet Darmstadt,
// Universitaet Basel, and Universita della Svizzera italiana, Lugano. This
// source code is subject to the GNU General Public License version 3 and
// provided WITHOUT ANY WARRANTY, see <http://www.bembel.eu> for further
// information.
#ifndef _testframe_
#define _testframe_

inline int run_test(std::string name, int done, int testnum, int in) {
  std::cout << "-- Test " << done << "/" << testnum << "\t" << name << "\t"
            << std::flush;

  std::cout << "\t"
            << ((in == 0) ? "\033[1;32mOK!\033[0m" : "\033[1;31mFailed!\033[0m")
            << std::endl;
  return in;
}

#endif