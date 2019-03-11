// This file is part of Bembel, the higher order C++ boundary element library.
// It was written as part of a cooperation of J. Doelz, H. Harbrecht, S. Kurz,
// M. Multerer, S. Schoeps, and F. Wolf at Technische Universtaet Darmstadt,
// Universitaet Basel, and Universita della Svizzera italiana, Lugano. This
// source code is subject to the GNU General Public License version 3 and
// provided WITHOUT ANY WARRANTY, see <http://www.bembel.eu> for further
// information.
#ifndef __BEMBEL_C_BEMTEST_
#define __BEMBEL_C_BEMTEST_

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <Eigen/Dense>
#include <Eigen/IterativeLinearSolvers>
#include <Eigen/Sparse>
#include <algorithm>
#include <array>
#include <iostream>
#include <unsupported/Eigen/IterativeSolvers>
#include <vector>
#include "BEMRHS.h"
#include "Discretization.hpp"
#include "Geometry.hpp"
#include "Grids.hpp"
#include "H2_level2.h"
#include "H_print_geometry.h"
#include "HierarchicalMatrix.hpp"
#include "Logger.hpp"
#include "PDEproblem.hpp"
#include "Rhs.hpp"
#include "Spline.hpp"
#include "Stopwatch.hpp"
#include "assert.h"
#include "cluster_tree.h"
#include "constants.h"
#include "discretization.h"
#include "error.h"
#include "gauss_square.h"
#include "glue.h"
#include "gram.h"
#include "hmatrixfactory.h"
#include "hmatrixsettings.h"
#include "init_gridpoints.h"
#include "init_points.h"
#include "meshdata.h"
#include "myvector.h"
#include "pdeproblem.h"
#include "phi.h"
#include "pot.h"
#include "tictoc.h"

#include "spline/testframe.h"

namespace Bembel {
namespace Test {
std::vector<Spl::Patch> mkScreen();
std::vector<Spl::Patch> mkSphere();

int test_glue();
int test_projector();
int test_phi();
int test_phiphi();
int test_divdiv();
int test_divconf();
int test_HMatrixWrappert();
}  // namespace Test
}  // namespace Bembel
#endif
