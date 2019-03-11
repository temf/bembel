// This file is part of Bembel, the higher order C++ boundary element library.
// It was written as part of a cooperation of J. Doelz, H. Harbrecht, S. Kurz,
// M. Multerer, S. Schoeps, and F. Wolf at Technische Universtaet Darmstadt,
// Universitaet Basel, and Universita della Svizzera italiana, Lugano. This
// source code is subject to the GNU General Public License version 3 and
// provided WITHOUT ANY WARRANTY, see <http://www.bembel.eu> for further
// information.
#ifndef __C_BEMBEL_SPLTEST_
#define __C_BEMBEL_SPLTEST_

#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <Eigen/Dense>
#include <iostream>
#include <random>
#include <string>
#include <vector>
#include "math.h"
#include "spline/basis.h"
#include "spline/deBoor.h"
#include "spline/helper.h"
#include "spline/knots.h"
#include "spline/localise.h"
#include "spline/patch.h"

#include "testframe.h"

int test_bernstein();
int test_deBoor();
int test_knots();
int test_localise();
int test_patch();
int test_import();

#endif