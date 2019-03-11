// This file is part of Bembel, the higher order C++ boundary element library.
// It was written as part of a cooperation of J. Doelz, H. Harbrecht, S. Kurz,
// M. Multerer, S. Schoeps, and F. Wolf at Technische Universtaet Darmstadt,
// Universitaet Basel, and Universita della Svizzera italiana, Lugano. This
// source code is subject to the GNU General Public License version 3 and
// provided WITHOUT ANY WARRANTY, see <http://www.bembel.eu> for further
// information.
#ifndef __BEMBEL_C_TICTOC__
#define __BEMBEL_C_TICTOC__

#include <stdio.h>
#include <stdlib.h>
#include "sys/time.h"

namespace Bembel {

typedef struct {
  double dtime;
  struct timeval time1;
  struct timeval time2;
} tictoc;

inline void tictoc_start(tictoc *t) {
  gettimeofday(&t->time1, NULL);

  return;
}

inline double tictoc_stop(tictoc *t, char text[55]) {
  gettimeofday(&t->time2, NULL);
  t->dtime = t->time2.tv_sec - t->time1.tv_sec +
             1e-6 * (t->time2.tv_usec - t->time1.tv_usec);

  printf("%20.7g secs: %s\n", t->dtime, text);

  return t->dtime;
}
}  // namespace Bembel
#endif
