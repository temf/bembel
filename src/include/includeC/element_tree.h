// This file is part of Bembel, the higher order C++ boundary element library.
// It was written as part of a cooperation of J. Doelz, H. Harbrecht, S. Kurz,
// M. Multerer, S. Schoeps, and F. Wolf at Technische Universtaet Darmstadt,
// Universitaet Basel, and Universita della Svizzera italiana, Lugano. This
// source code is subject to the GNU General Public License version 3 and
// provided WITHOUT ANY WARRANTY, see <http://www.bembel.eu> for further
// information.
#ifndef __BEMBEL_C_ELEMENT_TREE_
#define __BEMBEL_C_ELEMENT_TREE_

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <algorithm>
#include "constants.h"
#include "myvector.h"
#include "sparse.h"
#include "topology.h"

namespace Bembel {

typedef struct _et_node et_node;
typedef struct _et_root et_root;

struct _et_node {
  vector3 midpoint; /**< Circumcenter                            */

  et_node *father; /**< Father                                  */

  et_node *son[4]; /**< Sons                                    */

  double radius; /**< Circumradius                            */

  int level; /**< Level of the node                       */

  int number; /**< Index of the node on this level         */

  int patch; /**< Parameter domain on which the node lives*/

  int index_s; /**< Index on the parameter domain in s      */

  int index_t; /**< Index on the parameter domain in t      */

  int vertex[4]; /**< Indices of the four vertices            */
  /*
   * direction is ll, lr, ul, ur
   */

  int vertexi[4]; /**< index of the basis function to which    */
  /*
   * the ansatz function belongs
   */

  int edges[4]; /**< Indices of the first node not being a   */
                /*
                 * vertex on the edge, direction is edge
                 */
  /*
   * 0-1, 0-2, 1-3, 2-3, only initialized
   */
  /*
   * if the node is a leaf and if there is
   */
  /*
   * a node on the corresponding edge
   */
};

struct _et_root {
  int nop;
  et_node **patch;
};

int init_element_tree(et_root *, vector3 *, int **, int, int, int);

int free_element_tree(et_root *);

}  // namespace Bembel
#endif
