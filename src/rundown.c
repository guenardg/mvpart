/* **************************************************************************
 *
 *    Actual version maintained by:
 *    Guillaume Guenard <guillaume.guenard@gmail.com>
 *    Department de sciences biologiques,
 *    Universite de Montreal
 *    Montreal, QC, Canada
 *    from 2024-04-13
 *
 *    Previous versions:
 *    rpart by Terry M Therneau and Beth Atkinson <atkinson@mayo.edu>,
 *    R port of rpart by Brian Ripley <ripley@stats.ox.ac.uk>,
 *    Some routines added from vegan by Jari Oksanen <jari.oksanen@oulu.fi>,
 *    Extensions and adaptations of rpart to mvpart by Glenn De'ath
 *    <g.death@aims.gov.au>
 *
 *    This file is part of mvpart
 *
 *    mvpart is free software: you can redistribute it and/or modify
 *    it under the terms of the GNU General Public License as published by
 *    the Free Software Foundation, either version 3 of the License, or
 *    (at your option) any later version.
 *
 *    mvpart is distributed in the hope that it will be useful,
 *    but WITHOUT ANY WARRANTY; without even the implied warranty of
 *    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *    GNU General Public License for more details.
 *
 *    You should have received a copy of the GNU General Public License
 *    along with mvpart. If not, see <https://www.gnu.org/licenses/>.
 *
 * SCCS @(#)rundown.c   1.5 06/03/01
 * Run an observation down the tree, and return the prediction error,
 * for several CP values at once.
 *
 * *************************************************************************/

#include <stdio.h>
#include "rpart.h"
#include "node.h"
#include "rpartproto.h"

void rundown(
    struct node *tree,
    int obs,
    double *cp,
    double *xpred,
    double *xtemp) {
  
  int i;
  struct node *otree;
  /*
   *  Now, repeat the following: for the cp of interest, run down the tree
   *  until I find a node with smaller complexity.  The parent node will
   *  not have collapsed, but this split will have, so this is my
   *  predictor.
   */
  otree = tree;
  for(i = 0; i < rp.num_unique_cp; i++) {
    while(cp[i] < tree->complexity) {
      tree = branch(tree, obs);
      if(tree == 0)
        goto oops;
      otree = tree;
    }
    xpred[i] =  tree->response_est[0];
    xtemp[i] = (*rp_error)(rp.ydata[obs], tree->response_est);
  }
  return;
  
oops:;
  
  if(rp.usesurrogate < 2) { /*must have hit a missing value */
    for(; i < rp.num_unique_cp; i++)
      xpred[i] =  otree->response_est[0];
    xtemp[i] = (*rp_error)(rp.ydata[obs], otree->response_est);
    return;
  }
  /*
   *  I never really expect to get to this code.  It can only happen if
   *  the last cp on my list is smaller than the terminal cp of the
   *  xval tree just built.  This is impossible (I think).  But just in
   *  case I put a message here.
   */
  REprintf("Warning message--see rundown.c\n");
}
