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
 * SCCS @(#)rpcountup.c 1.5 06/06/01
 * Count up the number of nodes and splits in the final result.
 *
 * Gather the counts for myself, add in those of my children, and
 * pass the total back to my parent.
 *
 * *************************************************************************/

#include "rpart.h"
#include "node.h"
#include "rpartproto.h"

void rpcountup(
    struct node *me,
    int *nnode,
    int *nsplit,
    int *ncat) {
  
  int node2, split2;
  int cat2;
  int i, j, k;
  struct split *ss;
  
  if((me->complexity <= rp.alpha) || (me->leftson == 0)) { /*no kids */
    *nnode=1;
    *nsplit=0;
    *ncat =0;
  }
  else {
    i=0;
    j=0;
    k=0;
    for(ss = me->primary; ss != 0; ss = ss->nextsplit) {
      i++;
      if(rp.numcat[ss->var_num] > 0)
        k++;
    }
    for(ss = me->surrogate; ss != 0; ss = ss->nextsplit) {
      j++;
      if(rp.numcat[ss->var_num] > 0)
        k++;
    }
    rpcountup(me->leftson, nnode,  nsplit, ncat);
    rpcountup(me->rightson, &node2, &split2, &cat2);
    *nnode += 1 + node2;
    *nsplit += i + j + split2;
    *ncat += k + cat2;
  }
  
  return;
}
