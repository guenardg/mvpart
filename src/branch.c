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
 * SCCS  @(#)branch.c   1.5 06/06/01 
 * The four routines for anova splitting
 * Walk an observation 'one more split' down the tree.  If there are no
 * more splits, return 0, otherwise return the address of the new node.
 * A return of zero also comes about if surrogates aren't being used, and I
 * hit a missing value.
 *
 * *************************************************************************/

#include "rpart.h"
#include "node.h"
#include "rpartproto.h"

struct node *branch(
    struct node *tree,
    int obs) {
  
  int i, j, dir;
  struct node *me;
  struct split *tsplit;
  double **xdata;
  int **sorts;
  if(tree->leftson == 0)
    return(0);
  
  me = tree;
  xdata = rp.xdata;
  sorts = rp.sorts;
  /*
  ** choose left or right son
  **   this may use lots of surrogates before we're done
  */
  tsplit = me->primary;
  j = tsplit->var_num;
  if(rp.numcat[j] == 0) { /* continuous */
    for(i = 0; i < rp.n; i++) {
      if(sorts[j][i] == obs) {  /* found the match */
        if(xdata[j][i] < tsplit->spoint)
          dir = tsplit->csplit[0];
        else
          dir = -tsplit->csplit[0];
        goto down;
      }
    }
  }
  else {
    dir = (tsplit->csplit)[(int)xdata[j][obs] - 1];
    if(dir != 0)
      goto down;
  }
  
  if(rp.usesurrogate == 0)
    return(0);
  /*
  ** use the surrogates
  */
  for(tsplit = me->surrogate; tsplit != 0; tsplit = tsplit->nextsplit) {
    j = tsplit->var_num;
    if(rp.numcat[j] == 0) {
      for (i = 0; i < rp.n; i++) {
        if(sorts[j][i] == obs) {
            if(xdata[j][i] < tsplit->spoint)
              dir =  tsplit->csplit[0];
            else
              dir = -tsplit->csplit[0];
            goto down;
        }
      }
    }
    else {
      dir = (tsplit->csplit)[(int)xdata[j][obs] - 1];
      if (dir != 0)
        goto down;
    }
  }
  
  if(rp.usesurrogate < 2)
    return(0);
  /*
  ** split it by default
  */
  dir = me->lastsurrogate;

down:
  if(dir==LEFT)
    return(me->leftson);
  else
    return(me->rightson);
}
