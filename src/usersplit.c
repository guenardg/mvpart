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
 * SCCS @(#)usersplit.c 1.4 06/06/01
 * These routines interface via the rpart_callback routine to
 * provide for user-written split functions.
 *
 * *************************************************************************/

#include <stdio.h>
#include "rpart.h"
#include "rpartS.h"
#include "node.h"
#include "rpartproto.h"

static int n_return;     /* number of return values from the eval fcn */
static double *uscratch; /* variously used scratch vector */

int usersplit_init(
    int n,
    double *y[],
    int maxcat,
    char **error,
    double *parm,
    int *size,
    int who,
    double *wt) {
  
  if(who == 1) {
    /* If who==0 we are being called internally via xval, and don't
    **   need to rerun the initialization.
    ** Call-back to the C code to get the number of columns for y and
    **   the length of the return vector
    **  the scratch vector needed is of length max(2n, nreturn +1)
    */
    rpart_callback0(&n_return);

    if((n_return + 1) > 2*n)
      uscratch = (double*) ALLOC(n_return + 1, sizeof(double));
    else
      uscratch = (double*) ALLOC(2*n, sizeof(double));
  }
  *size = n_return;
  return(0);
}
/*
** The user evaluation function
*/
void usersplit_eval(
    int n,
    double *y[],
    double *value,
    double *risk,
    double *wt) {
  
  int i;
  
  rpart_callback1(n, y, wt, uscratch);
  *risk = uscratch[0];
  for (i = 0; i < n_return; i++)
    value[i] = uscratch[i + 1];
}
/*
** Call the user-supplied splitting function.
*/
void usersplit(
    int n,
    double *y[],
    double *x,
    int nclass,
    int edge,
    double *improve,
    double *split,
    int *csplit,
    double myrisk,
    double *wt) {
  
  int i, j, k;
  int m;
  int left_n,  right_n;
  int where = 0;
  double best;
  double *dscratch;
  double ftemp;
  /*
  ** If it's categorical, and all are tied, don't bother to callback.
  **    (Completely tied continuous is caught earlier than this).
  ** (This isn't common, but callbacks are expensive).
  */
  if(nclass > 0) {
    ftemp = x[0];
    for(i = 1; i < n; i++)
      if (x[i] != ftemp)
        break;
    if(i == n) {
      *improve =0.0;
      return;
    }
  }
  
  /*
  ** get the vector of "goodness of split"
  **  on return uscratch contains the goodness for each split
  **  followed by the 'direction'
  */
  rpart_callback2(n, nclass, y, wt, x, uscratch);
  
  if(nclass == 0) {
    /*
    **  Find the split point that has the best goodness, subject
    **   to the edge criteria, and tied x's
    **  Remember, uscratch[0] contains the goodnes for x[0] left,
    **   and all others right, so has n-1 real elements.
    **  The 'direction' vector is returned pasted onto the end of
    **   uscratch.
    */
    dscratch = uscratch + n - 1;
    best = 0;
    
    for(i = edge - 1; i < n - edge; i++)
      if((x[i] < x[i + 1]) && (uscratch[i] > best)) {
        best = uscratch[i];
        where = i;
      }
    
    if(best > 0) {   /* found something */
      csplit[0] = dscratch[where];
      *split = (x[where] + x[where+1]) /2;
    }
  }
  else {
    /*
    **  Categorical -- somewhat more work to be done here to
    **   guarrantee the edge criteria.
    **  The return vector uscratch has first the number of categories
    **   that were found (call it m), then m-1 goodnesses, then m labels
    **   in order, and the assurance that the best split is one of
    **   those that use categories in that order.
    */
    for(i = 0; i < nclass; i++)
      csplit[i] = 0;
    best = 0;
    m = uscratch[0];
    dscratch = uscratch + m;
    
    where = -1;
    left_n = 0;
    for(i = 1; i < m; i++) {
      k = dscratch[i - 1];   /* the next group of interest */
      for(j = 0; j < n; j++)
        if(x[j] == k)
          left_n++;
      right_n = n - left_n;
      if(right_n < edge)
        break;
      if((where < 0) || (uscratch[i] > best)) {
        best = uscratch[i];
        where = i;
      }
    }
    /*
    ** Now mark the groups as to left/right
    **   If there was no way to split it with at least 'edge' in each
    **   group, best will still =0.
    */
    if(best > 0) {
      for(i = 0; i < m ; i++) {
        k = dscratch[i];     /* the next group of interest */
        if(i < where)
          csplit[k - 1] = LEFT;
        else
          csplit[k - 1] = RIGHT;
      }
    }
  }
  *improve =  best;
  
  return;
}
/*
**  We don't do in-C cross validation for user splits, so there
**    is no prediction routine.
**  (Because of the structure of the calls, it's faster to make
**    use of xpred.rpart for user-written split routines).
*/
double usersplit_pred(
    double *y,
    double *yhat) {

  return(0.0);
}
