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
 * SCCS @(#)rpartexp2.c 1.1 07/05/01
 * Cuts down a list of death times, to avoid ones that differ by
 * very, very little.
 *  n   number of death times
 *  y   list of death times, sorted
 *  eps     machine precision
 * output
 *   keep    1=keep this one, 0=don't
 * *************************************************************************/

#include "rpart.h"
#include "rpartproto.h"

void rpartexp2(
    int *n2,
    double *y,
    double *eps,
    int *keep) {
  
  int n;
  double delta;
  int i, j;
  double lasty;
  n = *n2;
  
  /* let delta = eps * interquartile range */
  i = n/4;
  j = (3*n)/4;
  delta = (*eps)*(y[j] - y[i]);
  /*
  ** make sure that each y is at least "delta" greater than
  ** the last y that we have decided to keep
  */
  lasty = y[0];
  keep[0] = 1;
  for(i = 1; i < n; i++) {
    if((y[i] - lasty) <= delta)
      keep[i] =0;
    else {
      keep[i] =1;
      lasty = y[i];
    }
  }
  
  return;
}
