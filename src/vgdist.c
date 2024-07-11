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
 * Author: Jari Oksanen : modified Glenn De'ath
 * Distance measures for vegetation sciences.
 * The measures here were recommended by Peter Minchin, since
 * they have a good rank-order relation with gradient distance.
 * The standard distances are found in standard R library 
 * mva in function dist (distance.c).
 *
 * *************************************************************************/

#include <R.h>
#include <math.h>

/* Indices */

#define MANHATTAN 1
#define EUCLIDEAN 2
#define CANBERRA 3
#define BRAY 4
#define KULCZYNSKI 5
#define GOWER 6
#define MAXIMUM 7
#define BINARY 8

/* Distance functions */

double g_manhattan(
    double *x,
    int nr,
    int nc,
    int i1,
    int i2) {
  
  double dist;
  int count, j;
  
  dist = 0.0;
  count = 0;
  for(j = 0; j < nc; j++) {
    if(R_FINITE(x[i1]) && R_FINITE(x[i2])) {
      dist += fabs(x[i1] - x[i2]);
      count++;
    }
    i1 += nr;
    i2 += nr;
  }
  if(count == 0)
    dist = NA_REAL;
  return dist;
}

/* Gower is like Manhattan, but data were standardized to 
   range 0..1 for rows before call and dist is divided
   by the number of non-zero pairs 
*/

double g_gower(
    double *x,
    int nr,
    int nc,
    int i1,
    int i2) {
  
  double dist;
  int count, j;
  
  dist = 0.0;
  count = 0;
  for(j = 0; j < nc; j++) {
    if(R_FINITE(x[i1]) && R_FINITE(x[i2])) {
      dist += fabs(x[i1] - x[i2]);
      count++;
    }
    i1 += nr;
    i2 += nr;
  }
  dist /= (double)count;
  return dist;
}

double g_euclidean(
    double *x,
    int nr,
    int nc,
    int i1,
    int i2) {
  
  double dist, dev;
  int count, j;

  count = 0;
  dist = 0.0;
  for(j = 0; j < nc; j++) {
    if(R_FINITE(x[i1]) && R_FINITE(x[i2])) {
      dev = x[i1] - x[i2];
      dist += dev*dev;
      count++;
    }
    i1 += nr;
    i2 += nr;
  }
  if(count == 0)
    dist = NA_REAL;
  return sqrt(dist);
}

double g_canberra(
    double *x,
    int nr,
    int nc,
    int i1,
    int i2) {
  
  double numer, denom, dist;
  int count, j;

  count = 0;
  dist = 0.0;
  for(j = 0; j < nc; j++) {
    if(R_FINITE(x[i1]) && R_FINITE(x[i2])) {
      if((x[i1] != 0) || (x[i2] != 0)) {
        count++;
        denom = x[i1] + x[i2];
        if(denom > 0.0) {
          numer = fabs(x[i1] - x[i2]);
          dist += numer/denom;
        }
        else
          dist += R_PosInf;
      }
    }
    i1 += nr;
    i2 += nr;
  }
  if(count == 0)
    dist = NA_REAL;
  dist /= (double)count;
  return dist;
}

double g_bray(
    double *x,
    int nr,
    int nc,
    int i1,
    int i2) {
  
  double dist, total;
  int count, j;
  
  total = 0.0;
  count = 0;
  dist = 0;
  for(j = 0; j < nc; j++) {
    if(R_FINITE(x[i1]) && R_FINITE(x[i2])) {
      dist += fabs(x[i1] - x[i2]);
      total += x[i1] + x[i2];
      count++;
    }
    i1 += nr;
    i2 += nr;
  }
  if(count == 0)
    dist = NA_REAL;
  dist /= total;
  return dist;
}

double g_kulczynski(
    double *x,
    int nr,
    int nc,
    int i1,
    int i2) {
  
  double sim, dist, t1, t2;
  int count, j;
  
  t1 = 0.0;
  t2 = 0.0;
  count = 0;
  sim = 0.0;
  for(j = 0; j < nc; j++) {
    if(R_FINITE(x[i1]) && R_FINITE(x[i2])) {
      sim += (x[i1] < x[i2]) ? x[i1] : x[i2];
      t1 += x[i1];
      t2 += x[i2];
      count++;
    }
    i1 += nr;
    i2 += nr;
  }
  if(count == 0)
    dist = NA_REAL;
  dist = 1 - sim/t1/2 - sim/t2/2;
  return dist;
}

double g_maximum(
    double *x,
    int nr,
    int nc,
    int i1,
    int i2) {
  
  double dev, dist;
  int count, j;
  
  count = 0;
  dist = -10^10;
  for(j = 0; j < nc; j++) {
    if(R_FINITE(x[i1]) && R_FINITE(x[i2])) {
      dev = fabs(x[i1] - x[i2]);
      if(dev > dist)
        dist = dev;
      count++;
    }
    i1 += nr;
    i2 += nr;
  }
  if(count == 0)
    dist = NA_REAL;
  return dist;
}

double g_binary(
    double *x,
    int nr,
    int nc,
    int i1,
    int i2) {
  
  int count, dist;
  int j;
  
  count = 0;
  dist = 0;
  for(j = 0; j < nc; j++) {
    if(R_FINITE(x[i1]) && R_FINITE(x[i2])) {
      if(x[i1] || x[i2]) {
        count++;
        if(!(x[i1] && x[i2]))
          dist++;
      }
    }
    i1 += nr;
    i2 += nr;
  }
  if(count == 0)
    dist = NA_REAL;
  return (double)dist/count;
}

/* Driver */
static double (*distfun)(double*, int, int, int, int);

void gdistance(
    double *x,
    int *nr,
    int *nc,
    double *d,
    int *diag,
    int *method) {
  
  int dc, i, j, ij;
  switch(*method) {
  case MANHATTAN :
    distfun = g_manhattan;
    break;
  case EUCLIDEAN :
    distfun = g_euclidean;
    break;
  case CANBERRA :
    distfun = g_canberra;
    break;
  case BRAY :
    distfun = g_bray;
    break;
  case KULCZYNSKI :
    distfun = g_kulczynski;
    break;
  case GOWER :
    distfun = g_gower;
    break;
  case MAXIMUM :
    distfun = g_maximum;
    break;
  case BINARY :
    distfun = g_binary;
    break;
  }
  
  dc = (*diag) ? 0 : 1; 
  ij = 0;
  for (j = 0; j <= *nr; j++)
    for (i = j + dc; i < *nr; i++)
      d[ij++] = distfun(x, *nr, *nc, i, j);
}
