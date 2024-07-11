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
 * xdiss.c
 *
 * *************************************************************************/

void xdists(
    double *dist,
    int *nn,
    double *dcrit,
    int *usemin,
    double *eps,
    double *bigd) {
  
  double maxd, mind, sum, big = *bigd, dki, dkj, adj, mink;
  int i, j, k, ki, kj, pos, n = *nn, ndist, nbig, nbigold, nloop, count;
  
  ndist = n*(n - 1)/2;
  nbig = 0;
  maxd = 0;
  adj = 0;
  for(i = 0; i < ndist; i++)
  if(dist[i] > *dcrit - *eps) {
    dist[i] =  -1;
    nbig += 1;
  }
  nloop = 0;

  while(nbig > 0) {
    nbigold = nbig;
    nloop += 1;
    if(nloop == 1)
      adj = *dcrit;
    else if(nloop>1)
      adj = maxd;
    maxd = 0.0;
    mind = 1e9;
    
    for(j = 1; j < n; j++)
      for(i = 0; i < j; i++) {
        pos = n*i - i*(i + 1)/2 + j - i - 1;
        if(dist[pos] < 0) {
          mink = 1e9;
          sum = 0.0;
          count = 0;
          for(k = 0; k < n; k++) {
            if((k != i) & (k != j)) {
              if(k > i)
                ki = n*i - i*(i+1)/2 + k - i - 1;
              else
                ki = n*k - k*(k+1)/2 + i - k - 1;
              if(k > j)
                kj = n*j-j*(j+1)/2+k-j-1;
              else
                kj = n*k-k*(k+1)/2+j-k-1;
              dki = dist[ki];
              dkj = dist[kj];
              if((dki >= 0) & (dkj >= 0) & (dki < big) & (dkj < big)) {
                if((dki + dkj) < mink)
                  mink = dki + dkj;
                sum += dki + dkj;
                count += 1;
              }
            }
          }
          if(count > 0) {
            if(*usemin > 0)
              dist[pos] = mink + big;
            else dist[pos] = sum/count + big;
            if(dist[pos] > maxd)
              maxd = dist[pos];
            if(dist[pos] < mind)
              mind = dist[pos];
          }
        }
      }
    
    maxd = maxd-big;
    mind = mind-big;
    nbig = 0;
    for(i = 0; i < ndist; i++)
      if(dist[i] < 0)
        nbig += 1;
    if(nbigold == nbig)
      nbig=0;
    
    for(i = 0; i < ndist; i++) {
      if(dist[i] > big) {
        if(*usemin > 0)
          dist[i] = dist[i] - big;
        else
          dist[i] = dist[i] - big - mind + adj + *eps;
      }
    }
    maxd = maxd - mind + adj;
  }
}
