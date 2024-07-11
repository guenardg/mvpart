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
 * SCCS @(#)pred_rpart.c   1.6  06/06/01
 * Do rpart predictions given the matrix form of the tree.
 * Input
 *   dimx        : # of rows and columns in the new data
 *   nnode       : # of nodes in the tree
 *   nsplit      : # of split structures
 *   dimc        : dimension of the categorical splits matrix
 *   nnum        : node number for each row of 'nodes'
 *   nodes       : matrix of node info
 *                 row 0= count, 1=index of primary, 2=#competitors,
 *                     3= number of surrogates
 *   vnum        : variable number of each split
 *   split       : matrix of split info
 *                 row 0=useage count, 1= #categories if >1, otherwise
 *                     the split parity, 2= utility, 3= index to csplit
 *                     or numeric split point
 *   csplit      : matrix of categorical split info
 *   usesur      : at what level to use surrogates
 *   xdata       : the new data
 *   xmiss       : shows missings in the new data
 *
 * Output
 *   where       : the "final" row in nodes for each observation
 *
 * *************************************************************************/

#include "rpartS.h"
#include "rpart.h"
#include "rpartproto.h"

void pred_rpart(
    int *dimx,
    int *nnode,
    int *nsplit,
    int *dimc,
    int *nnum,
    int *nodes2,
    int *vnum,
    double *split2,
    int *csplit2,
    int *usesur,
    double *xdata2,
    int *xmiss2,
    int *where) {
  
  int i,j;
  int n;
  int ncat;
  int node, nspl, var, dir;
  int lcount, rcount;
  int npos;
  double temp;
  int *nodes[4];
  double *split[4];
  int **csplit = NULL, **xmiss;
  double **xdata;
  n = dimx[0];
  for(i = 0; i < 4; i++) {
    nodes[i] = &(nodes2[(*nnode)*i]);
    split[i] = &(split2[(*nsplit)*i]);
  }
  if(dimc[1] > 0) {
    csplit = (int**)ALLOC((int)dimc[1], sizeof(int*));
    for(i = 0; i < dimc[1]; i++)
      csplit[i] = &(csplit2[i * dimc[0]]);
  }
  xmiss = (int**)ALLOC((int)dimx[1], sizeof(int*));
  xdata = (double**)ALLOC((int)dimx[1], sizeof(double*));
  for(i = 0; i < dimx[1]; i++) {
    xmiss[i] = &(xmiss2[i*dimx[0]]);
    xdata[i] = &(xdata2[i*dimx[0]]);
  }
  for(i = 0; i < n; i++) {
    node = 1;   /*current node of the tree */
    
next:
    
    for(npos = 0; nnum[npos] != node; npos++);  /*position of the node */
    /* walk down the tree */
    nspl = nodes[3][npos] - 1;    /*index of primary split */
    if(nspl >= 0) {               /* not a leaf node */
      var = vnum[nspl] - 1;
      if(xmiss[var][i] == 0) {    /* primary var not missing */
        ncat = split[1][nspl];
        temp = split[3][nspl];
        if(ncat >= 2)
          dir = csplit[(int)xdata[var][i] - 1][(int)temp - 1];
        else if(xdata[var][i] < temp)
          dir = ncat;
        else
          dir = -ncat;
        if(dir != 0) {
          if(dir == -1)
            node = 2*node;
          else
            node = 2*node + 1;
          goto next;
        }
      }
      if(*usesur > 0) {
        for (j = 0; j < nodes[2][npos]; j++) {
          nspl = nodes[1][npos] + nodes[3][npos] + j;
          var = vnum[nspl] - 1;
          if(xmiss[var][i] == 0) {     /* surrogate not missing */
            ncat = split[1][nspl];
            temp = split[3][nspl];
            if(ncat >= 2)
              dir = csplit[(int)xdata[var][i] - 1][(int)temp - 1];
            else if(xdata[var][i] < temp)
              dir = ncat;
            else
              dir = -ncat;
            if(dir != 0) {
              if(dir== -1)
                node = 2*node;
              else
                node = 2*node +1;
              goto next;
            }
          }
        }
      }
      
      if(*usesur > 1) { /* go with the majority */
        for(j = 0; nnum[j] != (2*node); j++);
        lcount = nodes[0][j];
        for(j = 0; nnum[j] != (1 + 2*node); j++);
        rcount = nodes[0][j];
        if(lcount != rcount) {
          if(lcount > rcount)
            node = 2*node;
          else
            node = 2*node + 1;
          goto next;
        }
      }
    }
    where[i] = npos + 1;
  }
}
