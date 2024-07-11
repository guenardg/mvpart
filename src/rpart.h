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
 * SCCS @(#)rpart.h 1.10 06/06/01
 * Commom variables for the rpart routine.
 *
 * *************************************************************************/

#ifdef __GLIBC__
# define _POSIX_C_SOURCE 200809L
#endif

#include <R.h>
#undef error
#define LEFT  (-1)     /*used for the variable "extra" in nodes */
#define RIGHT  1
#define MISSING 0

#ifdef MAINRP
#define EXTERN
#else
#define EXTERN extern
#endif
/* As a sop to S, I need to keep the total number of external symbols
**  somewhat smaller.  So, pack most of them all into a structure.
*/
EXTERN struct {
  double complexity;
  double alpha;
  double iscale;      /* used to check improvement==0, with error */
  double **ydata;
  double  **xdata;
  double  *xtemp;
  double *wt;
  double **ytemp;
  double *wtemp;      /* temp vector of weights */
  double *lwt;
  double *rwt;        /* scratch double vectors, of length ncat */
  double *vcost;      /* variable costs */
  int *numcat;        /* variable type: 0=cont, 1+  =#categories */
  int **sorts;        /* allocated on the fly */
  int n;              /* total number of subjects */
  int num_y;          /* number of y variables */
  int nvar;           /* number of predictors */
  int method;         /* method -- needed for dists -- GD 11/03 */
  int maxpri;
  int maxsur;         /* max # of primary or surrogate splits to use */
  int usesurrogate;
  int num_unique_cp;
  int min_node;       /* minimum size for any terminal node */
  int min_split;      /* minimum size before we attempt a split */
  int num_resp;       /* length of the response vector */
  int sur_agree;      /* 0 =  my style, 1=CART style */
  int maxnode;        /* controls the maximum depth of the tree */
  int *tempvec;       /* to be allocated by the mainline, of length n */
  int *which;
  int *csplit;
  int *left;
  int *right;
  int dissim;
}  rp;

EXTERN struct cptable *cptable_tail;

// Called to initialize a splitting function:
EXTERN int (*rp_init)(int,double**,int,char**,double*,int*,int,double*);
/* GG: Had to fix "s_xpred.c" and "xval.c" because of an incompatible pointer
 * issue (double* instead of int*) while calling the init function. */

// Set to the splitting function:
EXTERN void (*rp_choose)(int,double**,double*,int,int,double*,double*,int*,
             double,double*);

// Set to the evaluation routine:
EXTERN void (*rp_eval)(int,double**,double*,double*,double*);

// Set to the prediction error routine:
EXTERN double (*rp_error)(double*,double*);

/*
** The user inputs his complexity parameter as a percentage. and the
**   printout is also scaled in this way.  The book and the computations all
**   have an easier time with absolute cp.  So complex = what the user
**   typed and alpha = complex * (risk of top node) = what is used
**   internally.
** The variable 'complex' in node.h is also on the internal scale.
**
** Categorical variables must be coded as 1,2,3, ..., and there may be
**  missing categories.  The upper limit is determined on the fly.
*/
