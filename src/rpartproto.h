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
 * SCCS  @(#)rpartproto.h  1.5 06/06/01
 * prototypes for all of the rpart functions
 * This helps the ansi compiler do tight checking.
 *
 * *************************************************************************/

struct node* branch(
    struct node *tree, int obs);

void bsplit(
    struct node *me, int nodenum);

void choose_surg(
    int nodenum, int *y, double *x, int *order, int ncat, double *agreement,
    double *split, int *csplit, double ltot, double rtot, double *adj);

void fix_cp(
    struct node *me, double parent_cp);

void free_tree(
    struct node *node, int freenode);

void graycode_init0(
    int maxcat);

void graycode_init1(
    int numcat, int *count);

void graycode_init2(
    int numcat, int *count, double *val);

int graycode();

struct split* insert_split(
    struct split **listhead, int ncat, double improve, int max);

void make_cp_list(
    struct node *me, double parent, struct cptable *cptable_head);

struct cptable* make_cp_table(
    struct node *me, double parent, int nsplit);

void mysort(
    int start, int stop, double *x, int *cvec);

void nodesplit(
    struct node *me, int nodenum);

int partition(
    int nodenum, struct node *splitnode, double *sumrisk);

void pred_rpart(
    int *dimx, int *nnode, int *nsplit, int *dimc, int *nnum, int *nodes2,
    int *vnum, double *split2, int *csplit2, int *usesur, double *xdata2,
    int *xmiss2, int *where);

int rpart(
    int n, int nvarx, int *ncat, int method, int  maxpri, double *parms,
    double *ymat, double *xmat, int dissim, int *missmat,
    struct cptable *cptable, struct node **tree, char **error, int *which,
    int xvals, int *x_grp, double *wt, double *opt, int ny, double *cost);

void rpart_callback0(
    int *nr);

void rpart_callback1(
    int n, double *y[], double *wt, double *z);

void rpart_callback2(
    int n, int ncat, double *y[], double *wt, double *x, double *good);

void rpcountup(
    struct node *me, int *nnode, int *nsplit, int *ncat);

void rplabel(
    int *nsplit, int *index, double *splits, int *ncat, int *csplit,
    char **cutleft, char **cutright);

void rpmatrix(
    struct node *me, int *nodecount, int *splitcount, int *catcount,
    int *numcat, double **dsplit, int **isplit, int **csplit, double **dnode,
    int **inode, int id);

void rundown(
    struct node *tree, int obs, double *cp, double *xpred, double *xtemp);

void rundown2(
    struct node *tree, int obs, double *cp, double *xpred);

void s_to_rp(
    int *n, int *nvarx, int *ncat, int *method, double *opt, double *parms,
    int *xvals, int *x_grp, double *y, double *xmat, int *dissim, int *missmat,
    char **error, double *wt, int *ny, double *cost);

void s_to_rp2(
    int *n, int *nsplit, int *nnode, int *ncat, int *numcat, int *maxcat,
    int *xvals, int *which, double *cptable, double *dsplit, int *isplit,
    int *csplit, double *dnode, int *inode);

void s_xpred(
    int *sn, int *nvarx, int *ncat, int *method, double *opt, double *parms,
    int *xvals, int *x_grp, double *ymat, double *xmat, int *missmat,
    double *predict, int *ncp, double *cp, char **error, double *wt, int *ny,
    double *cost);

void surrogate(
    struct node *me, int nodenum);

void xval(
    int n_xval, struct cptable *cptable_head, int *x_grp, int maxcat,
    char **error, double * parms);
