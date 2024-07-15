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
 * SCCS @(#)func_table.h    1.5 06/06/01
 * The table of implimented splitting functions
 * 
 *   init_split   - Will be called before a tree is started.  May do very
 *                   little, but empirical Bayes like methods will need to set
 *                   some global variables.
 *   choose_split - function to find the best split
 *   eval         - function to calculate the response estimate and risk
 *   error        - Function that returns the prediction error.
 *   num_y        - Number of columns needed to represent y (usually 1)
 *
 * *************************************************************************/

#ifdef __GLIBC__
# define _POSIX_C_SOURCE 200809L
#endif

/*
 * Prototype for the split initialization functions:
 * [name]init(int,double**,int,char**,double*,int*,int,double*)
 */

extern int anovainit(
    int n, double *y[], int maxcat, char **error, double *parm, int *size,
    int who, double *wt);

extern int mrtinit(
    int n, double *y[], int maxcat, char **error, double *parm, int *size,
    int who, double *wt);

extern int distinit(
    int n, double *y[], int maxcat, char **error, double *parm, int *size,
    int who, double *wt);

extern int poissoninit(
    int n, double *y[], int maxcat, char **error, double *parm, int *size,
    int who, double *wt);

extern int giniinit(
    int n, double *y[], int maxcat, char **error, double *parm, int *size,
    int who, double *wt);

extern int usersplit_init(
    int n, double *y[], int maxcat, char **error, double *parm, int *size,
    int who, double *wt);

/*
 * Prototype for the split choosing functions:
 * _name_(int,double**,double*,int,int,double*,double*,int*,double,double*)
 */

extern void anova(
    int n, double *y[], double *x, int nclass, int edge, double *improve,
    double *split, int *csplit, double myrisk, double *wt);

extern void mrt(
    int n, double *y[], double *x, int nclass, int edge, double *improve,
    double *split, int *csplit, double myrisk, double *wt);

extern void poisson(
    int n, double *y[], double *x, int nclass, int edge, double *improve,
    double *split, int *csplit, double myrisk, double *wt);

extern void gini(
    int n, double *y[], double *x, int nclass, int edge, double *improve,
    double *split, int *csplit, double myrisk, double *wt);

extern void dist(
    int n, double *y[], double *x, int nclass, int edge, double *improve,
    double *split, int *csplit, double myrisk, double *wt);

extern void usersplit(
    int n, double *y[], double *x, int nclass, int edge, double *improve,
    double *split, int *csplit, double myrisk, double *wt);

/*
 * Prototype for the evaluation functions:
 * [name][ss, dev, or eval](int,double**,double*,double*,double*)
 */

extern void anovass(
    int n, double *y[], double *value, double *risk, double *wt);

extern void mrtss(
    int n, double *y[], double *value, double *risk, double *wt);

extern void poissondev(
    int n, double *y[], double *value, double *risk, double *wt);

extern void ginidev(
    int n, double *y[], double *value, double *risk, double *wt);

extern void distss(
    int n, double *y[], double *value, double *risk, double *wt);

extern void usersplit_eval(
    int n, double *y[], double *value, double *risk, double *wt);

/*
 * Prototype for the error functions:
 * [name]pred(double*,double*)
 */

extern double anovapred(double *y, double *yhat);

extern double mrtpred(double *y, double *yhat);

extern double poissonpred(double *y, double *yhat);

extern double ginipred(double *y, double *yhat);

extern double distpred(double *y, double *yhat);

extern double usersplit_pred(double *y, double *yhat);

// The function table:

static struct {
  int (*init_split)(int,double**,int,char**,double*,int*,int,double*);
  void (*choose_split)(int,double**,double*,int,int,double*,double*,int*,double,
        double*);
  void (*eval)(int,double**,double*,double*,double*);
  double (*error)(double*,double*);
} func_table[] =
  {{anovainit,      anova,     anovass,        anovapred},
   {mrtinit,        mrt,       mrtss,          mrtpred},
   {poissoninit,    poisson,   poissondev,     poissonpred},
   {giniinit,       gini,      ginidev,        ginipred},
   {distinit,       dist,      distss,         distpred},
   {usersplit_init, usersplit, usersplit_eval, usersplit_pred}};

#define NUM_METHODS 6   /*size of the above structure */
