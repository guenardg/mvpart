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
 * GG 2024-05-01
 * C routines registration
 *
 * *************************************************************************/

#include <R.h>
#include <Rinternals.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* .C calls */
extern void formatgC(int*, double*, char**, char**);
extern void gdistance(double*, int*, int*, double*, int*, int*);
extern void pred_rpart(int*, int*, int*, int*, int*, int*, int*, double*, int*,
                       int*, double*, int*, int*);
extern void rpartexp2(int*, double*, double*, int*);
extern void s_to_rp(int*, int*, int*, int*, double*, double*, int*, int*,
                    double*, double*, int*, int*, char**, double*, int*,
                    double*);
extern void s_to_rp2(int*, int*, int*, int*, int*, int*, int*, int*, double*,
                     double*, int*, int*, double*, int*);
extern void xdists(double*, int*, double*, int*, double*, double*);
extern void s_xpred(int*, int*, int*, int*, double*, double*, int*, int*,
                    double*, double*, int*, double*, int*, double*, char**,
                    double*, int*, double*);

/* .Call calls */
extern SEXP init_rpcallback(SEXP, SEXP, SEXP, SEXP, SEXP);

static const R_CMethodDef CEntries[] = {
  {"formatgC",        (DL_FUNC) &formatgC,         4},
  {"gdistance",       (DL_FUNC) &gdistance,        6},
  {"pred_rpart",      (DL_FUNC) &pred_rpart,      13},
  {"rpartexp2",       (DL_FUNC) &rpartexp2,        4},
  {"s_to_rp",         (DL_FUNC) &s_to_rp,         16},
  {"s_to_rp2",        (DL_FUNC) &s_to_rp2,        14},
  {"xdists",          (DL_FUNC) &xdists,           6},
  {"s_xpred",         (DL_FUNC) &s_xpred,         18},
  {NULL,              NULL,                        0}
};

static const R_CallMethodDef CallEntries[] = {
  {"init_rpcallback", (DL_FUNC) &init_rpcallback,  5},
  {NULL,              NULL,                        0}
};

void R_init_mvpart(DllInfo *dll)
{
  R_registerRoutines(dll, CEntries, CallEntries, NULL, NULL);
  R_useDynamicSymbols(dll, FALSE);
}
