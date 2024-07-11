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
 * SCCS @(#)node.h  1.3 12/13/99
 * definition of a node in the tree
 *
 * *************************************************************************/

#ifdef __GLIBC__
# define _POSIX_C_SOURCE 200809L
#endif

#ifdef MAINRP
#define EXTERN
#else
#define EXTERN extern
#endif

struct split {
  double improve;
  double adj;               /* for surrogates only, adjusted agreement */
  double spoint;            /* only used if it is continuous */
  struct split *nextsplit;
  int var_num;
  int count;
  int csplit[2];    /* the actual length will be longer for a categorical */
};                  /* predictor with >2 levels */

struct node {
  double  risk;             /* risk for the node */
  double complexity;        /* complexity at which it will collapse */
  double sum_wt;            /* sum of the weights for the node  */
  struct split *primary;
  struct split *surrogate;
  struct node *rightson;
  struct node *leftson;
  int num_obs;
  int lastsurrogate;
  double response_est[2];   /* actual length depends on splitting rule */
                            /*-- changed form 2 to 400 GD ?? --temp */
};

EXTERN int nodesize;

struct cptable {
  double cp;
  double risk;
  double xrisk;
  double xstd;
  int nsplit;
  struct cptable *forward;
  struct cptable *back;
};

/**************************************************************************
**
**  Split:
**      variable number of the split; 0 = no more surrogates (or primaries)
**
**      split point: the actual split point for a continuous
**
**      improve:  For primary splits, the iprovement index returned by the
**                 bsplit routine.  This is the measure that determines the
**                 winning split.
**                For surrogate splits, this holds the error rate, i.e., the
**                 % incorrect guesses of the primary by using this surrogate.
**
**      count: The number of observations split using this variable.  For the
**             first primary, this will = the number of non-missing values.
**             For surrogates, it will be the number missing in the primary
**             and all earlier surrogates but not missing on this one.  (For
**             all primaries but the first, the number is theoretical).
**
**  adj:  Let "maj" be the %agreement for going with the majority,
**                and "agree" the %agreement for this surrogate.  The
**                adjusted value is (agree - maj)/(1-maj); the amount of
**                the potential improvement actually realized.  The denominator
**                for both percents depends on the sur_agree option.
**
**      csplit[0]:   For a continuous variable, we also need to know the
**                    direction of the split.  We use this "extra" variable
**                    as 1: <x to the left, -1: <x to the right.
**
**      csplit[]:    For a categorical, the labels are LEFT, RIGHT, and
**                    0=missing.  (Even if a particular category is not empty,
**                    there may be no subjects from that category present
**                    at a particular split further down the tree).
**
**
**  Node:
**      num_obs: Number of observations in the node.
**
**      response_est: From the eval routine.  Estimate of the response, if
**                      this node were terminal.
**
**      risk: From the eval routine. Estimate of risk, if this node were
**                      terminal.
**
**      complexity: On the way down, it holds equation 5.18.  On the way up
**              it holds a provisional C.P. (The actual C.P. for each node
**              will be the minimum of this number and the provisional C.P.
**              of all nodes above it.  One more pass downward can establish
**              the proper C.P.).
**
**      lastsurrogate: Which direction to send obs for which the primary and
**              all the surrogates are missing.  (The child with the greatest
**      sum of weights).
*/
