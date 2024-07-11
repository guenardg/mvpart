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
 * SCCS @(#)make_cp_table.c 1.2 02/08/98
 * 
 * Given a cptable list already initialized with the unique cp's in it,
 * fill in the columns for risk and number of splits.
 *
 * The basic logic is: for each terminal node on the tree, start myself
 * down at the bottom of the list of complexity parameters. For each
 * unique C.P. until my parent collapses, the node I'm in adds into that
 * line of the CP table. So walk up the CP list, adding in, until my
 * parent would collapse; then report my position in the cp list to the
 * parent and quit.
 *
 * parent: complexity of the parent node
 *
 * *************************************************************************/

/*

*/
#include "rpart.h"
#include "node.h"
#include "rpartproto.h"

struct cptable* make_cp_table(
    struct node *me,
    double parent,
    int nsplit) {
  
  struct cptable *cplist;
  
  if (me->leftson) {  /* if there are splits below */
    /*
    ** The 2 lines below are perhaps devious
    **  1) Since the return value depends on ones parent, both calls will
    **       return the same thing.
    **  2) I send 0 to the left to keep the current split (me) from
    **       being counted twice, once by each child.
    */
    make_cp_table(me->leftson, me->complexity, 0);
    cplist = make_cp_table(me->rightson, me->complexity, nsplit+1);
  }
  else cplist = cptable_tail;
  
  while (cplist->cp < parent) {
    cplist->risk += me->risk;
    cplist->nsplit += nsplit;
    cplist = cplist->back;
  }
  
  return(cplist);
}
