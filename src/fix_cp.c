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
 * SCCS @(#)fix_cp.c    1.2 02/08/98
 * When partition is done, each node is labeled with the complexity
 * appropriate if it were the top of the tree. Actually, the complexity
 * should be min(me, any-node-above-me). This routine fixes that.
 *
 * *************************************************************************/

#include "rpart.h"
#include "node.h"
#include "rpartproto.h"

void fix_cp(
    struct node *me,
    double parent_cp) {
  
  if(me->complexity > parent_cp)
    me->complexity = parent_cp;
  
  if(me->leftson != 0) {
    fix_cp(me->leftson, me->complexity);
    fix_cp(me->rightson, me->complexity);
  }
}
