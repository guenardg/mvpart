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
 * SCCS @(#)free_tree.c 1.2 02/08/98
 * free up all of the memory associated with a tree
 *
 * *************************************************************************/

#include "rpart.h"
#include "node.h"
#include "rpartS.h"
#include "rpartproto.h"

void free_tree(
    struct node *node,
    int freenode) {
  
  struct split *s1, *s2;
  
  if(node->rightson != 0)
    free_tree(node->rightson, 1);
  if(node->leftson  != 0)
    free_tree(node->leftson,  1);
  
  for(s1 = node->surrogate; s1 != 0;) {
    s2 = s1;
    s1 = s1->nextsplit;
    Free(s2);
  }
  for(s1 = node->primary; s1!=0;) {
    s2 = s1;
    s1 = s1->nextsplit;
    Free(s2);
  }
  if (freenode == 1) Free(node);
}
