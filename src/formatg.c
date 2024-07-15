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
 * SCCS @(#)formatg.c  1.2 06/06/01
 * Write a bunch of numbers in a desired C format, (almost always %g)
 * GG: previous version used sprintf from stdio.h, which is regarded as unsafe
 * on MacOS platforms: issue solved by changing to snprintf from cstdio.
 *
 * *************************************************************************/

#include "rpartS.h"

#ifdef WIN32
#include <ctype.h>
#endif

void formatgC(
    int *n,
    double *x,
    char **format,
    char **out) {
  
  int i;
  int len;
  
#ifdef WIN32
  char *p;
#endif
  
  for(i = 0; i < *n; i++) {
    len = strlen(out[i]);
    snprintf(out[i], len, format[i], x[i]);
    
#ifdef WIN32
    /* change e+/-00n to e+/-0n etc */
    p = out[i];
    if(tolower(p[len - 5]) == 'e' && (p[len - 4] == '+' || p[len - 4] == '-') &&
       p[len - 3] == '0' && isdigit(p[len - 2]) && isdigit(p[len - 1])) {
       p[len - 3] = p[len - 2];
       p[len - 2] = p[len - 1];
       p[len - 1] = '\0';
    }
#endif
    
  }
}
