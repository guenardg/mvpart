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
 * SCCS @(#)rpartS.h    1.6 02/24/00
 * The S.h file defines a few things that I need, and hundreds that I don't.
 * In particular, on some architectures, it defines a variable "time"
 * which of course conflicts with lots of my C-code, 'time' being a natural
 * variable name for survival models and thus used in the poisson routines.
 * 
 * Thanks to Brian Ripley for suggesting a machine independent way of
 * fixing this.
 *
 * The S_alloc function changed it's argument list from version 4 to 5, and
 * the ALLOC macro allows me to have common C code for the two versions,
 * with only this file "survS.h" changed.
 *
 * *************************************************************************/

/*
**   
*/

#ifdef __GLIBC__
# define _POSIX_C_SOURCE 200809L
#endif

#define time timexxx
#include "R.h"
#undef time
#undef error

/*
** Memory defined with S_alloc is removed automatically
**  That with "CALLOC" I have to remove myself.  Use the
**  latter for objects that need to persist between the 
**  s_to_rp1 and s_to_rp2 calls
*/

#define ALLOC(a,b) S_alloc(a,b)
#define CALLOC(a,b) R_chk_calloc((size_t)(a), b)
