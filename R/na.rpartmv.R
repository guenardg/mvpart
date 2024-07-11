## **************************************************************************
##
##    Actual version maintained by:
##    Guillaume Guénard <guillaume.guenard@gmail.com>
##    Department de sciences biologiques,
##    Université de Montréal
##    Montreal, QC, Canada
##    from 2024-04-13
##
##    Previous versions:
##    rpart by Terry M Therneau and Beth Atkinson <atkinson@mayo.edu>,
##    R port of rpart by Brian Ripley <ripley@stats.ox.ac.uk>,
##    Some routines added from vegan by Jari Oksanen <jari.oksanen@oulu.fi>,
##    Extensions and adaptations of rpart to mvpart by Glenn De'ath
##    <g.death@aims.gov.au>
##    Because the current version of rpart by cannot fit a multivariate
##    regression tree, this package has a function called rpartmv that is a
##    clone of a previous version of rpart which has been adapted by G. De'ath
##    to feature the multivariate regression tree as a possible model.
##
##    This file is part of mvpart
##
##    mvpart is free software: you can redistribute it and/or modify
##    it under the terms of the GNU General Public License as published by
##    the Free Software Foundation, either version 3 of the License, or
##    (at your option) any later version.
##
##    mvpart is distributed in the hope that it will be useful,
##    but WITHOUT ANY WARRANTY; without even the implied warranty of
##    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
##    GNU General Public License for more details.
##
##    You should have received a copy of the GNU General Public License
##    along with mvpart. If not, see <https://www.gnu.org/licenses/>.
##
##    **Internal function na.rpartmv**
##
## **************************************************************************
##

na.rpartmv <- function(x) {
  Terms <- attr(x, 'terms')
  if(!is.null(Terms)) yvar <- attr(Terms, "response") else yvar <- 0
  if(yvar == 0) {
    xmiss <- is.na(x)
    keep <-  (xmiss %*% rep(1,ncol(xmiss))) < ncol(xmiss)
  } else {
    xmiss <- is.na(x[-yvar])
    ymiss <- is.na(x[[yvar]])
    if(is.matrix(ymiss))
      keep <- ((xmiss %*% rep(1,ncol(xmiss))) < ncol(xmiss)) &
      ((ymiss %*% rep(1,ncol(ymiss))) == 0)
    else
      keep <- ((xmiss %*% rep(1,ncol(xmiss))) < ncol(xmiss)) & !ymiss
  }
  if(all(keep)) x else {
    temp <- seq(keep)[!keep]
    names(temp) <- row.names(x)[!keep]
    #the methods for this group are all the same as for na.omit
    class(temp) <- c("na.rpartmv", "omit")
    structure(x[keep,], na.action=temp)
  }
}
