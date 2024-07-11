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
##    **Internal function formatg**
##
## **************************************************************************
##

formatg <- function(x, digits= unlist(options('digits')),
                    format= paste("%.", digits, "g", sep='')) {
  if (!is.numeric(x)) stop("x must be a numeric vector")
  
  n <- length(x)
  #
  # the resultant strings could be up to 8 characters longer,
  #   assume that digits = 4,  -0.dddde+104 is a worst case, where
  #   dddd are the 4 significant digits.
  #
  dummy  <- paste(rep(" ", digits + 8L), collapse='')
  .C("formatgC",
     as.integer(n),
     as.double(x),
     rep(format,n),
     out= rep(dummy, n),
     NAOK=TRUE,
     PACKAGE="mvpart")$out -> temp
  if (is.matrix(x)) matrix(temp, nrow=nrow(x)) else matrix(temp,nrow=1)
}
