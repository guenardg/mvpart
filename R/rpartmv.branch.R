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
##    **Internal function rpartmv.branch**
##
## **************************************************************************
##

rpartmv.branch <- function(x, y, node, branch) {
  # Draw a series of horseshoes, left son, up, over, down to right son
  #   NA's in the vector cause lines() to "lift the pen"
  is.left <- (node%%2 == 0)        #left hand sons
  node.left <- node[is.left]
  parent <- match(node.left/2, node)
  sibling <- match(node.left + 1, node)
  temp <- (x[sibling] - x[is.left])*(1 - branch)/2
  xx <- rbind(x[is.left], x[is.left] + temp, x[sibling] - temp, x[sibling], NA)
  yy <- rbind(y[is.left], y[parent], y[parent], y[sibling], NA)
  list(x=xx, y=yy)
}
