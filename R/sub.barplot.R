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
##    **Internal function sub.barplot**
##
## **************************************************************************
##

sub.barplot <- function (x, y, z, keep = rep(TRUE , length(x)),
                         row.scale = FALSE , xadj = 1, yadj = 1, bord = TRUE,
                         line = TRUE , col = col) {
  par(xpd = TRUE)
  drawbar <- function(x, y, z, line.adj = 0, xwid, ywid, border = TRUE,
                      line = TRUE, colbar = 10L:12L) {
    xx <- c(x - xwid/2,x + xwid/2,x + xwid/2,x - xwid/2)
    yy <- c(y - ywid,y - ywid, y,y)
    nbar <- length(z)
    xbwid <- xwid/nbar
    for(i in 1:nbar) {
      xb <- x - xwid/2 + xbwid * c(i - 1,i,i,i - 1)
      yb <- y - ywid + z[i]*ywid*c(0,0,1,1) - line.adj*ywid
      polygon(xb, yb, col = colbar[i])
    }
    if(border)
      polygon(xx, yy, col = 1, density = 0)
    if(line)
      lines(xx[1L:2L], yy[1L:2L] - line.adj * ywid)
  }
  xrnge <- range(x)
  yrnge <- range(y)
  n <- length(x)
  xdiv <- max(sum(keep),6)
  xwid <- (xadj * diff(xrnge))/xdiv
  ywid <- (yadj * diff(yrnge))/12
  x <- x[keep]
  y <- y[keep]
  z <- z[keep,]
  nkeep <- sum(keep)
  nz <- ncol(z)
  if(length(col) < nz)
    col <- rep(col, nz, length = nz)
  if(any(z < 0)) {
    z <- z/diff(range(z))
    ladj <- min(z)
  } else {
    if(row.scale) {
      z <- z/apply(z, 1, sum)
    } else z <- z/max(z)
    ladj <- 0
  }
  for(i in 1:nkeep)
    drawbar(x[i], y[i], z[i,], line.adj=ladj, xwid, ywid, border=bord,
            line=line, colbar=col)
}
