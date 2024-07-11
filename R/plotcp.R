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
##    **Function plotcp**
##
## **************************************************************************
##
#' 
#' Plot a Complexity Parameter Table for an Rpartmv Fit
#' 
#' Gives a visual representation of the cross-validation results in an
#' \code{\link{rpartmv-class}} object.
#' 
#' @name plotcp
#' 
#' @param x An object of class \code{rpartmv}.
#' @param xvse Multiplier for xvse * SE above the minimum of the curve.
#' @param minline Whether a horizontal line is drawn 1SE above the minimum of
#' the curve.
#' @param lty Type of lines used in plot, default = 3 (a dotted line).
#' @param col Color of lines, default = 1 (black).
#' @param upper What is plotted on the top axis: the size of the tree (the
#' number of leaves), the number of splits or nothing.
#' @param tab Used for multiple cross-validation.
#' @param resub.err Uses the re-substitution error for the calculation of SEs.
#' @param adj.df Adjusts the degrees of freedom of the re-substitution error
#' estimate in the calculation of SEs.
#' @param ... Additional graphical parameters
#' 
#' @returns \code{NULL} (invisibly).
#' 
#' @details The set of possible cost-complexity prunings of a tree from a nested
#' set. For the geometric means of the intervals of values of \code{cp} for
#' which a pruning is optimal, a cross-validation has (usually) been done in
#' the initial construction by \code{\link{rpartmv}}. Member \code{cptable} of
#' the \code{\link{rpartmv-class}} contains the mean and standard deviation of
#' the errors in the cross-validated prediction against each of the geometric
#' means, and these are plotted by this function. A good choice of \code{cp} for
#' pruning is often the leftmost value for which the mean lies below the
#' horizontal line.
#' 
#' @seealso \code{\link{rpartmv}}, \code{\link{printcp}},
#' \code{\link{rpartmv-class}}, and \code{\link{mvpart}} (for an example).
#' 
#' @keywords tree
#' 
#' @importFrom graphics abline axis box par points segments text
#' 
#' @export
plotcp <- function(x, xvse = 1, minline = TRUE , lty = 3, col = 1,
                   upper = c("size", "splits", "none"), tab, resub.err = TRUE ,
                   adj.df = FALSE , ...) {
  if(!inherits(x, "rpartmv"))
    stop("Not legitimate rpartmv object")
  upper <- match.arg(upper)
  p.rpartmv <- x$cptable
  if(xv <- (ncol(p.rpartmv) > 3L)) {
    xstd <- p.rpartmv[,5L]
    xerror <- p.rpartmv[,4L]
  }
  error <- p.rpartmv[,3L]
  nsplit <- p.rpartmv[,2L]
  ns <- seq(along = nsplit)
  cp0 <- p.rpartmv[,1L]
  cp <- sqrt(cp0 * c(Inf, cp0[-length(cp0)]))
  if (xv) {
    ylo <- min(c(xerror - xstd, error)) - 0.05
    yhi <- max(c(xerror + xstd, error)) + 0.05
  } else {
    ylo <- min(error) - 0.05
    yhi <- max(error) + 0.05
  }
  ylim <- c(ylo, yhi)
  plot(ns, error, axes=FALSE , xlab="cp", ylab="X-val Relative Error",
       ylim=ylim, type="n", ...)
  if(xv) {
    inpt <- (xerror == min(xerror))
    points(ns[inpt], xerror[inpt], col = "red", pch = 16L, cex = 2)
    inpt <- min(ns[xerror < min(xerror + xvse * xstd)])
    points(ns[inpt], xerror[inpt], col = "orange", pch = 16L, cex = 2)
    points(ns, xerror, type = "b", col = "blue", ...)
    segments(ns, xerror - xstd, ns, xerror + xstd, col = "blue", ...)
  }
  if(resub.err)
    points(ns, error, type="b", lty=1L, col = "darkgreen", ...)
  box()
  axis(2L, ...)
  axis(1L, at=ns, labels=as.character(signif(cp, 2)), ...)
  if (!missing(tab)) {
    xp <- as.numeric(names(tab))
    segments(ns[match(xp, nsplit + 1)], yhi, ns[match(xp,nsplit + 1)],
             yhi - 0.5*(tab/sum(tab))*(yhi - ylo), col=col + 1, lwd=2L, ...)
  }
  switch(
    upper,
    size = {
      axis(3L, at=ns, labels=as.character(nsplit + 1), ...)
      mtext("Size of tree", side=3L, line=3, cex=par()$cex, ...)
    },
    splits = {
      axis(3L, at=ns, labels=as.character(nsplit), ...)
      mtext("Number of splits", side=3L, line=3, ...)
    },
  )
  if(xv) {
    minpos <- min(seq(along = xerror)[xerror == min(xerror)])
    if(minline) {
      abline(h = (xerror + xvse * xstd)[minpos], lty=1L, col=col, xpd=FALSE)
      text(ns[2], (xerror + 0.5*xvse*xstd)[minpos],
           paste("Min +", xvse, "SE"), col=col, ...)
    }
  }
  invisible()
}
