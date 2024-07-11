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
##    **Function rpconvert**
##
## **************************************************************************
##
#'
#' Update an Rpartmv Object
#' 
#' Objects of class rpartmv changed (slightly) in their internal format in order
#' to accommodate the changes for user-written split functions.  This routine
#' updates an old object to the new format.
#' 
#' @name rpconvert
#' 
#' @param x An \code{\link{rpartmv-class}} object.
#' 
#' @returns An updated rpartmv object
#' 
#' @seealso \code{\link{rpartmv}}
#' 
#' @keywords tree
#' 
#' @export
rpconvert <- function(x) {
  if (!inherits(x, "rpartmv"))
    stop("x does not appear to be an rpartmv object")
  ff <- x$frame
  if(is.null(ff$splits)) {
    # this appears to be a new style one already
    warning("x not converted")
    return(x)
  }
  ff$splits <- NULL
  ff$wt <- ff$n
  xlev <- attr(x, "xlevels")
  if(length(xlev) > 0L) {
    zz <- as.numeric(names(xlev))
    names(xlev) <- attr(x$terms, "term.labels")[zz]
    attr(x, "xlevels") <- xlev
  }
  if(x$method == "class") {
    temp <- cbind(ff$yval, ff$yval2, ff$yprob)
    dimnames(temp) <- NULL
    ff$yval2 <- temp
    ff$yprob <- NULL
    x$frame <- ff
    temp <- rpartmv.class(c(1,1,2,2), NULL, wt=c(1,1,1,1))#dummy call
    x$functions <- list(summary=temp$summary, print=temp$print, text=temp$text)
  } else if (x$method=="anova") {
    x$frame <- ff
    temp <- rpartmv.anova(1:5, NULL, wt=rep(1,5))#dummy call
    x$functions <- list(summary=temp$summary, text = temp$text)
  } else {
    #either exp or poisson (they have the same summary/text pair)
    ff$yval2 <- cbind(ff$yval, ff$yval2)
    x$frame <- ff
    temp <- rpartmv.poisson(1:5, NULL, wt=rep(1,5))#dummy call
    x$functions <- list(summary=temp$summary, text = temp$text)
  }
  class(x) <- "rpartmv"
  x
}
