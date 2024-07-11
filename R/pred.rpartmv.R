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
##    **Internal function pred.rpartmv**
##
## **************************************************************************
##

pred.rpartmv <- function(fit, x) {
  frame <- fit$frame
  if(nrow(frame) == 1) { # handle root-only tree separately
      temp <- rep(1, nrow(x))
  } else {
    nc <- frame[, c("ncompete","nsurrogate")]
    frame$index <- 1 +
      c(0, cumsum((frame$var != "<leaf>") + nc[[1]] + nc[[2]]))[-(nrow(frame) + 1)]
    frame$index[frame$var == "<leaf>"] <- 0
    vnum <- match(dimnames(fit$split)[[1]], dimnames(x)[[2]])
    if(any(is.na(vnum)))
      stop("Tree has variables not found in new data")
    .C("pred_rpart",
      as.integer(dim(x)),
      as.integer(dim(frame)[1]),
      as.integer(dim(fit$splits)),
      as.integer(if(is.null(fit$csplit)) rep(0,2) else dim(fit$csplit)),
      as.integer(row.names(frame)),
      as.integer(unlist(frame[,c("n","ncompete","nsurrogate","index")])),
      as.integer(vnum),
      as.double(fit$splits),
      as.integer(fit$csplit -2),
      as.integer((fit$control)$usesurrogate),
      as.double(x),
      as.integer(is.na(x)),
      where = integer(dim(x)[1]),
      NAOK = TRUE,
      PACKAGE = "mvpart")$where -> temp
  }
  names(temp) <- rownames(x)
  temp
}
