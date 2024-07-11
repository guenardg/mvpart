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
##    **Internal function rpartmv.matrix**
##
## **************************************************************************
##

rpartmv.matrix <- function(frame) {
  if(!inherits(frame, "data.frame"))
    return(as.matrix(frame))
  frame$"(weights)" <- NULL
  terms <- attr(frame, "terms")
  if(is.null(terms)) {
    predictors <- names(frame)
  } else {
    a <- attributes(terms)
    predictors <- as.character(a$variables)[-1L] # R change
    removals <- NULL
    if((TT <- a$response) > 0L) {
      removals <- TT
      frame[[predictors[TT]]] <- NULL
    }
    if(!is.null(TT <- a$offset)) {
      removals <- c(removals, TT)
      frame[[predictors[TT]]] <- NULL
    }
    if(!is.null(removals))
      predictors <- predictors[ - removals]
    labels <- a$term.labels
    if(abs(length(labels)-length(predictors))>0)
      predictors <- predictors[match(labels,predictors)]
  }
  factors <- sapply(frame, function(x) !is.null(levels(x)))
  characters <- sapply(frame, is.character)
  if(any(factors | characters)) {
    # change characters to factors
    for (preds in predictors[characters])
      frame[[preds]] <- as.factor(frame[[preds]])
    factors <- factors | characters
    column.levels <- lapply(frame[factors], levels)
    
    # Now make them numeric
    for (preds in predictors[factors])
      frame[[preds]] <- as.numeric(frame[[preds]])
    x <- as.matrix(frame)
    attr(x, "column.levels") <- column.levels
  }
  else x <- as.matrix(frame[predictors])
  class(x) <- "rpartmv.matrix"
  x
}
