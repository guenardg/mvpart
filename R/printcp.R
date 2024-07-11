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
##    **Function printcp**
##
## **************************************************************************
##
#'
#' Displays CP table for an Rpartmv Object
#' 
#' Displays the \code{cp} table for fitted \code{rpartmv} object.
#' 
#' @name printcp
#' 
#' @param x An \code{\link{rpartmv-class}} object.
#' @param digits The number of digits of numbers to print.
#' 
#' @returns The \code{\link{rpartmv-class}} object's cp table.
#' 
#' @details Prints a table of optimal pruning based on a complexity parameter.
#'
#' @seealso \code{\link{summary.rpartmv}}
#' 
#' @keywords tree
#' 
#' @examples ## Load the car data set:
#' data("car.test.frame")
#' 
#' ## Estimating a regression tree with a single descriptor:
#' rt.auto1 <- rpartmv(Mileage ~ Weight, data=car.test.frame)
#' 
#' ## Print the model:
#' rt.auto1
#' 
#' ## Print the cp table:
#' printcp(rt.auto1)
#' 
#' @importFrom stats naprint
#' 
#' @export
printcp <- function(x, digits = getOption("digits") - 2L) {
  if (!inherits(x, 'rpartmv'))
    stop("Must be an rpartmv x")
  cat(
    switch(
      x$method,
      anova = "\nRegression tree:\n" ,
      class = "\nClassification tree:\n" ,
      poisson="\nRates regression tree:\n",
      exp = "\nSurvival regression tree:\n"
    )
  )
  if(!is.null(cl <- x$call)) {
    dput(cl)
    cat("\n")
  }
  frame <- x$frame
  leaves <- frame$var == "<leaf>"
  used <- unique(frame$var[!leaves])
  if(!is.null(used)) {
    cat("Variables actually used in tree construction:\n")
    print(sort(as.character(used)), quote=FALSE)
    cat("\n")
  }
  cat(
    "Root node error: ", format(frame$dev[1], digits=digits), '/',
    frame$n[1], ' = ',
    format(frame$dev[1]/frame$n[1], digits=digits),
    '\n\n', sep=''
  )
  n <- x$frame$n
  omit <- x$na.action
  if (length(omit)) {
    cat("n=", n[1], " (", naprint(omit), ")\n\n", sep="")
  } else cat("n=", n[1], "\n\n")
  print(x$cptable, digits=digits)
  invisible(x$cptable)
}
