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
##    **Function scaler**
##
## **************************************************************************
##
#'
#' Row and Column Scaling of a Data Matrix
#' 
#' The function provides some popular (and effective) standardization methods
#' for community ecologists.
#' 
#' @name scaler
#' 
#' @param x Data matrix.
#' @param col Character vector of column standardization.
#' @param row Character vector of row standardization.
#' 
#' @details The function offers the following data matrix standardization:
#' \itemize{
#'   \item \code{mean1}: scale to mean of 1.
#'   \item \code{max1}: scale to maximum of 1.
#'   \item \code{ssq1}: scale to sums of squares equal 1.
#'   \item \code{range01}: scale range to 0-1.
#'   \item \code{zsc}: standardize to z-scores (mean=0, sd=1).
#'   \item \code{pa}: scale to presence/absence scale (0/1).
#'   \item \code{rank}: scale to rank order (1=lowest).
#' }
#' 
#' Standardization is performed first on columns then on rows. "pa" applies to
#' the whole matrix and can be specified using row or col.
#' 
#' @returns Returns the standardized matrix.
#' 
#' @author Jari Oksanen -- modified Glenn De'ath (Dec 03)
#' 
#' @note Common transformations can be made with standard R functions.
#' 
## @examples
## data(spider)
## spid.data <- scaler(spider, col = "max", row="mean1")
#' 
#' @keywords multivariate manip
#' 
#' @export
scaler <- function(
    x, col = c("mean1","max1","min0","ssq1","range01","zsc","pa","rank")[2L],
    row = c("mean1","max1","min0","ssq1","range01","zsc","pa","rank")[1L]) {
  fun <- function(x, method, MARGIN) {
    x <- as.matrix(x)
    switch(
      method,
      mean1 = {
        tmp <- apply(x, MARGIN, mean, na.rm = TRUE)
        x <- sweep(x, MARGIN, tmp, "/")
      }, max1 = {
        tmp <- apply(x, MARGIN, max, na.rm = TRUE)
        x <- sweep(x, MARGIN, tmp, "/")
      }, min0 = {
        tmp <- apply(x, MARGIN, min, na.rm = TRUE)
        x <- sweep(x, MARGIN, tmp, "-")
      }, ssq1 = {
        tmp <- apply(x^2, MARGIN, sum, na.rm = TRUE)
        tmp <- sqrt(tmp)
        x <- sweep(x, MARGIN, tmp, "/")
      }, range01 = {
        tmp <- apply(x, MARGIN, min, na.rm = TRUE)
        ran <- apply(x, MARGIN, max, na.rm = TRUE)
        ran <- ran - tmp
        x <- sweep(x, MARGIN, tmp, "-")
        x <- sweep(x, MARGIN, ran, "/")
      }, zsc = {
        if(MARGIN == 1L) x <- t(scale(t(x))) else x <- scale(x)
      }, pa = {
        tmp <- dim(x)
        nam <- dimnames(x)
        x <- as.numeric(x > 0)
        dim(x) <- tmp
        dimnames(x) <- nam
      }, rank = {
        x <- apply(x, MARGIN, rank)
        if (MARGIN == 1) 
          x <- t(x)
      }
    )
    x
  }
  METHODS <- c("mean1","max1","min0","ssq1","range01","zsc","pa","rank")
  if (is.null(col) & is.null(row)) 
    cat("Scalings are",METHODS,"\n")
  if(!is.null(col)) {
    for (i in 1:length(col)){
      method <- match.arg(col[i], METHODS)
      x <- fun(x, method = method, MARGIN = 2L)
    }
  }
  if(!is.null(row)) {
    for(i in 1:length(row)) {
      method <- match.arg(row[i], METHODS)
      x <- fun(x, method = method, MARGIN = 1L)
    }
  }
  x
}
