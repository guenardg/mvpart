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
##    **Function gdist**
##
## **************************************************************************
##
#'
#' Dissimilarity Measures
#' 
#' The function computes useful dissimilarity indices which are known to have a
#' good rank-order relation with gradient separation and are thus efficient in
#' community ordination with multidimensional scaling.
#' 
#' @name gdist
#' 
#' @param x Data matrix.
#' @param method Dissimilarity index. Choices are: "manhattan", "euclidean",
#' "canberra", "bray", "kulczynski", "gower", "maximum", "binary", "chisq",
#' "chord", "beta0", "beta1", or "beta2"
#' @param keepdiag Compute amd keep diagonals.
#' @param full Return the square dissimilarity matrix.
#' @param sq Square the dissimilarities -- useful for distance-based
#' partitioning.
#' 
#' @details
#' Infamous ''double zeros'' are removed in Canberra dissimilarity.
#' 
#' Euclidean and Manhattan dissimilarities are not good in gradient separation
#' without proper standardization but are still included for comparison and
#' special needs.
#' 
#' Some of indices become identical or rank-order similar after some
#' standardization.
#' 
#' @returns Should be interchangeable with \code{\link[stats]{dist}} and returns
#' a distance object of the same type.
#' 
#' @author Jari Oksanen  -- modified Glenn De'ath (Dec 03)
#' 
#' @references
#' Faith, D.P, Minchin, P.R. and Belbin, L. (1987) Compositional dissimilarity
#' as a robust measure of ecological distance. \emph{Vegetatio} 69, 57-68.
#' 
#' @note The function is an alternative to \code{\link[stats]{dist}} adding some
#' ecologically meaningful indices.  Both methods should produce similar types
#' of objects which can be interchanged in any method accepting either.
#' Manhattan and Euclidean dissimilarities should be identical in both methods,
#' and Canberra dissimilarity may be similar.
#' 
## @examples
## data(spider)
## spider.dist <- gdist(spider[,1L:12L])
#' 
#' @keywords multivariate
#' 
#' @useDynLib mvpart, .registration = TRUE
#' 
#' @export
gdist <- function(x, method = "bray", keepdiag = FALSE , full = FALSE,
                  sq = FALSE) {
  METHODS <- c("manhattan","euclidean","canberra","bray","kulczynski","gower",
               "maximum","binary","chisq","chord","beta0","beta1","beta2")
  method <- pmatch(method, METHODS)
  if(is.na(method))
    stop("invalid distance method")
  N <- nrow(x <- as.matrix(x))
  if(method == 6L) x <- scaler(x, col=c("min0","max1"))
  if(method == 9L) {
    rr <- apply(x, 1L, sum)
    cc <- apply(x, 2L, sum)
    x <- diag(1/sqrt(rr)) %*% x %*% diag(1/sqrt(cc))
    method <- 2L
  }
  else if(method == 10L) {
    mns <- sqrt(apply(x^2, 1L, sum))
    x <- x/(mns*sqrt(2))
    method <- 2L
  }
  else if(method > 10L) method <- method - 2L
  .C("gdistance",
     x = as.double(x),
     nr = N,
     nc = ncol(x),
     d = double((N*(N - 1L))/2L),
     keepdiag = as.integer(FALSE),
     method = as.integer(method),
     PACKAGE = "mvpart")$d -> d
  attr(d, "Size") <- N
  class(d) <- "dist"
  if (full) d <- distfull(d)
  if (sq) d <- d^2
  d
}
