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
##    **Function xdiss**
##
## **************************************************************************
##
#'
#' Extended Dissimilarity Measures
#' 
#' The function computes extended dissimilarity indices which are for long
#' gradients have better good rank-order relation with gradient separation and
#' are thus efficient in community ordination with multidimensional scaling.
#' 
#' @name xdiss
#' 
#' @param data Data matrix.
#' @param dcrit Dissimilarities < \code{dcrit} are considered to have no species
#' in common and are recalculated.
#' @param dauto Automatically select tuning parameters -- recommended.
#' @param method Dissimilarity index.
#' @param use.min Minimum dissimilarity of pairs of distances used --
#' recommended.
#' @param dinf Internal parameter -- leave as is usually.
#' @param eps Internal parameter -- leave as is usually.
#' @param replace.neg Internal parameter -- leave as is usually.
#' @param big Internal parameter -- leave as is usually.
#' @param sumry Whether to print summary of extended dissimilarities.
#' @param full Return the square dissimilarity matrix.
#' @param sq Square the dissimilarities -- useful for distance-based
#' partitioning.
#' 
#' @details The function knows the same dissimilarity indices as
#' \code{\link{gdist}}.
#' 
#' @returns Returns an object of class distance with attributes "Size" and "ok".
#' "ok" is TRUE if rows are not disconnected (De'ath 1999).
#' 
#' @author Glenn De'ath
#' 
#' @references
#' De'ath, G. (1999) Extended dissimilarity: a method of robust estimation of
#' ecological distances from high beta diversity data. Plant Ecology 144(2):
#' 191-199.
#' 
#' Faith, D.P, Minchin, P.R. and Belbin, L. (1987) Compositional dissimilarity
#' as a robust measure of ecological distance. Vegetatio 69: 57-68.
#' 
#' @examples  ## Load the spider data set:
#' data(spider)
#' 
#' ## Getting the table of the first six species:
#' Y <- data.matrix(spider[1L:6L,1L:12L])
#' 
#' ## Calculate the dissimilarity measures:
#' xdiss(Y, method="manhattan")  ## Default
#' xdiss(Y, method="euclidean")
#' xdiss(Y, method="canberra")
#' xdiss(Y, method="bray")
#' xdiss(Y, method="kulczynski")
#' xdiss(Y, method="gower")
#' xdiss(Y, method="maximum")
#' xdiss(Y, method="binary")
#' xdiss(Y, method="chisq")
#' xdiss(Y, method="chord")
#' 
#' @keywords multivariate
#' 
#' @useDynLib mvpart, .registration = TRUE
#' 
#' @export
xdiss <- function(data, dcrit = 1, dauto = TRUE , dinf = 0.5, method = "man",
                  use.min = TRUE, eps = 0.0001, replace.neg = TRUE, big = 10000,
                  sumry = TRUE, full = FALSE, sq = FALSE) {
  scale.row <- function(data, p = 1) {
    tmp <- apply(data, 1L, sum, na.rm = TRUE)
    if(any(t0 <- (tmp == 0)))
      cat(sum(t0), " rows with sum = 0 !!  -- these rows untransformed\n")
    if(p == 1) {
      data[!t0,] <- data[!t0,]/apply(data[!t0,], 1L, sum, na.rm=TRUE)
    } else if(p == 2)
      data[!t0,] <- data[!t0,]/(apply(data[!t0,]^2, 1L, sum, na.rm=TRUE))^0.5
    data
  }
  METHODS <- c("manhattan","euclidean","canberra","bray","kulczynski","gower", 
               "maximum","binary","chisq","chord")
  method <- METHODS[pmatch(method, METHODS)]
  if(is.na(method))
    stop("invalid distance method")
  if(any(data < 0))
    data <- apply(data, 2L, function(x) x - min(x))
  n <- dim(data)[1L]
  if(method == "chisq" | method == "gower" | method == "maximum") {
    method <- "manhattan"
    cat("This dissimilarity is not suitable -- switching to Manhattan\n")
  }
  if(method == "manhattan") {
    cat("Using Extended Dissimilarity : Manhattan (Site Standardised by Mean)\n")
    data <- scale.row(as.matrix(data), p=1)/2
    d <- gdist(data, method="man")
  }
  else if(method == "chord") {
    cat("Using Extended Dissimilarity : Chord \n")
    data <- scale.row(as.matrix(data), p=2)/sqrt(2)
    d <- gdist(data, method="euc")
  }
  else if(method == "euclidean") {
    cat("Using Extended Dissimilarity : Euc (Site Standardised by SS) \n")
    data <- scale.row(as.matrix(data), p=2)/sqrt(2)
    d <- gdist(data, method="euc")
  }
  else if(method == "bray") {
    cat("Using Extended Dissimilarity : Bray \n")
    d <- gdist(data, method="bra")
  }
  else if(method == "canberra") {
    cat("Using Extended Dissimilarity : Canberra \n")
    d <- gdist(data, method="can")
  }
  else if(method == "binary") {
    cat("Using Extended Dissimilarity : Binary \n")
    d <- gdist(data, method="bin")
  }
  else if(method == "kulczynski") {
    cat("Using Extended Dissimilarity : Kulczynski \n")
    d <- gdist(data, method="kul")
  }
  cat("Maximum distance = ", round(max(d), 4), "\n")
  if(dauto) {
    dcrit <- max(apply(distfull(d) + diag(1, n), 1L, min))
    dcrit <- dcrit * (1 - dinf) + dinf
    cat("Critical distance = ", signif(dcrit, 4), "\n")
  }
  cat("% Distances > Crit Dist = ", round(100*mean(d > dcrit), 2), "\n")
  use.min <- ifelse(use.min, 1, 0)
  storage.mode(d) <- "double"
  storage.mode(n) <- "integer"
  storage.mode(dcrit) <- "double"
  storage.mode(use.min) <- "integer"
  storage.mode(eps) <- "double"
  storage.mode(big) <- "double"
  .C("xdists",
     d = d,
     n,
     dcrit,
     use.min,
     eps,
     big,
     PACKAGE="mvpart")$d -> dnew
  if(any(dnew == -1, na.rm = TRUE )) {
    attr(dnew, "ok") <- FALSE
  } else attr(dnew, "ok") <- TRUE 
  if(any(dnew == -1))
    cat("WARNING : Data disconnected\n")
  if(replace.neg)
    dnew[dnew == -1] <- max(dnew)
  if(sumry) {
    cat("Summary of Extended Dissimilarities\n")
    print(summary(dnew))
  }
  attr(dnew, "Size") <- n
  class(dnew) <- "dist"
  if (full) dnew <- distfull(dnew)
  if (sq) dnew <- dnew^2
  dnew
}
