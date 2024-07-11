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
##    **Declaration of user-defined generics**
##
## **************************************************************************
##
#' 
#' Generic Methods
#' 
#' A set of non-standard generic functions used by package mvpart.
#' 
#' @name rpartmv-generics
#' 
#' @param object An object of a class for which the method is implemented
#' (e.g., \code{rpartmv-class} object).
#' @param ... Optional arguments to be passed internally to other functions.
#' 
#' @return The return value depends on the implementation of the methods for a
#' particular S3 class.
#' 
#' @seealso \code{prune.rpartmv} \code{rsq.rpartmv}, and \code{meanvar.rpartmv}
#' 
NULL
#' 
#' @describeIn rpartmv-generics
#' 
#' Prune Method
#' 
#' Method to prune branch for a regression or classification tree
#' 
#' @export
prune <- function(object, ...)  UseMethod("prune")
#' 
#' @describeIn rpartmv-generics
#' 
#' R-square Method
#' 
#' Method to generate the R-square variation plot of an object.
#' 
#' @export
rsq <- function(object, ...) UseMethod("rsq")
#' 
#' @describeIn rpartmv-generics
#' 
#' Meanvar method
#' 
#' Method to generate the mean-variance plot for an object.
#' 
#' @export
meanvar <- function(object, ...) UseMethod("meanvar")
