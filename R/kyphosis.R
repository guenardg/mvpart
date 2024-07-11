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
##    **Description of the kyphosis data set**
##
## **************************************************************************
##
#'
#' Data on Children who have had Corrective Spinal Surgery
#' 
#' The \code{kyphosis} data frame has 81 rows and four columns representing data
#' on children who have had corrective spinal surgery.
#' 
#' @docType data
#' 
#' @name kyphosis
#' 
#' @usage data(kyphosis)
#' 
#' @format This data frame contains the following columns:
#' \describe{
#'   \item{\code{Kyphosis}}{
#'     a factor with levels
#'     \code{absent} 
#'     \code{present}
#'     indicating if a kyphosis (a type of deformation)
#'     was present after the operation.
#'   }
#'   \item{\code{Age}}{
#'     in months
#'   }
#'   \item{\code{Number}}{
#'     the number of vertebrae involved 
#'   }
#'   \item{\code{Start}}{
#'     the number of the first (topmost) vertebra operated on.
#'   }
#' }
#' 
#' @details Provided with the original package
#' 
#' @source Chambers and Hastie (1992).
#' 
#' @references
#'  John M. Chambers and Trevor J. Hastie eds. (1992) Statistical Models in S,
#'  Wadsworth and Brooks/Cole, Pacific Grove, CA 1992.
#' 
#' @examples ## Load the data set:
#' data("kyphosis")
#' 
#' ## Number of rows and columns in the data table:
#' dim(kyphosis)
#' 
#' ## First six rows of the data table:
#' head(kyphosis)
#' 
#' ## Summary of the data table:
#' summary(kyphosis)
#' 
#' ## A few plots of the data within the data frame:
#' plot(Number ~ Age, data=kyphosis, las=1L)
#' plot(Start ~ as.factor(Number), data=kyphosis, las=1L)
#' plot(Age/12 ~ Kyphosis, data=kyphosis, las=1L, ylab="Age (years)")
#' plot(Number ~ Kyphosis, data=kyphosis, las=1L)
#' plot(Age/12 ~ Start, data=kyphosis, las=1L, ylab="Age (years)")
#' 
#' @keywords datasets
#' 
NULL
