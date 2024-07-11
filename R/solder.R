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
##    **Description of the solder data set**
##
## **************************************************************************
##
#'
#' Soldering of Components on Printed-Circuit Boards
#' 
#' The \code{solder} data frame has 720 rows and six columns, representing a
#' balanced subset of a designed experiment varying five factors on the
#' soldering of components on printed-circuit boards.
#' 
#' @docType data
#' 
#' @name solder
#' 
#' @usage data(solder)
#' 
#' @format This data frame contains the following columns:
#' \describe{
#'   \item{\code{Opening}}{
#'     a factor with levels
#'     \code{L} 
#'     \code{M} 
#'     \code{S}
#'     indicating the amount of clearance around the mounting pad.
#'   }
#'   \item{\code{Solder}}{
#'     a factor with levels
#'     \code{Thick} 
#'     \code{Thin}
#'     giving the thickness of the solder used.
#'   }
#'   \item{\code{Mask}}{
#'     a factor with levels
#'     \code{A1.5} 
#'     \code{A3} 
#'     \code{B3} 
#'     \code{B6}
#'     indicating the type and thickness of mask used.
#'   }
#'   \item{\code{PadType}}{
#'     a factor with levels
#'     \code{D4} 
#'     \code{D6} 
#'     \code{D7} 
#'     \code{L4} 
#'     \code{L6} 
#'     \code{L7} 
#'     \code{L8} 
#'     \code{L9} 
#'     \code{W4} 
#'     \code{W9}
#'     giving the size and geometry of the mounting pad.
#'   }
#'   \item{\code{Panel}}{
#'     \code{1:3} indicating the panel on a board being tested.
#'   }
#'   \item{\code{skips}}{
#'     a numeric vector giving the number of visible solder skips.
#'   }
#' }
#' 
#' @source Chambers and Hastie (1992).
#' 
#' @references
#' John M. Chambers and Trevor J. Hastie eds. (1992) Statistical Models in S,
#' Wadsworth and Brooks/Cole, Pacific Grove, CA 1992.
#' 
#' @examples ## Load the data set:
#' data("solder")
#' 
#' ## Number of rows and columns in the data table:
#' dim(solder)
#' 
#' ## First six rows of the data table:
#' head(solder)
#' 
#' ## Summary of the data table:
#' summary(solder)
#' 
#' ## A few plots of the data within the data frame:
#' plot(skips ~ Opening, data=solder, las=1L)
#' plot(skips ~ Solder, data=solder, las=1L)
#' plot(skips ~ Mask, data=solder, las=1L)
#' plot(skips ~ PadType, data=solder, las=1L)
#' plot(skips ~ as.factor(Panel), data=solder, las=1L, xlab="Panel number")
#' 
#' @keywords datasets
#' 
NULL
