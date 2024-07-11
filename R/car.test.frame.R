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
##    **Description of the car.test.frame data set**
##
## **************************************************************************
##
#'
#' Automobile Data from `Consumer Reports' 1990
#' 
#' The \code{car.test.frame} data frame has 60 rows and eight columns, giving
#' data on makes of cars taken from the April, 1990 issue of magazine Consumer
#' Reports.
#' 
#' @docType data
#' 
#' @name car.test.frame
#' 
#' @usage data(car.test.frame)
#' 
#' @format This data frame contains the following columns:
#' \describe{
#'   \item{\code{Price}}{
#'     a numeric vector giving the list price in US dollars of a standard model
#'   }
#'   \item{\code{Country}}{
#'     of origin, a factor with levels
#'     \code{France} 
#'     \code{Germany} 
#'     \code{Japan} 
#'     \code{Japan/USA} 
#'     \code{Korea} 
#'     \code{Mexico} 
#'     \code{Sweden} 
#'     \code{USA} 
#'   }
#'   \item{\code{Reliability}}{
#'     a numeric vector coded \code{1} to \code{5}.
#'   }
#'   \item{\code{Mileage}}{
#'     fuel consumption miles per US gallon, as tested. 
#'   }
#'   \item{\code{Type}}{
#'     a factor with levels
#'     \code{Compact} 
#'     \code{Large} 
#'     \code{Medium} 
#'     \code{Small} 
#'     \code{Sporty} 
#'     \code{Van} 
#'   }
#'   \item{\code{Weight}}{
#'     kerb weight in pounds.
#'   }
#'   \item{\code{Disp.}}{
#'     the engine capacity (displacement) in litres.
#'   }
#'   \item{\code{HP}}{
#'     the net horsepower of the vehicle.
#'   }
#' }
#' 
#' @source Consumer Reports, April, 1990, pp. 235--288 quoted in Chambers and
#' Hastie (1992).
#' 
#' @references
#' John M. Chambers and Trevor J. Hastie eds. (1992) Statistical Models in S,
#' Wadsworth and Brooks/Cole, Pacific Grove, CA 1992, pp. 46--47.
#' 
#' @examples ## Load the data set:
#' data("car.test.frame")
#' 
#' ## Number of rows and columns in the data table:
#' dim(car.test.frame)
#' 
#' ## First six rows of the data table:
#' head(car.test.frame)
#' 
#' ## Summary of the data table:
#' summary(car.test.frame)
#' 
#' ## A few plots of the data within the data frame:
#' plot(Mileage ~ Weight, data=car.test.frame, las=1L)
#' plot(Disp. ~ HP, data=car.test.frame, log="xy", las=1L)
#' plot(Price/1000 ~ Type, data=car.test.frame, las=1L,
#'      ylab="Price (1000$)", ylim=c(5,30))
#' 
#' @keywords datasets car
#' 
NULL
