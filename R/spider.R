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
##    **Description of the spider data set**
##
## **************************************************************************
##
#'
#' Spider Data
#' 
#' Data set on abundances of spiders and environmental predictors. All variables
#' are rated on a 0-9 scale.
#' 
#' @docType data
#' 
#' @name spider
#' 
#' @usage data(spider)
#' 
#' @format A data frame with 28 observations with 12 species and six
#' environmental predictors.
#' 
#' @details Provided with the original package. The first 12 columns appears to
#' represent spider species since their names look like abbreviated binomial
#' names (for instance, the first name, arct.lute, might be Arctosa lutetiana,
#' whereas the second name, pard.lugu, might be Pardosa lugubris). The last six
#' columns appear to be abiotic (two of them: water, sand) and biotic (four of
#' them: moss, reft, twigs, herbs) environmental descriptors.
#' 
#' @source Van der Aart and Smeeck-Enserink 1975.
#' 
#' @references
#' Van der Aart, P. J. and N. Smeeck-Enserink. 1975. Correlations between
#' distributions of hunting spiders (Lycosidae, Ctenidae) and environmental
#' characteristics in a dune area. Netherlands Journal of Zoology 25: 1-45.
#' 
#' These data were analysed using multivariate trees in De'ath, G. 2002.
#' Multivariate Regression Trees: A New Technique for Modelling
#' Species-Environment Relationships. Ecology 83(4): 1103-1117
#' 
#' @examples ## Load the data set:
#' data("spider")
#' 
#' ## Number of rows and columns in the data table:
#' dim(spider)
#' 
#' ## First six rows of the data table:
#' head(spider)
#' 
#' ## Splitting the table into the responses (species, Y)
#' ## and the descriptors (X):
#' Y <- spider[,1L:12L]
#' X <- spider[,13L:18L]
#' 
#' ## Summary of the species data table:
#' summary(Y)
#' 
#' ## Summary of the descriptors data table:
#' summary(X)
#' 
#' ## Showing the species data using a principal component analysis:
#' prY <- princomp(Y)
#' 
#' ## Data table of the species loading for the first two principal components:
#' data.frame(
#'   AX1 = prY$loadings[,1L],
#'   AX2 = prY$loadings[,2L]
#' ) -> dat
#' 
#' par(mar=c(4.25,4.25,1.25,1.25))
#' plot(AX2 ~ AX1, data=dat, asp=1, xlim=c(-0.4,0.8), ylim=c(-0.4,0.8))
#' abline(v=0, lty=3L)
#' abline(h=0, lty=3L)
#' text(x=dat$AX1, y=dat$AX2 + 0.03, labels=colnames(Y), cex=0.75)
#' 
#' ## Descriptor: Water
#' par(mfrow=c(3L,2L), mar=c(4.25,4.25,1.25,1.25))
#' plot(arct.lute ~ water, data = spider, ylim=c(0,9), pch=21L, bg="red")
#' plot(pard.lugu ~ water, data = spider, ylim=c(0,9), pch=21L, bg="orange")
#' plot(zora.spin ~ water, data = spider, ylim=c(0,9), pch=21L, bg="yellow")
#' plot(aulo.albi ~ water, data = spider, ylim=c(0,9), pch=21L, bg="green")
#' plot(troc.terr ~ water, data = spider, ylim=c(0,9), pch=21L, bg="blue")
#' plot(alop.cune ~ water, data = spider, ylim=c(0,9), pch=21L, bg="purple")
#' 
#' ## Descriptor: Sand
#' par(mfrow=c(3L,2L), mar=c(4.25,4.25,1.25,1.25))
#' plot(arct.lute ~ sand, data = spider, ylim=c(0,9), pch=21L, bg="red")
#' plot(pard.lugu ~ sand, data = spider, ylim=c(0,9), pch=21L, bg="orange")
#' plot(zora.spin ~ sand, data = spider, ylim=c(0,9), pch=21L, bg="yellow")
#' plot(aulo.albi ~ sand, data = spider, ylim=c(0,9), pch=21L, bg="green")
#' plot(troc.terr ~ sand, data = spider, ylim=c(0,9), pch=21L, bg="blue")
#' plot(alop.cune ~ sand, data = spider, ylim=c(0,9), pch=21L, bg="purple")
#' 
#' @keywords spider
#' 
NULL
