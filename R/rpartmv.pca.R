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
##    **Function rpartmv.pca**
##
## **************************************************************************
##
#' 
#' Principal Components Plot of a Multivariate Rpartmv Object
#' 
#' Plots a PCA of the rpartmv object on the current graphics device.
#' 
#' @name rpartmv.pca
#' 
#' @param tree A fitted object of class \code{rpartmv} containing a multivariate
#' regression tree.
#' @param pts If \code{TRUE}, large points representing the leaf means are
#' plotted.
#' @param plt.allx If \code{TRUE}, small points representing individual cases
#' are plotted.
#' @param speclabs If \code{TRUE} the labels of the response variables are
#' plotted.
#' @param specvecs If \code{TRUE} the vectors of the response variables are
#' plotted provided \code{wgt.ave} is \code{FALSE}.
#' @param wgt.ave If \code{TRUE} use weighted averages of responses not vectors.
#' @param add.tree If \code{TRUE} add the tree structure to the plot.
#' @param cv1 Defines the principal component to plot horizontally -- but see
#' \code{interact}.
#' @param cv2 Defines the principal component to plot vertically -- but see
#' \code{interact}.
#' @param chulls If \code{TRUE} adds convex hulls to thr tree groups.
#' @param interact If \code{TRUE} the plot can be viewed in dimensions by
#' left-clicking to top-left, bottom-right or bottom-left (reset).
#' @param ... Arguments to be passed to or from other functions or methods.
#' 
#' @returns \code{NULL} (well, this is not really true).
#' 
#' @details This function plots a PCA biplot of the group means (leaves) of
#' multivariate regresssion objects of class \code{rpartmv}. The responses and
#' group means and indivdual cases can be shown on the plot. If responses are
#' positive (e.g., species-environment data) weighted averages of responses can
#' be plotted.
#' 
#' The PCA biplot plot is produced on the current graphics device.
#' 
#' @seealso \code{\link{rpartmv}}
#' 
#' @examples ## Load the spider data set:
#' data(spider)
#' 
#' ## Splitting the table into the responses (species, Y)
#' ## and the descriptors (X):
#' Y <- data.matrix(spider[,1L:12L])
#' X <- spider[,13L:18L]
#' ## Note: the multivariate response needs to procesed using function data.matrix
#' ## before being used as a multivariate response by function rpartmv.
#' 
#' ## Multivariate regression tree on the spider species response matrix and all
#' ## the descriptors in the spider data set.
#' 
#' ## Using the default settings:
#' spider1 <- mvpart(Y ~ herbs + reft + moss + sand + twigs + water, data=X)
#' 
#' ## Non-interactive PCA plot:
#' rpartmv.pca(spider1)
#' 
#' ## Interactive PCA plot:
#' rpartmv.pca(spider1, wgt.ave=TRUE, interact=TRUE)
#'
#' @keywords tree
#' 
#' @importFrom grDevices chull
#' @importFrom stats cor
#' 
#' @export
rpartmv.pca <- function(tree, pts = TRUE, plt.allx = TRUE, speclabs = TRUE,
                      specvecs = TRUE, wgt.ave = FALSE, add.tree = TRUE,
                      cv1 = 1L, cv2 = 2L, chulls = TRUE, interact = FALSE,
                    ...) {
  if(tree$method != "mrt")
    stop("Only for multivariate trees !! \n")
  if(nrow(tree$frame) < 4L)
    stop("Only 2 terminal nodes -- PCA not done !! \n")
  old.par <- par(mar=rep(2,4L), xpd=TRUE)
  on.exit(par(old.par))
  frame <- tree$frame
  ncf <- ncol(frame)
  data <- tree$y
  ny <- ncol(data)
  treegrps <- tree$where
  specs <- dimnames(data)[[2L]]
  leaves <- frame$var == "<leaf>"
  n <- length(leaves)
  ln <- sum(leaves)
  lnot <- sum(!leaves)
  key <- dimnames(frame)[[1L]]
  node <- as.numeric(key)
  even.node <- node[even <- node%%2 == 0L]
  num <- length(specs)
  node.means <- as.matrix(frame[,ncf])
  tnode.means <- node.means[leaves,]
  dimnames(node.means) <- list(key, specs)
  mat <- amat <- node.means - node.means[rep(1L,n),]
  mat <- mat[leaves,]
  temp <- mat[rep(1L:ln, frame[leaves,2L]),]
  z <- svd(temp)
  maxd <- sum(z$d > 1e-06)
  d <- diag(z$d[1L:maxd])
  xall <- z$u[,1L:maxd,drop = FALSE] %*% d
  x <- amat %*% (z$v)[,1L:maxd, drop=FALSE]
  xlv <- x[leaves,]
  if(!wgt.ave) {
    y <- z$v[,1L:maxd,drop=FALSE]
  } else {
    specvecs <- FALSE
    rc <- apply(tnode.means * frame$n[leaves], 2L, sum)
    wgt <- diag(1/rc) %*% t(tnode.means * frame$n[leaves])
    y <- wgt %*% xlv
  }
  label <- 4L:(3L + num)
  dstat <- signif(frame[leaves, "yval2"], digits=options()$digits)
  ln <- dim(dstat)[1L]
  stat <- vector("character", length = ln)
  for(i in 1L:ln) stat[i] <- paste(dstat[i,], collapse = ", ")
  ymax <- max(dstat)
  ymin <- min(0,min(dstat))
  treegrps <- as.numeric(factor(treegrps))
  (scale(
    as.matrix(data),
    center = TRUE,
    scale = FALSE) %*% z$v
  )[,1L:maxd,drop=FALSE] -> xx
  xrb <- rbind(x, xx)
  if(plt.allx) {
      mxx <- sqrt(apply(xrb[,c(cv1,cv2)]^2, 1L, sum))
  } else mxx <- sqrt(apply(x[,c(cv1,cv2)]^2, 1L, sum))
  cvar <- round((100 * z$d[1L:maxd]^2)/sum(z$d[1L:maxd]^2), digits = 2L)
  cvar2 <- round(diag(cor(xall, xx[order(tree$where),]))[1L:maxd], 3L)
  dlabs <- paste("   Dim ", c(1L:maxd), " ", cvar, "% : [", cvar2, "]")
  myy <- sqrt(apply(y[,c(cv1,cv2)]^2, 1L, sum))
  sc <- ifelse(wgt.ave, 1, max(mxx)/max(myy))
  repeat{
    plot(c(sc*y[,cv1],xx[,cv1]), c(sc*y[,cv2],xx[,cv2]), axes=FALSE,
         xlab="", ylab="", type="n", asp=1)
    cxy <- par("cxy")
    sze <- par()$fin/par()$din
    adj <- ifelse(pts, cxy[2L]*sze[2L], 0)
    if(specvecs)
      segments(sc*y[,cv1], sc*y[,cv2], rep(0,nrow(y)), rep(0,nrow(y)),
               col="gray", lty=1L)
    mtext(dlabs[cv1], side=1L, las=0, adj=0, line=0, cex=0.85*par()$cex)
    mtext(dlabs[cv2], side=2L, las=0, adj=0, line=0, cex=0.85*par()$cex)
    if(add.tree) {
      pp <- match(c(even.node,even.node + 1L), node)
      nn <- length(even.node)
      from <- pp[1L:nn]
      to <- pp[(nn + 1L):(2L*nn)]
      segments(x[from,cv1], x[from,cv2], x[to,cv1], x[to,cv2])
    }
    if(chulls) {
      unitg <- sort(unique(treegrps))
      for(i in 1L:length(unitg)) {
        hpts <- chull(xx[unitg[i] == treegrps,c(cv1,cv2)])
        hpts <- c(hpts,hpts[1])
        lines(xx[unitg[i] == treegrps, c(cv1,cv2)][hpts,], col=i + 1L)
      }
    }
    if(plt.allx) {
      unitg <- sort(unique(treegrps))
      for(i in 1L:length(unitg))
        points(xx[unitg[i] == treegrps,cv1], xx[unitg[i] == treegrps,cv2],
               pch=21L, col=1L, bg=i + 1L, cex=1.2*par()$cex)
    }
    if(pts) {
      lvnode <- sort(node[leaves])
      for(i in 1L:length(lvnode))
        points(xlv[,cv1][lvnode[i] == lvnode], xlv[,cv2][lvnode[i] == lvnode],
               pch=21L, cex=2*par()$cex, col=1L, bg=i + 1L)
    }
    if(speclabs)
      text(sc*y[,cv1], sc*y[,cv2] + 0.5*adj*specvecs*(y[,cv2] > 0),
           specs, col="black", cex=par()$cex)
    points(0, 0, pch=3L, cex=par()$cex*2.5, col=1)
    if(interact) {
      z <- locator(1L)
      if(length(z$x)) {
        if (z$x > 0 & z$y < 0) {
          if(cv1 < maxd) {
            cv1 <- cv1 + 1L
          } else cv1 <- 1L
        } else if(z$x < 0 & z$y > 0) {
          if(cv2 < maxd) {
            cv2 <- cv2 + 1L
          } else cv2 <- 2L
        } else if(z$x < 0 & z$y < 0) {
          cv1 <- 1L
          cv2 <- 2L
        }
      } else break
    } else break
  }
  invisible(list(y = sc * y, xlv = xlv, xx = xx))
}
