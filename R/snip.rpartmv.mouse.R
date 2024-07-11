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
##    **Internal function snip.rpartmv.mouse**
##
## **************************************************************************
##

#' @importFrom graphics identify lines

"snip.rpartmv.mouse" <- function(tree, parms) {
  #~     parms=paste(".rpartmv.parms", dev.cur(), sep = ".")) {
  xy <- rpartmvco(tree,parms)
  uniform <- parms$uniform
  nspace <- parms$nspace
  nbranch <- parms$nbranch
  minbranch <- parms$minbranch
  toss <- NULL
  ff <- tree$frame
  if(exists(parms, envir=.GlobalEnv)) {
    parms <- get(parms, envir=.GlobalEnv)
    branch <- parms$branch
  }
  else branch <- 1
  node <- as.numeric(row.names(tree$frame))
  draw <- rpartmv.branch(xy$x,xy$y, node, branch)
  lastchoice <- 0
  while(length(choose <- identify(xy, n=1, plot=FALSE)) > 0) {
    if(ff$var[choose] == '<leaf>') {
      cat("Terminal node -- try again\n")
      next
    }
    if(choose != lastchoice) {
      # print out some info on the click
      cat("node number:", node[choose], " n=", ff$n[choose], "\n")
      cat("    response=", format(ff$yval[choose]))
      if(is.null(ff$yval2)) cat ("\n")
      else if(is.matrix(ff$yval2))
        cat(" (", format(ff$yval2[choose,]), ")\n")
      else cat(" (", format(ff$yval2[choose]), ")\n")
      cat("    Error (dev) = ", format(ff$dev[choose]), "\n")
      lastchoice <- choose
    } else {
      # second click-- erase all of the descendants
      #   (stolen from snip.tree)
      id  <- node[choose]
      id2 <- node
      while(any(id2 > 1)) {
        id2 <- floor(id2/2)
        temp  <- (match(id2, id, nomatch=0) > 0)
        id <- c(id, node[temp])
        id2[temp] <- 0
      }
      temp <- match(id, node[ff$var != '<leaf>'], nomatch=0)
      lines(c(draw$x[,temp]), c(draw$y[,temp]), col=0)
      toss <- c(toss, node[choose])
    }
  }
  toss
}
