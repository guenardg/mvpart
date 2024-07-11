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
##    **Function snip.rpartmv**
##
## **************************************************************************
##
#'
#' Snip Subtrees of an Rpartmv Object
#' 
#' Creates a "snipped" \code{\link{rpartmv-class}} object, containing the nodes
#' that remain after selected sub trees have been snipped off. The user can snip
#' nodes using the toss argument, or interactively by clicking the mouse button
#' on specified nodes within the graphics window.
#' 
#' @name snip.rpartmv
#' 
#' @param x An \code{\link{rpartmv-class}} object.
#' @param toss An integer vector containing indices (node numbers) of all
#' sub trees to be snipped off. If missing, user selects branches to snip off as
#' described below.
#' 
#' @returns A \code{\link{rpartmv-class}} object containing the nodes that
#' remain after specified or selected sub trees have been snipped off.
#' 
#' @details A dendrogram of \code{rpartmv} is expected to be visible on the
#' graphics device, and a graphics input device (e.g., a mouse) is required.
#' Clicking (the selection button) on a node displays the node number, sample
#' size, response yvalue, and Error (dev). Clicking a second time on the same
#' node snips that sub tree off and visually erases the sub tree. This process
#' may be repeated an number of times. Warnings result from selecting the root
#' or leaf nodes.  Clicking the exit button will stop the snipping process and
#' return the resulting \code{\link{rpartmv-class}} object.
#' 
#' Warning: visually erasing the plot is done by over-plotting with the
#' background colour. This will do nothing if the background is transparent
#' (often true for screen devices).
#' 
#' See the documentation for the specific graphics device for details on
#' graphical input techniques.
#' 
#' @seealso \code{\link{plot.rpartmv}}
#' 
#' @keywords tree
#' 
#' @export
snip.rpartmv <- function(x, toss) {
  if(!inherits(x,"rpartmv"))
    stop("Not an rpartmv object")
  if(missing(toss) || length(toss)==0) {
    toss <- snip.rpartmv.mouse(x)
    if (length(toss)==0) return(x)
  }
  where <- x$where
  ff   <- x$frame
  id    <- as.numeric(row.names(ff))
  index <- ff$index
  ff.n  <- length(id)
  toss <- unique(toss)
  toss.idx <- match(toss, id, nomatch=0) #the rows of the named nodes
  if(any(toss.idx ==0)) {
    warning(paste("Nodes", toss[toss.idx==0], "are not in this tree"))
    toss <- toss[toss.idx>0]
    toss.idx <- toss.idx[toss.idx>0]
  }
  #    if (any(toss==1))  {
  #   # a special case that causes grief later
  #   warning("Can't prune away the root node and still have a tree!")
  #        return(NULL)
  #   }
  # Now add all of the descendants of the selected nodes
  #   We do this be finding all node's parents.
  #        (Division by 2 gives the parent of any node.)
  #   At each step we make id2 <- parent(id2), and augment 'toss' with
  #     found children.  The loop should take <  log_2(maxdepth)/2 steps
  id2 <- id
  while(any(id2>1)) {
    id2 <- floor(id2/2)
    xx <- (match(id2, toss, nomatch=0) >0)
    toss <- c(toss, id[xx])
    id2[xx] <- 0
  }
  # Now "toss" contains all of the nodes that should not be splits
  temp <- match(floor(toss/2) , toss, nomatch=0)  #which are leaves?
  newleaf <- match(toss[temp==0], id)             # row numbers, leaves
  keepit <- (1:ff.n)[is.na(match(id,toss))]  # row numbers to be let be
  # Compute the parent row for each row in the splits structure
  #  Then "thin out" the splits and csplit components
  n.split <- rep((1:ff.n), ff$ncompete + ff$nsurrogate+ 1*(ff$var!='<leaf>'))
  split <- x$splits[match(n.split, keepit, nomatch=0) >0, ,drop=FALSE]
  temp <- split[,2] >1      #which rows point to categoricals?
  if(any(temp)) {
    x$csplit <- x$csplit[split[temp,4], , drop=FALSE]
    split[temp,4] <- 1
    if(is.matrix(x$csplit)) split[temp,4] <- 1:nrow(x$csplit)
  } else x$csplit <- NULL
  x$splits <- split
  # Thin out unneeded rows in the frame component
  ff$ncompete[newleaf] <- ff$nsurrogate[newleaf] <- 0
  ff$var[newleaf]     <- "<leaf>"
  x$frame <- ff[sort(c(keepit, newleaf)),]
  # Now do the 'parents' loop one more time, to fix up the "where"
  #   vector
  # This pass requires log_2(depth) iterations
  #
  id2 <- id[x$where]         #the list of old leaf nodes
  id3 <- id[sort(c(keepit, newleaf))]
  temp <- match(id2, id3, nomatch=0)
  while (any(temp==0)) {
    id2[temp==0] <- floor(id2[temp==0]/2)
    temp <- match(id2, id3, nomatch=0)
  }
  x$where <- match(id2, id3)
  x
}
