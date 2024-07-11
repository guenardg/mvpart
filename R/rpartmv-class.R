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
##    **Class rpartmv: definition and methods**
##
## **************************************************************************
##
#' 
#' Class and Methods for Recursive Partitioning and Regression Trees
#' 
#' Class and methods to handle an \code{rpartmv} model.
#' 
#' @docType class
#' 
#' @name rpartmv-class
#' @aliases labels.rpartmv print.rpartmv prune.rpartmv predict.rpartmv
#' residuals.rpartmv rsq.rpartmv
#' 
#' @param object An object of class \code{rpartmv}. This is assumed to be the
#' result of some function that produces an object with the same named
#' components as that returned by the \code{rpartmv} function.
#' @param x Identical as argument \code{object}.
#' @param digits the number of digits to be used for numeric values. All of the
#' \code{rpartmv} functions that call labels explicitly set this value, with
#' \code{options("digits")} as the default.
#' @param minlength The minimum length for abbreviation of character or factor
#' variables. If 0 no abbreviation is done; if 1 then single letters are used
#' with \code{"a"} for the first level, \code{"b"} for the second and so on. If
#' the value is greater than 1, the \code{abbreviate} function is used.
#' @param pretty An argument included for backwards compatibility:
#' \code{pretty=0} implies \code{minlength=0}, \code{pretty=NULL} implies
#' \code{minlength=1}, and \code{pretty=TRUE} implies \code{minlength=4}. For
#' \code{text.rpartmv}: an integer denoting the extent to which factor levels in
#' split labels will be abbreviated.  A value of (0) signifies no abbreviation.
#' A \code{NULL}, the default, signifies using elements of letters to represent
#' the different factor levels.
#' @param collapse Whether the returned set of labels is always of the same
#' length as the number of nodes in the tree (logical). If \code{collapse=TRUE}
#' (default), the returned value is a vector of labels for the branch leading
#' into each node, with \code{"root"} as the label for the top node. If
#' \code{FALSE}, the returned value is a two column matrix of labels for the
#' left and right branches leading out from each node, with \code{"leaf"} as the
#' branch labels for terminal nodes.
#' @param spaces The number of spaces to indent nodes of increasing depth.
#' @param cp Prune all nodes with a complexity less than \code{cp} from the
#' printout. Ignored if unspecified.
#' @param uniform If \code{TRUE}, uniform vertical spacing of the nodes is used;
#' this may be less cluttered when fitting a large plot onto a page. The default
#' is to use a non-uniform spacing proportional to the error in the fit.
#' @param branch Controls the shape of the branches from parent to child node.
#' Any number from 0 to 1 is allowed.  A value of 1 (the default) gives square
#' shouldered branches, a value of 0 give V shaped branches, with other values
#' being intermediate.
#' @param compress If \code{FALSE}, the leaf nodes will be at the horizontal
#' plot coordinates of \code{1:nleaves}. If \code{TRUE}, the routine attempts a
#' more compact arrangement of the tree. The compaction algorithm assumes
#' \code{uniform = TRUE}; surprisingly, the result is usually an improvement
#' even when that is not the case.
#' @param nspace The amount of extra space between a node with children and a
#' leaf, as compared to the minimal space between leaves. Applies to compressed
#' trees only. The default is the value of \code{branch}.
#' @param margin An extra percentage of white space to leave around the borders
#' of the tree. (Long labels sometimes get cut off by the default computation).
#' @param minbranch Set the minimum length for a branch to \code{minbranch}
#' times the average branch length. This parameter is ignored if
#' \code{uniform=TRUE}. Sometimes a split will give very little improvement, or
#' even (in the classification case) no improvement at all. A tree with branch
#' lengths strictly proportional to improvement leaves no room to squeeze in
#' node labels.
#' @param bar Length of bar at root (default = 0.03) -- used instead of char "|"
#' @param file Write the output to a given file name (Full listings of a tree
#' are often quite long).
#' @param newdata A data frame containing the values at which predictions are
#' required. The predictors referred to in the right side of
#' \code{formula(object)} must be present by name in \code{newdata}. If missing,
#' the fitted values are returned.
#' @param type A character string denoting the type of predicted or residual
#' values returned. For predictions, if the \code{rpartmv} object is a
#' classification tree, then the default is to return \code{prob} predictions,
#' a matrix whose columns are the probability of the first, second, etc. class.
#' (This agrees with the default behavior of tree).  Otherwise, a vector of the
#' results is returned. For residuals, in the case of a regression or anova
#' tree, all three residual definitions reduce to \code{y - fitted}. This is the
#' residual returned for \code{user} method trees as well. For classification
#' trees the \code{usual} residuals are the missclassification losses L(actual,
#' predicted) where L is the loss matrix.  With default losses this residual is
#' 0/1 for correct/incorrect classification. The \code{pearson} residual is
#' (1 - fitted)/sqrt(fitted(1 - fitted)) and the \code{deviance} residual is
#' sqrt(minus twice logarithm of fitted). For \code{poisson} and \code{exp} (or
#' survival) trees, the \code{usual} residual is the observed - expected number
#' of events. The \code{pearson} and \code{deviance} residuals are as defined in
#' McCullagh and Nelder.
#' @param splits logical flag. If \code{TRUE} (default), then the splits in the
#' tree are labelled with the criterion for the split.
#' @param which labels splits 1 = centre, 2 = left, 3 = right , 4 = both.
#' @param label A column name of \code{x$frame};  values of this will label the
#' nodes.  For the \code{"class"} method, \code{label="yval"} results in the
#' factor levels being used, \code{"yprob"} results in the probability of the
#' winning factor level being used, and 'specific yval level' results in the
#' probability of that factor level.
#' @param FUN The name of a labelling function, e.g. \code{text}.
#' @param all.leaves Logical. If \code{TRUE}, all nodes are labelled, otherwise
#' just terminal nodes.
#' @param tadj Adjustment of text above (or below) splits.
#' @param stats If \code{TRUE} adds statistics to nodes.
#' @param use.n If \code{TRUE} adds N to labels, i.e., the number of events in
#' level 1, number of events in level 2, etc., for \code{class}, \code{n} for
#' \code{anova}, and number of events / \code{n} for \code{poisson} and
#' \code{exp}.
#' @param bars If \code{TRUE} adds barplots for multivariate regression trees.
#' @param legend If \code{TRUE} adds legends for multivariate regression trees.
#' @param xadj varies the horizontal size of barplots for multivariate
#' regression trees.
#' @param yadj varies the vertical size of barplots for multivariate
#' regression trees.
#' @param bord Adds borders (boxes) to barplots for multivariate regression
#' trees.
#' @param big.pts Adds color coded points to nodes. Useful to track groups to
#' PCA plot (see \code{rpartmv.pca}).
#' @param uniform Uniform spacing of tree branches; default is \code{FALSE}.
#' @param xlab The x-axis label for the plot.
#' @param ylab The y-axis label for the plot.
#' @param ... Optional arguments to be passed internally to other functions.
#' 
#' @return The returned value depends on the method:
#' \describe{
#'   \item{labels.rpartmv}{ A vector of split labels (\code{collapse=TRUE}) or
#'   matrix of left and right splits (\code{collapse=FALSE}) for the supplied
#'   \code{rpartmv} object.  This function is called by printing methods for
#'   \code{rpartmv} and is not intended to be called directly by the users. }
#'   \item{plot.rpartmv}{The coordinates of the nodes are returned as a list,
#'   with components \code{x} and \code{y}.}
#'   \item{predict.rpartmv}{A new object is obtained by dropping \code{newdata}
#'   down the object. For factor predictors, if an observation contains a level
#'   not used to grow the tree, it is left at the deepest possible node and
#'   \code{frame$yval} at the node is the prediction.
#'   
#'   If \code{type="vector"}:\cr
#'   a vector of predicted responses. For regression trees this is the mean
#'   response at the node, for Poisson trees it is the estimated response rate,
#'   and for classification trees it is the predicted class.
#'   
#'   If \code{type="prob"}:\cr
#'   (for a classification tree) a matrix of class probabilities.
#'   
#'   If \code{type="matrix"}:\cr
#'   a matrix of the full responses (\code{frame$yval2} if this exists,
#'   otherwise \code{frame$yval}). For regression trees, this is the mean
#'   response, for Poisson trees it is the response rate and the number of
#'   events at that node in the fitted tree, and for classification trees it is
#'   the concatonation of the predicted class, the class counts at that node in
#'   the fitted tree, and the class probabilities.
#'   
#'   If \code{type="class"}:\cr
#'   (for a classification tree) a factor of classifications based on the
#'   responses.
#'   }
#'   \item{fitted.rpartmv}{A vector of residuals of type \code{type} from an
#'   \code{rpartmv} object.}
#'   \item{text.rpartmv}{\code{NULL} (invisibly).}
#'   \item{meanvar.rpartmv}{An invisible list containing the following vectors: x,
#'   fitted value at terminal nodes (\code{yval}); y, deviance of node divided
#'   by number of observations at node; and label, the node numbers.}
#' }
#' 
#' @details
#' \describe{
#'   \item{print.rpartmv}{This function is a method for the generic function
#'   \code{print} for class \code{"rpartmv"}.  It can be invoked by calling print
#'   for an object of the appropriate class, or directly by calling
#'   \code{print.rpartmv} regardless of the class of the object. Side effects: A
#'   semi-graphical layout of the contents of \code{x$frame} is printed.
#'   Indentation is used to convey the tree topology. Information for each node
#'   includes the node number, split, size, deviance, and fitted value.  For the
#'   \code{"class"} method, the class probabilties are also printed.}
#'   \item{plot.rpartmv}{This function is a method for the generic function
#'   \code{plot}, for objects of class \code{rpartmv}. The y-coordinate of the top
#'   node of the tree will always be 1. Side effects: an unlabeled plot is
#'   produced on the current graphics device.}
#'   \item{summary.rpartmv}{This function is a method for the generic function
#'   summary for class \code{"rpartmv"}.  It can be invoked by calling
#'   \code{summary} for an object of the appropriate class, or directly by
#'   calling \code{summary.rpartmv} regardless of the class of the object.}
#'   \item{predict.rpartmv}{This function is a method for the generic function
#'   predict for class \code{rpartmv}. It can be invoked by calling \code{predict}
#'   for an object of the appropriate class, or directly by calling
#'   \code{predict.rpartmv} regardless of the class of the object. }
#'   \item{rsq.rpartmv}{Two plots are produced. The labels are only appropriate
#'   for the \code{"anova"} method.}
#' }
#' 
#' @format A \code{rpartmv-class} object contains:
#' \describe{
#'   \item{frame}{Yet undocumented.}
#'   \item{where}{Yet undocumented.}
#'   \item{call}{The function call.}
#'   \item{terms}{The model's term.}
#'   \item{cptable}{A table of the error with respect to the number of
#'   splits.}
#'   \item{splits}{A table of the splits.}
#'   \item{method}{The method used, which is specified by argument \code{method}
#'   in function \code{\link{rpartmv}}.}
#'   \item{dissim}{Type of (implicit) dissimilarities.}
#'   \item{parms}{Yet undocumented.}
#'   \item{control}{List of control parameters used by \code{\link{rpartmv}},
#'   which is generated by function \code{\link{rpartmv.control}}.}
#'   \item{functions}{A list of functions used by other methods or functions
#'   taking an \code{rpartmv-class} object as an argument.}
#'   \item{y}{Yet undocumented.}
#'   \item{ordered}{Yet undocumented.}
#' }
#' 
#' @references
#' McCullagh P. and Nelder, J. A. (1989) Generalized Linear Models. London:
#' Chapman and Hall.
#' 
#' @importFrom grDevices dev.cur
#' @importFrom graphics legend lines par points segments text polygon
#' @importFrom stats delete.response model.frame naprint naresid
#' 
#' @seealso \code{\link{rpartmv}} for examples.
#' 
NULL
#' 
#' @describeIn rpartmv-class
#' 
#' Print an Rpartmv Object
#' 
#' This function prints an \code{rpartmv} object. It is a method for the generic
#' function \code{print} of class \code{rpartmv}.
#' 
#' @method print rpartmv
#' 
#' @export
print.rpartmv <- function(x, minlength = 0, spaces = 2, cp,
                          digits = getOption("digits"), ...) {
  if(!inherits(x, "rpartmv"))
    stop("Not legitimate rpartmv object")
  if (!is.null(x$frame$splits))
    x <- rpconvert(x)  #help for old objects
  if(!missing(cp))
    x <- prune.rpartmv(x, cp=cp)
  frame <- x$frame
  ylevel <- attr(x, "ylevels")
  node <- as.numeric(row.names(frame))
  depth <- tree.depth(node)
  indent <- paste(rep(" ", spaces*32), collapse="")
  #32 is the maximal depth
  if(length(node) > 1L) {
    indent <- substring(indent, 1, spaces * seq(depth))
    indent <- paste(c("", indent[depth]), format(node), ")", sep = "")
  } else indent <- paste(format(node), ")", sep = "")
  tfun <- (x$functions)$print
  if(!is.null(tfun)) {
    if(is.null(frame$yval2)) {
      yval <- tfun(frame$yval, ylevel, digits)
    } else yval <- tfun(frame$yval2,  ylevel, digits)
  } else yval <- format(signif(frame$yval, digits = digits))
  term <- rep(" ", length(depth))
  term[frame$var == "<leaf>"] <- "*"
  z <- labels(x, digits=digits, minlength=minlength, ...)
  n <- frame$n
  z <- paste(indent, z, n, format(signif(frame$dev, digits = digits)), yval,
             term)
  omit <- x$na.action
  if(length(omit)) {
    cat("n =", n[1L], " (", naprint(omit), ")\n\n", sep="")
  } else cat("n=", n[1L], "\n\n")
  #This is stolen, unabashedly, from print.tree
  if(x$method == "class") {
    cat("node), split, n, loss, yval, (yprob)\n")
  } else cat("node), split, n, deviance, yval\n")
  cat("      * denotes terminal node\n\n")
  cat(z, sep = "\n")
  return(invisible(x))
  #end of the theft
}
#' 
#'@describeIn rpartmv-class
#' 
#' Create Split Labels For an Rpartmv Object
#' 
#' This function provides labels for the branches of an \code{rpartmv} tree.
#' 
#' @method labels rpartmv
#' 
#' @export
labels.rpartmv <- function(object, digits = 4, minlength = 1, pretty,
                           collapse = TRUE, ...) {
  if(missing(minlength) && !missing(pretty)) {
    if(is.null(pretty)) minlength <-1
    else if(is.logical(pretty)) {
      if (pretty) minlength <- 4
      else        minlength <- 0
    }
    else minlength <- 0
  }
  ff <- object$frame
  n  <- nrow(ff)
  if(n == 1) return("root")  #special case of no splits
  is.leaf <- (ff$var == "<leaf>")
  whichrow <- !is.leaf
  vnames <- ff$var[whichrow]  #the variable names for the primary splits
  index <- cumsum(c(1, ff$ncompete + ff$nsurrogate + 1*(!is.leaf)))
  irow  <- index[c(whichrow, FALSE)]     #we only care about the primary split
  ncat  <- object$splits[irow, 2]
  # Now to work: first create labels for the left and right splits,
  #  but not for leaves of course
  #
  lsplit <- rsplit <- vector(mode='character', length= length(irow))
  if(any(ncat <2)) {  # any continuous vars ?
    jrow <- irow[ncat <2]
    cutpoint <- formatg(object$splits[jrow,4], digits)
    temp1 <- (ifelse(ncat<0, "< ", ">="))[ncat <2]
    temp2 <- (ifelse(ncat<0, ">=", "< "))[ncat <2]
    lsplit[ncat<2] <- paste(temp1, cutpoint, sep='')
    rsplit[ncat<2] <- paste(temp2, cutpoint, sep='')
  }
  if(any(ncat > 1)) { # any categorical variables ?
    xlevels <- attr(object, 'xlevels')
    #
    # jrow will be the row numbers of factors within lsplit and rsplit
    # crow the row number in "csplit"
    # and cindex the index on the "xlevels" list
    #
    jrow <- (seq(along=ncat))[ncat>1]
    crow <- object$splits[irow[ncat>1],4]    #row number in csplit
    cindex <- (match(vnames, names(xlevels)))[ncat >1]
    # Now, abbreviate the levels
    if(minlength ==1) {
      if(any(ncat>52))
        warning(paste("More than 52 levels in a predicting factor,",
                      "truncated for printout"))
      xlevels <- lapply(
        xlevels,
        function(z) {
          k <- length(z)
          k <- pmin(1L:k, 52L)
          c(letters, LETTERS)[k]
        }
      )
    } else if(minlength >1)
      xlevels <- lapply(xlevels, abbreviate, minlength, ...)
    # Now tuck in the labels
    # I'll let some other clever person vectorize this
    for(i in 1:(length(jrow))) {
      j <- jrow[i]
      splits <- object$csplit[crow[i],]
      # splits will contain 1=left, 2=right, 3= neither
      ltemp <- (1:length(splits))[splits== 1]
      rtemp <- (1:length(splits))[splits== 3]
      if(minlength==1) {
        lsplit[j] <- paste((xlevels[[cindex[i]]])[ltemp], collapse='')
        rsplit[j] <- paste((xlevels[[cindex[i]]])[rtemp], collapse='')
      } else {
        lsplit[j] <-paste((xlevels[[cindex[i]]])[ltemp], collapse=',')
        rsplit[j] <-paste((xlevels[[cindex[i]]])[rtemp], collapse=',')
      }
    }
  }
  if(!collapse) {  #called by no routines that I know of
    ltemp <- rtemp <- rep("<leaf>", n)
    ltemp[whichrow] <- lsplit
    rtemp[whichrow] <- rsplit
    return(cbind(ltemp, rtemp))
  }
  lsplit <- paste(ifelse(ncat<2, "", "="), lsplit, sep='')
  rsplit <- paste(ifelse(ncat<2, "", "="), rsplit, sep='')
  #
  # Now match them up to node numbers
  #   The output will have one label per row of object$frame, each
  #   corresponding the the line segement joining this node to its parent
  #
  varname <- (as.character(vnames))
  node <- as.numeric(row.names(ff))
  parent <- match(node %/% 2, node[whichrow])
  odd <- (as.logical(node %%2))
  labels <- vector('character', length=n)
  labels[odd] <- paste(varname[parent[odd]], rsplit[parent[odd]], sep="")
  labels[!odd]<- paste(varname[parent[!odd]],lsplit[parent[!odd]], sep="")
  labels[1] <- "root"
  labels
}
#' 
#' @describeIn rpartmv-class
#' 
#' Cost-complexity Pruning of an Rpartmv Object
#' 
#' Determines a nested sequence of subtrees of the supplied \code{rpartmv}
#' object by recursively \code{snipping} off the least important splits, based
#' on the complexity parameter (\code{cp}).
#' 
#' @method prune rpartmv
#' 
#' @export
prune.rpartmv <- function(object, cp, ...) {
  ff <- object$frame
  id <- as.integer(row.names(ff))
  toss <- id[ff$complexity <= cp &  ff$var!='<leaf>']#not a leaf
  if(length(toss) == 0L)
    return(object)   #all the tree is retained
  newx <- snip.rpartmv(object, toss)
  ## Now cut down the CP table
  temp <- pmax(object$cptable[,1], cp)
  keep <- match(unique(temp), temp)
  newx$cptable <- object$cptable[keep,,drop=FALSE]
  newx$cptable[max(keep),1] <- cp
  newx
}
#' 
#' @describeIn rpartmv-class
#' 
#' Plot an Rpartmv Object
#' 
#' Plots an rpartmv object on the current graphics device.
#' 
#' @method plot rpartmv
#' 
#' @export
plot.rpartmv <- function (x, uniform = FALSE, branch = 1, compress = FALSE,
                          nspace, margin = 0.0, minbranch = 0.3, bar = 0.03,
                          ...) {
  if (!inherits(x, "rpartmv")) 
    stop("Not an rpartmv object")
  if (!is.null(x$frame$splits)) 
    x <- rpconvert(x)
  if (compress & missing(nspace)) 
    nspace <- branch
  if (!compress) 
    nspace <- -1
  dev <- dev.cur()
  
  ##  added 15/04/13 
  parms <- list(uniform = uniform, branch = branch, nspace = nspace,
                minbranch = minbranch)
  temp <- rpartmvco(x,parms)
  xx <- temp$x
  yy <- temp$y
  temp1 <- range(xx) + diff(range(xx)) * c(-margin, margin)
  temp2 <- range(yy) + diff(range(yy)) * c(-margin, margin)
  plot(temp1, temp2, type = "n", axes = FALSE, xlab = "", ylab = "", 
       ...)
  node <- as.numeric(row.names(x$frame))
  temp <- rpartmv.branch(xx, yy, node, branch)
  if (branch > 0) 
    lines(c(xx[1], xx[1]), c(yy[1], yy[1] + bar*diff(range(yy))), ...)
  lines(c(temp$x), c(temp$y))
  invisible(list(x = xx, y = yy))
}
#' 
#' @describeIn rpartmv-class
#' 
#' Summarize an Rpartmv Object
#' 
#' Returns a detailed listing of an \code{rpartmv} object.
#' 
#' @method summary rpartmv
#' 
#' @export
summary.rpartmv <- function(object, cp = 0, digits = getOption("digits"), file,
                            ...) {
  
  if(!inherits(object, "rpartmv"))
    stop("Not legitimate rpartmv object")
  
  # If this is an older-style rpartmv object, convert it
  #  either way, rename it to "x" to save typing
  if(!is.null(object$frame$splits)) {
    x <- rpconvert(object)
  } else x <- object
  
  if(!missing(file)) {
    sink(file)
    on.exit(sink())
  }
  
  if(!is.null(x$call)) {
    cat("Call:\n")
    dput(x$call)
  }
  
  omit <- x$na.action
  n <- x$frame$n
  if(length(omit)) {
    cat("  n=", n[1L], " (", naprint(omit), ")\n\n", sep="")
  } else cat("  n =", n[1L], "\n\n")
  
  print(x$cptable, digits=digits)
  ff <- x$frame
  ylevel <- attr(x,'ylevels')
  id <- as.integer(row.names(ff))
  parent.id <- ifelse(id==1,1, floor(id/2))
  parent.cp <- ff$complexity[match(parent.id, id)]
  rows <- (1L:length(id))[parent.cp > cp]
  if(length(rows) > 0L) {
    rows <- rows[order(id[rows])]
  } else rows <- 1
  is.leaf <- (ff$var=='<leaf>')
  index <- cumsum(c(1, ff$ncompete + ff$nsurrogate + 1*(!is.leaf)))
  
  if(!all(is.leaf)) {  #skip these lines for a "no splits" tree
    sname <- dimnames(x$splits)[[1L]]
    cuts <- vector(mode='character', length=nrow(x$splits))
    temp <- x$splits[,2L]
    for(i in 1L:length(cuts)) {
      if(temp[i] == -1) {
        cuts[i] <-paste("<", format(signif(x$splits[i,4L], digits=digits)))
      } else if(temp[i] == 1) {
        cuts[i] <-paste("<", format(signif(x$splits[i,4L], digits=digits)))
      } else {
        paste(
          "splits as ",
          paste(c("L", "-", "R")[x$csplit[x$splits[i,4L], 1L:temp[i]]],
                collapse='', sep=''), collapse=''
        ) -> cuts[i]
      }
    }
    # S-PLUS 4.0 can't handle null vectors here
    if(any(temp < 2))
      cuts[temp < 2] <- format(cuts[temp < 2],justify="left")
    paste(
      cuts,
      ifelse(
        temp >= 2, ",",
        ifelse(temp == 1, " to the right,", " to the left, ")),
      sep = ''
    ) -> cuts
  }
  
  if(is.null(ff$yval2)) {
    x$functions$summary(
      ff$yval[rows], ff$dev[rows],
      ff$wt[rows], ylevel, digits) -> tprint
  } else {
    x$functions$summary(
      ff$yval2[rows,], ff$dev[rows],
      ff$wt[rows], ylevel, digits) -> tprint
  }
  
  for(ii in 1L:length(rows)) {
    i <- rows[ii]
    nn <- ff$n[i]
    twt <- ff$wt[i]
    cat("\nNode number ", id[i], ": ", nn, " observations", sep='')
    if(ff$complexity[i] < cp || is.leaf[i]) {
      cat("\n")
    } else {
      cat(",    complexity param=", format(signif(ff$complexity[i], digits)),
          "\n", sep="")
    }
    cat(tprint[ii], "\n")
    if(ff$complexity[i] > cp && !is.leaf[i]) {
      sons <- 2*id[i] + c(0,1)
      sons.n <- ff$n[match(sons, id)]
      cat("  left son=", sons[1], " (", sons.n[1], " obs)",
          " right son=", sons[2], " (", sons.n[2], " obs)", sep='')
      j <- nn - (sons.n[1] + sons.n[2])
      if(j>1) {
        cat(", ", j, " observations remain\n", sep='')
      } else if(j==1) {
        cat(", 1 observation remains\n")
      } else cat("\n")
      cat("  Primary splits:\n")
      j <- seq(index[i], length=1+ff$ncompete[i])
      if(all(nchar(cuts[j]) < 25L)) {
        temp <- format(cuts[j], justify="left")
      } else  temp <- cuts[j]
      cat(paste("      ", format(sname[j], justify="left"), " ", temp,
                " improve=", format(signif(x$splits[j,3], digits)),
                ", (", nn - x$splits[j,1], " missing)", sep=''),
          sep="\n")
      if(ff$nsurrogate[i] > 0) {
        cat("  Surrogate splits:\n")
        j <- seq(1 +index[i] + ff$ncompete[i], length=ff$nsurrogate[i])
        agree <- x$splits[j,3]
        if(all(nchar(cuts[j]) < 25L)) {
          temp <- format(cuts[j], justify="left")
        } else  temp <- cuts[j]
        if(ncol(x$splits) == 5L) {
          adj   <- x$splits[j,5]
          cat(paste("      ", format(sname[j], justify="left"), " ", temp,
                    " agree=", format(round(agree, 3L)),
                    ", adj=" , format(round(adj, 3L)),
                    ", (", x$splits[j,1L], " split)", sep=''),
              sep="\n")
        } else {
          #an older style rpart object -- no adj value present
          cat(paste("      ", format(sname[j], justify="left"), " ", temp,
                    " agree=", format(round(agree, 3L)),
                    ", (", x$splits[j,1L], " split)", sep=''),
              sep="\n")
        }
      }
    }
  }
  
  cat("\n")
  invisible(x)
}
#' 
#' @describeIn rpartmv-class
#' 
#' Predictions from a Fitted Rpartmv Object
#' 
#' Returns a vector of predicted responses from a fitted \code{rpartmv} object.
#' 
#' @method predict rpartmv
#' 
#' @export
predict.rpartmv <- function (object, newdata = list(),
                             type = c("vector","prob","class","matrix","where"),
                             ...) {
  if(!inherits(object, "rpartmv"))
    stop("Not legitimate tree")
  mtype <- missing(type)
  type <- match.arg(type)
  if(missing(newdata)) {
    where <- object$where
  } else {
    if(is.null(attr(newdata, "terms"))) {
      Terms <- delete.response(object$terms)
      act <- (object$call)$na.action
      if(is.null(act)) 
        act <- na.rpartmv
      model.frame(
        Terms,
        newdata,
        na.action = act,
        xlev = attr(object, "xlevels")
      ) -> newdata
    }
    where <- pred.rpartmv(object, rpartmv.matrix(newdata))
  }
  if(type=="where") {
    return(where)
  } else {
    frame <- object$frame
    method <- object$method
    ylevels <- attr(object, "ylevels")
    nclass <- length(ylevels)
    if(mtype && nclass > 0) 
      type <- "prob"
    if(mtype && method == "mrt") 
      type <- "matrix"
    if(type == "vector" || (type == "matrix" && is.null(frame$yval2))) {
      pred <- frame$yval[where]
      names(pred) <- names(where)
    } else if(type == "matrix") {
      pred <- frame$yval2[where,]
      dimnames(pred) <- list(names(where), NULL)
    } else if(type == "class" && nclass > 0) {
      pred <- factor(ylevels[frame$yval[where]], levels = ylevels)
      names(pred) <- names(where)
    } else if(type == "prob" && nclass > 0) {
      pred <- frame$yval2[where, 1 + nclass + 1:nclass]
      dimnames(pred) <- list(names(where), ylevels)
    } else stop("Invalid prediction for rpartmv object")
    if(missing(newdata) && !is.null(object$na.action)) 
      pred <- naresid(object$na.action, pred)
    pred
  }
}
#' 
#' @describeIn rpartmv-class
#' 
#' Residuals From a Fitted Rpartmv Object
#' 
#' Method for \code{residuals} for an \code{rpartmv} object.
#' 
#' @method residuals rpartmv
#' 
#' @export
residuals.rpartmv <- function(object, type = c("usual","pearson","deviance"),
                              ...) {
  if(!inherits(object, "rpartmv"))
    stop("Not legitimate rpartmv object")
  y <- object$y
  if(is.null(y))
    y <- model.extract(model.frame(object), "response")
  frame <- object$frame
  type <- match.arg(type)
  if(is.na(match(type, c("usual", "pearson", "deviance"))))
    stop("Don't know about this type of residual")
  if(object$method=='class') {
    ylevels <- attr(object, "ylevels")
    nclass <- length(ylevels)
    if(type == "usual") {
      yhat <- frame$yval[object$where]
      loss <- object$parms$loss
    } else {
      yprob <- frame$yval2[object$where, 1 + nclass + 1:nclass]
      yhat <- yprob[cbind(seq(y), unclass(y))]
    }
    switch(
      type,
      usual = loss[cbind(y, yhat)],
      pearson = (1 - yhat)/yhat,
      deviance = -2 * log(yhat)
    ) -> resid
  } else if(object$method=='poisson' || object$method=='exp') {
    lambda <- (object$frame$yval)[object$where]
    time   <- y[,1L]  # observation time in new data
    events <- y[,2L]  # number of events, in new data
    expect <- lambda * time #expected number of events
    temp <- ifelse(expect==0, .0001, 0)  #failsafe for log(0)
    switch(
      type,
      usual = events - expect,
      pearson = (events - expect)/sqrt(temp),
      deviance = sign(events- expect)*
        sqrt(2*(events*log(events/temp) - (events-expect)))
    ) -> resid
  } else resid <- y - frame$yval[object$where]
  
  names(resid) <- names(y)
  #Expand out the missing values in the result
  if(!is.null(object$na.action))
    resid <- naresid(object$na.action, resid)
  resid
}
#' 
#' @describeIn rpartmv-class
#' 
#' Plots the Approximate R-Square for the Different Splits
#' 
#' Produces 2 plots.  The first plots the r-square (apparent and apparent - from
#' cross-validation) versus the number of splits. The second plots the Relative
#' Error(cross-validation) +/- 1-SE from cross-validation versus the number of
#' splits.
#' 
#' @method rsq rpartmv
#' 
#' @export
rsq.rpartmv <- function(object, ...) {
  p.rpartmv <- printcp(object)
  xstd <- p.rpartmv[,5L]
  xerror <- p.rpartmv[,4L]
  rel.error <- p.rpartmv[,3L]
  nsplit <- p.rpartmv[,2L]
  method <- object$method
  if(!method == "anova")
    cat("May not be applicable for this method\n")
  plot(nsplit, 1 - rel.error, xlab="Number of Splits", ylab="R-square",
       ylim=c(0,1), type="o")
  par(new=TRUE)
  plot(nsplit, 1 - xerror, type="o", ylim=c(0,1), lty=2L, xlab="", ylab="", ...)
  legend(0, 1, c("Apparent","X Relative"), lty=1L:2L)
  ylim <- c(min(xerror-xstd) -.1, max(xerror + xstd) + .1)
  plot(nsplit, xerror, xlab="Number of Splits", ylab="X Relative Error",
        ylim=ylim, type="o")
  segments(nsplit, xerror - xstd, nsplit, xerror + xstd)
  invisible(NULL)
}
#' 
#' @describeIn rpartmv-class
#' 
#' Place Text on a Dendrogram
#' 
#' Labels the current plot of the tree dendrogram with text.
#' 
#' @method text rpartmv
#' 
#' @export
text.rpartmv <- function(x, splits = TRUE, which = 4, label = "yval",
                         FUN = text, all.leaves = FALSE, pretty = NULL,
                         digits = getOption("digits") - 2, tadj = 0.65,
                         stats = TRUE, use.n = FALSE, bars = TRUE,
                         legend = FALSE, xadj = 1, yadj = 1, bord = FALSE,
                         big.pts = FALSE, uniform = FALSE, branch = 1,
                         nspace = -1, minbranch = 0.3, ...) {
  if(!inherits(x, "rpartmv")) 
    stop("Not legitimate rpartmv")
  if(!is.null(x$frame$splits)) 
    x <- rpconvert(x)
  frame <- x$frame
  col <- names(frame)
  method <- x$method
  ylevels <- attr(x, "ylevels")
  if(!is.null(ylevels <- attr(x, "ylevels"))) 
    col <- c(col, ylevels)
  if(is.na(match(label, col))) 
    stop("Label must be a column label of the frame component of the tree")
  cxy <- par("cxy")
  if(!is.null(srt <- list(...)$srt) && srt == 90) 
    cxy <- rev(cxy)
  parms <- list(uniform = uniform, branch = branch, nspace = nspace,
                minbranch = minbranch)
  xy <- rpartmvco(x,parms)
  node <- as.numeric(row.names(x$frame))
  is.left <- (node%%2 == 0)
  node.left <- node[is.left]
  parent <- match(node.left/2, node)
  bars <- bars & is.matrix(frame$yval2)
  text.adj <- ifelse(bars, yadj * diff(range(xy$y))/12, 0)
  if(splits) {
    left.child <- match(2 * node, node)
    right.child <- match(node * 2 + 1, node)
    rows <- labels(x, pretty = pretty)
    if(which == 1) {
      FUN(xy$x, xy$y + tadj * cxy[2], rows[left.child], ...)
    } else {
      if(which == 2 | which == 4) 
        FUN(xy$x, xy$y + tadj * cxy[2], rows[left.child], pos = 2, ...)
      if(which == 3 | which == 4) 
        FUN(xy$x, xy$y + tadj * cxy[2], rows[right.child], pos = 4, ...)
    }
  }
  leaves <- if (all.leaves) rep(TRUE, nrow(frame)) else frame$var == "<leaf>"
  if(stats) {
    if(is.null(frame$yval2)) {
      x$functions$text(
        yval = frame$yval[leaves], 
        dev = frame$dev[leaves], wt = frame$wt[leaves], 
        ylevel = ylevels, digits = digits, n = frame$n[leaves], 
        use.n = use.n) -> stat
    } else {
      x$functions$text(
        yval = frame$yval2[leaves,], dev = frame$dev[leaves],
        wt = frame$wt[leaves], ylevel = ylevels, digits = digits,
        n = frame$n[leaves], use.n = use.n) -> stat
    }
    FUN(xy$x[leaves], xy$y[leaves] - tadj * cxy[2] - text.adj, stat, adj = 0.5,
        ...)
  }
  if(bars) {
    bar.vals <- x$functions$bar(yval2 = frame$yval2)
    sub.barplot(xy$x, xy$y, bar.vals, leaves, xadj = xadj, 
                yadj = yadj, bord = bord, line = TRUE,
                col = c("lightblue","blue","darkblue"))
    rx <- range(xy$x)
    ry <- range(xy$y)
    if(!is.null(ylevels)) {
      bar.labs <- ylevels
    } else bar.labs <- dimnames(x$y)[[2L]]
    if(legend & !is.null(bar.labs))
      legend(min(xy$x) - 0.1*rx, max(xy$y) + 0.05*ry, bar.labs, 
             col=c("lightblue","blue","darkblue"), pch=15L, bty="n", ...)
  }
  if(big.pts)
    points(xy$x[leaves], xy$y[leaves], pch=16L, cex=3*par()$cex,
           col=2L:(sum(leaves) + 1L))
  invisible(NULL)
}
#' 
#' @describeIn rpartmv-class
#' 
#' Mean-Variance Plot for an Rpartmv Object
#' 
#' Creates a plot on the current graphics device of the deviance of the node
#' divided by the number of observations at the node.  Also returns the node
#' number.
#' 
#' @method meanvar rpartmv
#' 
#' @export
meanvar.rpartmv <- function(object, xlab = "ave(y)", ylab = "ave(deviance)",
                            ...) {
  if(!object$method=="anova")
    stop("Plot not useful for classification or poisson trees")
  frame <- object$frame
  frame <- frame[frame$var == "<leaf>",]
  x <- frame$yval
  y <- frame$dev/frame$n
  label <- row.names(frame)
  plot(x, y, xlab=xlab, ylab=ylab, type="n", ...)
  text(x, y, label)
  invisible(list(x=x, y=y, label=label))
}
