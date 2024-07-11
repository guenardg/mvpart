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
##    **Internal function for adapting rpartmv to different types of responses**
##
## **************************************************************************
##

#' @importFrom stats approx quantile

rpartmv.anova <- function(y, offset, parms, wt) {
  if(!is.null(offset)) y <- y-offset
  list(
    y = y, parms = 0, numresp = 1, numy = 1,
    summary = function(yval, dev, wt, ylevel, digits) {
      paste("  mean=", formatg(yval, digits),
            ", MSE=" , formatg(dev/wt, digits),
            sep='')
    },
    text = function(yval, dev, wt, ylevel, digits, n, use.n) {
      if(use.n) {
        paste(formatg(yval,digits), "\nn=", n, sep="")
      } else {
        paste(formatg(yval,digits))}
      }
  )
}

rpartmv.mrt <- function (y, offset, parms, wt) {
  if(!is.null(offset))
    y <- y - offset
  ny <- ifelse(is.null(dim(y)), 1, dim(y)[2L])
  list(
    y = y, parms = 0, numresp = ny, numy = ny,
    summary = function(yval, dev, wt, ylevel, digits) {
      paste("  Means=",
            apply(formatg(yval, digits - 3), 1L,
                  paste, collapse = ",", sep = ""), ", Summed MSE=",
            formatg(dev/wt, digits), sep = "")
    },
    text = function(yval, dev, wt, ylevel, digits, n, use.n) {
      if(use.n) {
        paste(formatg(dev, digits), " : n=", n, sep = "")
      } else {
        paste(formatg(dev, digits))
      }
    },
    bar = function(yval2) yval2
  )
}

rpartmv.poisson <- function(y, offset, parms, wt) {
  if(is.matrix(y)) {
    if(ncol(y) != 2) stop("response must be a 2 column matrix or a vector")
    if(!is.null(offset)) y[,1L] <- y[,1L]*exp(offset)
  } else {
    if(is.null(offset)) y <- cbind(1,y)
    else y <- cbind(exp(offset), y)
  }
  if(any(y[,1L] <= 0)) stop("Observation time must be >0")
  if(any(y[,2L] < 0))  stop("Number of events must be >=0")
  if(missing(parms)) {
    parms <- c(shrink=1, method=1)
  } else {
    parms <- as.list(parms)
    if(is.null(names(parms))) stop("You must input a named list for parms")
    parmsNames <- c("method", "shrink")
    indx <- pmatch(names(parms), parmsNames, nomatch = 0)
    if(any(indx == 0L)) {
      stop(paste("parms component not matched: ", names(parms)[indx == 0L]))
    } else names(parms) <- parmsNames[indx]
    if(is.null(parms$method)) {
      method <- 1L
    } else method <- pmatch(parms$method, c("deviance", "sqrt"))
    if(is.null(method)) stop("Invalid error method for Poisson")
    if (is.null(parms$shrink)) {
      shrink <- 2 - method
    } else shrink <- parms$shrink
    if(!is.numeric(shrink) || shrink < 0L)
      stop("Invalid shrinkage value")
    parms <- c(shrink=shrink, method=method)
  }
  list(
    y = y, parms = parms, numresp = 2, numy = 2,
    summary = function(yval, dev, wt, ylevel, digits) {
      paste(
        "  events=", formatg(yval[,2L]), ",  estimated rate=" ,
        formatg(yval[,1L], digits), " , mean deviance=",
        formatg(dev/wt, digits), sep = ""
      )
    },
    text = function(yval, dev, wt, ylevel, digits, n, use.n) {
      if(use.n) {
        paste(
          formatg(yval[,1L],digits), "\n", formatg(yval[,2L]), "/", n, sep=""
        )
      } else paste(formatg(yval[,1L],digits))
    }
  )
}

rpartmv.class <- function (y, offset, parms, wt) {
  if(!is.null(offset)) 
    stop("No offset allowed in classification models")
  fy <- as.factor(y)
  y <- as.integer(fy)
  numclass <- max(y[!is.na(y)])
  counts <- tapply(wt, factor(y, levels = 1:numclass), sum)
  counts <- ifelse(is.na(counts), 0, counts)
  numresp <- 1 + numclass
  if(missing(parms) || is.null(parms)) {
    list(
      prior = counts/sum(counts),
      loss = matrix(rep(1, numclass^2) - diag(numclass), numclass),
      split = 1
    ) -> parms
  } else if (is.list(parms)) {
    if(is.null(names(parms))) 
      stop("The parms list must have names")
    temp <- pmatch(names(parms), c("prior", "loss", "split"), nomatch = 0)
    if(any(temp == 0))
      stop(paste("parms component not matched:", (names(parms))[temp == 0]))
    names(parms) <- c("prior", "loss", "split")[temp]
    if(is.null(parms$prior)) {
      temp <- c(counts/sum(counts))
    } else {
      temp <- parms$prior
      if(sum(temp) != 1) 
        stop("Priors must sum to 1")
      if(any(temp < 0)) 
        stop("Priors must be >= 0")
      if(length(temp) != numclass) 
        stop("Wrong length for priors")
    }
    if(is.null(parms$loss)) {
      temp2 <- 1 - diag(numclass)
    } else {
      temp2 <- parms$loss
      if(length(temp2) != numclass^2) 
        stop("Wrong length for loss matrix")
      temp2 <- matrix(temp2, ncol = numclass)
      if(any(diag(temp2) != 0)) 
        stop("Loss matrix must have zero on diagonals")
      if(any(temp2 < 0)) 
        stop("Loss matrix cannot have negative elements")
      if(any(rowSums(temp2) == 0)) 
        stop("Loss matrix has a row of zeros")
    }
    if(is.null(parms$split)) {
      temp3 <- 1
    } else {
      temp3 <- pmatch(parms$split, c("gini", "information"))
      if(is.null(temp3)) 
        stop("Invalid splitting rule")
    }
    parms <- list(prior = temp, loss = matrix(temp2, numclass), split = temp3)
  } else stop("Parameter argument must be a list")
  list(
    y = y, parms = parms, numresp = numclass + 1, counts = counts,
    ylevels = levels(fy), numy = 1,
    print = function(yval, ylevel, digits) {
      if(is.null(ylevel)) {
        temp <- as.character(yval[,1L])
      } else temp <- ylevel[yval[,1L]]
      nclass <- (ncol(yval) - 1)/2
      if(nclass < 5L) {
        yprob <- format(yval[, 1 + nclass + 1:nclass],
                        digits = digits, nsmall = digits)
      } else yprob <- formatg(yval[, 1 + nclass + 1:nclass], digits = 2)
      if (is.null(dim(yprob)))
        yprob <- matrix(yprob, ncol = length(yprob))
      temp <- paste(temp, " (", yprob[, 1], sep = "")
      for(i in 2:ncol(yprob))
        temp <- paste(temp, yprob[,i], sep = " ")
      temp <- paste(temp, ")", sep = "")
      temp
    },
    summary = function(yval, dev, wt, ylevel, digits) {
      nclass <- (ncol(yval) - 1)/2
      group <- yval[,1L]
      counts <- yval[,1L + (1L:nclass)]
      yprob <- yval[, 1L + nclass + 1L:nclass]
      if(!is.null(ylevel)) group <- ylevel[group]
      temp1 <- formatg(counts, format = "%5g")
      temp2 <- formatg(yprob, format = "%5.3f")
      if(nclass > 1) {
        temp1 <- apply(matrix(temp1, ncol = nclass), 1L, paste, collapse = " ")
        temp2 <- apply(matrix(temp2, ncol = nclass), 1L, paste, collapse = " ")
      }
      paste(
        "  predicted class=", format(group, justify = "left"),
        "  expected loss=", formatg(dev/wt, digits), "\n", "    class counts: ",
        temp1, "\n", "   probabilities: ", temp2, sep = ""
      )
    },
    text = function(yval, dev, wt, ylevel, digits, n, use.n) {
      nclass <- (ncol(yval) - 1)/2
      group <- yval[, 1]
      counts <- yval[, 1 + (1:nclass)]
      if(!is.null(ylevel)) group <- ylevel[group]
      temp1 <- formatg(counts, digits)
      if(nclass > 1) {
        temp1 <- apply(matrix(temp1, ncol = nclass), 1L, paste, collapse = "/")
      }
      if(use.n) {
        out <- paste(format(group, justify = "left"), "\n", temp1, sep = "")
      } else {
        out <- format(group, justify = "left")
      }
      out
    },
    bar = function(yval2) yval2[,2L:((ncol(yval2) + 1L)/2)]
  )
}

rpartmv.dist <- function (y, offset, parms, wt) {
  if(!is.null(offset)) 
    y <- y - offset
  ny <-  1
  list(
    y = y, parms = 0, numresp = ny, numy = ny,
    summary = function(yval, dev, wt, ylevel, digits) {
      paste(
        "  Means=",
        apply(formatg(yval, digits - 3), 1L, paste, collapse = ",", sep = ""),
        ", Summed MSE=", formatg(dev/wt, digits), sep = ""
      )
    },
    text = function(yval, dev, wt, ylevel, digits, n, use.n) {
      if(use.n) {
        paste(formatg(dev, digits), " : n=", n, sep = "")
      } else {
        paste(formatg(dev, digits))
      }
    },
    bar = function(yval2) yval2
  )
}

rpartmv.exp <- function(y, offset, parms, wt) {
  if (!inherits(y, "Surv"))
    stop("Response must be a survival object - use the Surv() function")
  ny <- ncol(y)
  n  <- nrow(y)
  status <- y[,ny]
  if(any(y[,1L] <= 0)) stop("Observation time must be >0")
  if(all(status == 0)) stop("No deaths in data set")
  time <- y[ ,ny - 1L]
  dtimes <- sort(unique(time[status == 1]))        # unique death times
  .C("rpartexp2",
     as.integer(length(dtimes)),
     as.double(dtimes),
     as.double(.Machine$double.eps),
     keep = integer(length(dtimes)),
     PACKAGE = "mvpart")$keep -> temp
  dtimes <- dtimes[temp == 1]
  if (length(dtimes) > 1000) dtimes <- quantile(dtimes, 0:1000/1000)
  # The last interval goes to the max time in the data set
  itable <- c(0, dtimes[-length(dtimes)], max(time)) # set of intervals
  drate2 <- function(n, ny, y, wt, itable) {
    # An alternative to the drate1 function
    # Why?  The pyears2 routine changed in 6/2001, with the inclusion
    #  of case weights.  We need the newer version.  If you have the
    #  older version of the survival library, the above will crash S.
    # How to tell -- list the pyears function, and see whether it's
    #  call to pyears2 has weights in the argument list.
    #
    time <- y[, ny - 1L]
    status <- y[,ny]
    ilength <- diff(itable)                   #lengths of intervals
    ngrp <- length(ilength)                   #number of intervals
    # The code below is as opaque as any I've written, all in the
    #  service of "no for loops".
    # First, 'index' gives the time interval (as defined by itable)
    #  in which the end of each observation's follow-up (time) lies.
    #  Then 'itime' will be the amount of time spent in that last
    #  interval, which is of course somewhat less than ilength.
    index <- unclass(cut(time, itable, include.lowest=TRUE))
    itime <- time - itable[index]
    if(ny ==3) {
      # there is both a start time and a stop time
      #  compute the amount of time NOT spent in the interval that
      #  the start time lies in.
      # stime <- y[,1]   #start time for each interval
      index2<- unclass(cut(y[,1], itable, include.lowest=TRUE))
      itime2<- y[,1] - itable[index2]
    }
    # Compute the amount of person-years in each of the intervals
    #   This is:  (width of interval) * (number of "time" elements that
    #                                     end in an interval farther right)
    #            + (ending times in this interval)
    # By construction, I know that there is at least 1 obs in each of the
    #  intervals, so tab1 is of a determined length
    tab1 <- table(index)
      temp <- rev(cumsum(rev(tab1)))  #cumsum, counting from the right
    pyears <- ilength * c(temp[-1], 0) +     tapply(itime, index, sum)
    if(ny==3) {
      # subtract off the time before "start"
      tab2 <- table(index2, levels=1:ngrp) #force the length of tab2
      temp <- rev(cumsum(rev(tab2)))
      py2  <-  ilength * c(0, temp[-ngrp]) +  tapply(itime2, index2, sum)
      pyears <- pyears - py2
    }
    deaths <- tapply(status, index, sum)
    rate <- deaths/pyears   #hazard rate in each interval
    rate
  }
  #
  # Now, compute the "new y" for each observation.
  #  This is a stretching of the time axis
  # The cumulative hazard over each interval is rate*length(interval),
  #  and is the basis of the rescaling.
  rate <- drate2(n, ny, y, wt, itable)
  cumhaz <- cumsum(c(0, rate*diff(itable)))
  newy <- approx(itable, cumhaz, time)$y
  if(ny==3) newy <- newy - approx(itable, cumhaz, y[,1])$y
  if(length(offset)==n)  newy <- newy * exp(offset)
  if(missing(parms)) {
    parms <- c(shrink=1, method=1)
  } else {
    parms <- as.list(parms)
    if(is.null(names(parms)))
      stop("You must input a named list for parms")
    parmsNames <- c("method", "shrink")
    indx <- pmatch(names(parms), parmsNames, nomatch= 0)
    if(any(indx==0))
      stop(paste("parms component not matched: ",
                 names(parms)[indx==0]))
    else names(parms) <- parmsNames[indx]
    if(is.null(parms$method)) method <- 1
    else method <- pmatch(parms$method, c("deviance", "sqrt"))
    if(is.na(method)) stop("Invalid error method for Poisson")
    if(is.null(parms$shrink)) shrink <- 2-method
    else shrink <- parms$shrink
    if(!is.numeric(shrink) || shrink < 0)
      stop("Invalid shrinkage value")
    parms <- c(shrink=shrink, method=method)
  }
  list(
    y = cbind(newy, y[,2]), parms=parms, numresp=2, numy=2,
    summary = function(yval, dev, wt, ylevel, digits) {
      paste("  events=", formatg(yval[,2]), ",  estimated rate=" ,
            formatg(yval[,1], digits), " , mean deviance=",
            formatg(dev/wt, digits), sep = "")
    },
    text = function(yval, dev, wt, ylevel, digits, n, use.n) {
      if(use.n) {
        paste(
          formatg(yval[,1], digits), "\n", formatg(yval[,2]),"/",n, sep=""
        )
      } else paste(formatg(yval[,1],digits))
    }
  )
}
