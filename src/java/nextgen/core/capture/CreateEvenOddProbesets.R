## Jesse Engreitz
## April 10, 2014

suppressPackageStartupMessages(library("optparse"))

option.list <- list(
  make_option(c("-i","--input"),type="character",help="Input design file; must have a header row with column names, including Start, End, and Probe_set"),
  make_option(c("-o","--output"),type="character",help="Output design file filtered to even/odd probesets for each existing Probe_set."))

opt <- parse_args(OptionParser(option_list=option.list))

getTilingPathFromFullDesign <- function(x, evenOdd=TRUE) {
  ## Start with the first probe and create a single tiling path by
  ## iteratively picking the next probe that does not overlap and has
  ## passed filters up to this point. 
  ##
  ## Find all greedy tiling paths
  ##
  ## Then, choose the tiling path with the most probes and lowest variance in distance
  ## between probes

  stopifnot(all(x$Probe_set== x$Probe_set[1]))
  x <- x[order(x$Start),]
  x$Probe_set <- as.character(as.matrix(x$Probe_set))
  
  tiling.paths <- list()
  
  first.end <- x$End[1]
  curr.start.index <- 1
 
  while (x$Start[curr.start.index] < first.end & curr.start.index <= nrow(x)) {
    ## Start a new tiling path
    indices.to.include <- c()
    curr.probe.end <- -1
    for (i in curr.start.index:nrow(x)) {
      if (x$Start[i] > curr.probe.end) {
        indices.to.include <- c(indices.to.include, i)
        curr.probe.end <- x$End[i]
      }
    }
    
    tiling.paths[[curr.start.index]] <- x[indices.to.include,]
    curr.start.index <- curr.start.index + 1
  }
  
  getDistancesBetween <- function(nums) {
    if (length(nums) == 1) return(0)
    return(nums[-1] - nums[-length(nums)])
  }
  
  n.probes <- unlist(lapply(tiling.paths, nrow))
  probe.dist.var <- unlist(lapply(tiling.paths, function(df) {
    var(getDistancesBetween(c(x$Start[1], x$End[nrow(x)], with(df, Start+(End-Start)/2))))
  }))

  which.path <- which(n.probes == max(n.probes))
  if (length(which.path) > 1) {
    which.path <- which.min(probe.dist.var)
  }

  result <- tiling.paths[[which.path]]

  if (!evenOdd | nrow(result) == 1) {
    result$Probe_set_size <- nrow(result)
  } else {
    ## Split into even and odd probesets
    which.odd <- seq(1,nrow(result),2)
    which.even <- seq(2,nrow(result),2)

    result$Probe_set[which.odd] <- paste(result$Probe_set[which.odd],"_odd",sep='')
    result$Probe_set[which.even] <- paste(result$Probe_set[which.even],"_even",sep='')

    result$Probe_set_size[which.odd] <- length(which.odd)
    result$Probe_set_size[which.even] <- length(which.even)
    result <- result[order(result$Probe_set),]
  }
  
  return(result)
}

                    
x <- read.delim(opt$input)
result <- tapply(1:nrow(x), x$Probe_set, function(i) getTilingPathFromFullDesign(x[i,]))
result <- do.call(rbind, result)
write.table(result, file=opt$output, sep='\t', quote=F, col.names=T, row.names=F)

