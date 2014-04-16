## Jesse Engreitz

suppressPackageStartupMessages(library("optparse"))

option.list <- list(
  make_option(c("-i", "--input"), type="character", help="Input PSL file to filter"),
#  make_option(c("-l", "--layout"), type="character", help="Full probe layout file"),
  make_option(c("-s", "--stats"), type="character", help="Output file for statistics on the number of probes removed per gene"),
  make_option(c("-r", "--remove"), type="character", help="Output file for list of probes to remove"),
#  make_option(c("-o", "--output"), type="character", help="Output filtered probe layout file"),
  make_option(c("-p", "--output"), type="character", help="Output filtered PSL file"),
  make_option(c("--region"), type="character", default=NULL, help="Allow multiple matches within this region (e.g., for Firre)")
                    )
opt <- parse_args(OptionParser(option_list=option.list))

## Load PSL file
input <- read.delim(opt$input, header=F)
colnames(input)[c(1,10,14,16,17)] <- c("matches","oligo.name","chr","target.start","target.end")
input.filtered <- subset(input, matches >= 30 & target.end-target.start-matches < matches)
input.ref <- input.filtered
cutoff <- 1


## Parse input in format chrX:1,000,000-2,000,000
parseRegionName <- function(region) {
  region <- gsub(',','',region)
  if (length(region) == 1) {
    split1 <- strsplit(region, ":")[[1]]
    split2 <- strsplit(split1[2], "-")[[1]]
    result <- list(chr=split1[1], start=as.numeric(split2[1]), end=as.numeric(split2[2]))
    return(result)
  } else {
    result <- data.frame(t(sapply(as.character(region), parseRegionName)))
    result$chr <- unlist(result$chr)
    result$start <- as.numeric(unlist(result$start))
    result$end <- as.numeric(unlist(result$end))
    return(result)
  }
}


## Filter out matches within the region if needed
if (!is.null(opt$region)) {
  region <- parseRegionName(opt$region)
  input.filtered <- subset(input.filtered, (chr != region$chr) | (target.end > region$end) | (target.start < region$start))
  cutoff <- 0  ## any probes that match outside this region will be removed
}

#to.keep <- factor(names(which(table(input.filtered$oligo.name) <= 1)), levels=levels(input.filtered$oligo.name))
#cat(paste("BLAT Filter: Removing",length(levels(to.keep))-length(to.keep),"of",length(levels(to.keep)),"probes.\n\n"))

to.remove <- factor(names(which(table(input.filtered$oligo.name) > cutoff)), levels=levels(input.ref$oligo.name))
cat(paste("BLAT Filter: Removing",length(to.remove),"of",length(levels(to.remove)),"probes.\n\n"))
write.table(to.remove, file=opt$remove, sep='\t', quote=F, col.names=F, row.names=F)

to.keep <- levels(input.ref$oligo.name)[!(levels(input.ref$oligo.name) %in% to.remove)]
input.filtered <- subset(input.ref, !(oligo.name %in% to.remove))

genesFromOligos <- function(x) {
  unlist(lapply(strsplit(x,"_"),function(l) paste(l[1:(length(l)-3)],collapse="_")))
}
all.genes <- genesFromOligos(levels(to.remove))
kept.genes <- genesFromOligos(to.keep)
stats <- merge(data.frame(table(all.genes)), data.frame(table(kept.genes)), by.x=1, by.y=1, all.x=T)
colnames(stats) <- c("Transcript","Total.Oligos","Passed.Filter")
write.table(stats, file=opt$stats, sep='\t', quote=F, col.names=T, row.names=F)

write.table(input.filtered, file=opt$output, sep='\t', quote=F, col.names=F, row.names=F)

#layout.filtered <- subset(layout, !(Probe_ID %in% to.remove))
#write.table(layout.filtered, file=opt$output, sep='\t', quote=F, col.names=T, row.names=F)
