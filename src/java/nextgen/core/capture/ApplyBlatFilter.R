## Jesse Engreitz

suppressPackageStartupMessages(library("optparse"))

option.list <- list(
  make_option(c("-i", "--input"), type="character", help="Input oligo design file"),
  make_option(c("-r", "--remove"), type="character", help="Common-separated file list of probes to remove"),
  make_option(c("--ignore"), default=NULL, type="character", help="Comma-separated list of transcripts to ignore the BLAT filter"),
  make_option(c("-s", "--stats"), type="character", help="Output file for filtering stats"),
  make_option(c("-o", "--output"), type="character", help="Output filtered oligo design file")
                    )
opt <- parse_args(OptionParser(option_list=option.list))


genesFromOligos <- function(x) {
  unlist(lapply(strsplit(x,"_"),function(l) paste(l[1:(length(l)-3)],collapse="_")))
}

layout <- read.delim(opt$input)
to.remove <- c()
for (file in strsplit(opt$remove,',')[[1]]) {
  to.remove <- c(to.remove, read.delim(file, stringsAsFactors=F)[,1])
}
genes.to.remove <- factor(genesFromOligos(to.remove))
to.remove <- factor(to.remove, levels=levels(layout$Probe_ID))

if (!is.null(opt$ignore)) {
  ignore <- factor(strsplit(opt$ignore,',')[[1]], levels=levels(genes.to.remove))
  cat(paste("Accepting all probes from the following transcripts:\n"))
  print(ignore)
  
  to.remove <- to.remove[!(genes.to.remove %in% ignore)]
}

cat(paste("Removing",length(to.remove),"probes.\n"))

layout.filtered <- subset(layout, !(Probe_ID %in% to.remove))

genes.all <- table(layout$Parent_sequence)
genes.filtered <- table(layout.filtered$Parent_sequence)
stats <- merge(data.frame(genes.all), data.frame(genes.filtered), by.x=1, by.y=1, all.x=T)
colnames(stats) <- c("Parent_sequence","Total.Oligos","Passed.Filter")
write.table(stats, file=opt$stats, sep='\t', quote=F, col.names=T, row.names=F)
              
#layout.filtered$Probeset_size <- sapply(layout$Probe_ID, function(x) genes.filtered[as.character(as.matrix(x))])

write.table(layout.filtered, file=opt$output, sep='\t', quote=F, col.names=T, row.names=F)
