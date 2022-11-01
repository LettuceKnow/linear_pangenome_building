#!/usr/bin/env Rscript
require(ggplot2)
args = commandArgs(trailingOnly=TRUE) #<depth output> <breadth output> <plot output> <ID to type table> <input...>

# Start plot
pdf(args[3])

# Read conversion table
idToType <- read.table(args[4], sep = "\t", header = T)

# Initialise
file <- args[5]
exons <- read.table(file, sep = "\t", header = T)
colnames(exons) <- c("parent",
                     "numberExons",
                     "lengthExons",
                     "totalCovered",
                     "fractionCovered",
                     "meanDepth")
exonData <- merge(exons, idToType, by.x = "parent", by.y = "ID")
pathSplit <- unlist(strsplit(file, "/"))
sample <- sub("\\.tsv$", "", pathSplit[length(pathSplit)])
if (identical(sample, character(0))) {
  sample <- pathSplit
}
print(sample)
ggplot(exonData, aes(x = meanDepth/median(meanDepth))) +
  geom_histogram(binwidth = 0.01) +
  xlim(0,3) +
  xlab("Normalised depth") +
  facet_grid(type ~ .) +
  ggtitle("Histogram for depth for all types", subtitle = sample)
ggplot(exonData, aes(x = fractionCovered)) +
  geom_histogram(binwidth = 0.01) +
  xlab("Breadth") +
  scale_y_log10() +
  facet_grid(type ~ .) +
  ggtitle("Histogram of breadth for all types", subtitle = sample)
depthData <- data.frame(parent = exonData$parent,
                        type = exonData$type,
                        numExons = exonData$numberExons,
                        lenExons = exonData$lengthExons,
                        sample = exonData$meanDepth/median(exonData$meanDepth))
colnames(depthData)[ncol(depthData)] <- sample
breadthData <- data.frame(parent = exonData$parent,
                          type = exonData$type,
                          numExons = exonData$numberExons,
                          lenExons = exonData$lengthExons,
                          sample = exonData$fractionCovered)
colnames(breadthData)[ncol(breadthData)] <- sample

for (file in args[c(6:length(args))]) {
  exons <- read.table(file, sep = "\t", header = T)
  colnames(exons) <- c("parent",
                       "numberExons",
                       "lengthExons",
                       "totalCovered",
                       "fractionCovered",
                       "meanDepth")
  exonData <- merge(exons, idToType, by.x = "parent", by.y = "ID")
  pathSplit <- unlist(strsplit(file, "/"))
  sample <- sub("\\.tsv$", "", pathSplit[length(pathSplit)])
  if (identical(sample, character(0))) {
    sample <- pathSplit
  }
  print(sample)
  print(ggplot(exonData, aes(x = meanDepth/median(meanDepth))) +
          geom_histogram(binwidth = 0.01) +
          xlim(0,3) +
          xlab("Normalised depth") +
          facet_grid(type ~ .) +
          ggtitle("Histogram for depth for all types", subtitle = sample))
  print(ggplot(exonData, aes(x = fractionCovered)) +
          geom_histogram(binwidth = 0.01) +
          xlab("Breadth") +
          scale_y_log10() +
          facet_grid(type ~ .) +
          ggtitle("Histogram of breadth for all types", subtitle = sample))
  depthData[sample] <- exonData$meanDepth/median(exonData$meanDepth)
  breadthData[sample] <- exonData$fractionCovered
}
dev.off()

write.table(depthData, file = args[1], sep = "\t", quote = F, row.names = F)
write.table(breadthData, file = args[2], sep = "\t", quote = F, row.names = F)

