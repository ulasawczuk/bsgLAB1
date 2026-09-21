library(genetics)

setwd("/Users/ulasawczuk/Downloads/szkola/BSG/lec1")

dataSNP <- read.table("TSICHR22RAW.raw", header = TRUE)
head(dataSNP)

# only 7th column and further

only_gen <- dataSNP[, 7:ncol(dataSNP)]
head(only_gen)

# a) how many variants and how many missing
