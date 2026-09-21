install.packages("genetics")
library(genetics)

setwd("/Users/ulasawczuk/Downloads/szkola/BSG/lab1")

dataSNP <- read.table("TSICHR22RAW.raw", header = TRUE)
head(dataSNP)

# only 7th column and further

only_gen <- dataSNP[, 7:ncol(dataSNP)]
head(only_gen)



# a) how many variants and how many missing

num_variants <- ncol(only_gen)
total_data_points <- nrow(only_gen) * ncol(only_gen)
missing_data <- sum(is.na(only_gen))
perc_missing <- (missing_data / total_data_points) * 100

cat("Total Variants:", num_variants, "\n")
cat("Percentage missing data:", perc_missing, "%\n\n")

# ANSWER TO QUESTION: Numebr of total variants: 20649, Percentage of missing data: 0.1986518 %



# b) percentage of monomorphic variants, remove them from database

# we just check for variance as monomorphic ones dont have any
monomorphic <- sapply(only_gen, function(x) { var(x, na.rm = TRUE) == 0 })
perc_monomorphic <- (sum(monomorphic) / num_variants) * 100

only_polymorphic <- only_gen[, !monomorphic]
num_remaining <- ncol(only_polymorphic)

cat("Percentage of monomorphic variants:", perc_monomorphic, "%\n")
cat("Variants remaining after exclusion:", num_remaining, "\n\n")

# ANSWER TO QUESTION: Percentage of monomorphic variants: 11.45818 %, Variants remaining after exclusion: 18283 



# c) genotype counts, allele counts and MAF of polymorphism rs8138488_C

rs81 <- only_polymorphic$rs8138488_C
gen_counts <- table(rs81)
cat("Genotype counts for rs8138488_C:\n")
print(gen_counts)

# coding is 0=AA, 1=AB, 2=BB

n0 <- sum(rs81 == 0, na.rm = TRUE)
n1 <- sum(rs81 == 1, na.rm = TRUE)
n2 <- sum(rs81 == 2, na.rm = TRUE)
total_alleles <- 2 * (n0 + n1 + n2)

count_A <- (2 * n0) + n1
count_B <- (2 * n2) + n1

maf_rs <- min(count_A, count_B) / total_alleles

cat("Allele count: A:", count_A, " B:", count_B, "%\n")
cat("MAF:", maf_rs, "\n\n")

# ANSWER TO QUESTION: counts for genotype rs8138488_C: rs81 0:41  1:47  2:14, Allele count: A: 129  B: 75 %,
# MAF: 0.3676471 



# d) MAF for all + histogram

calculate_maf <- function(x) {
  n0 <- sum(x == 0, na.rm = TRUE)
  n1 <- sum(x == 1, na.rm = TRUE)
  n2 <- sum(x == 2, na.rm = TRUE)
  total <- 2 * (n0 + n1 + n2)
  if(total == 0) return(NA)
  p_A <- ((2 * n0) + n1) / total
  p_B <- ((2 * n2) + n1) / total
  return(min(p_A, p_B))
}

mafs <- sapply(only_polymorphic, calculate_maf)
hist(mafs, main="Histogram of MAF for SNPs", xlab="Minor Allele Frequency", col="blue")

perc_below_05 <- mean(mafs < 0.05, na.rm = TRUE) * 100
perc_below_01 <- mean(mafs < 0.01, na.rm = TRUE) * 100

cat("Percentage of markers with MAF < 0.05:", perc_below_05, "%\n")
cat("Percentage of markers with MAF < 0.01:", perc_below_01, "%\n\n")

## ANSWER TO QUESTION: It does not follow a uniform distribution, it is very right-skewed. In human genetics it is 
# expectable to observe variants with low MAF as most mutations are rare. Percantage below 0.05: 14.22633 %, 
# below 0.01: 4.731171 %.


# e) heterozygosity H0 and histogram

calculate_H0 <- function(x) {
  n1 <- sum(x == 1, na.rm = TRUE)
  total_valid <- sum(!is.na(x))
  return(n1 / total_valid)
}

H0 <- sapply(only_polymorphic, calculate_H0)
hist(H0, main="Histogram of Observed Heterozygosity (SNPs)", xlab="H0", col="green")

## ANSWER TO QUESTION: The theoretical range for observed heterozygosity mathematically goes up to 1 (if an entire 
# sampled population consisted exclusively of heterozygous individuals), but realistically it matches expected 
# ranges, capping at approximately 0.5



# f) expected heterozygosity 

calculate_He <- function(maf) { return(1 - (maf^2 + (1 - maf)^2)) }
He <- sapply(mafs, calculate_He)

hist(He, main="Histogram of Expected Heterozygosity (SNPs)", xlab="He", col="orange")
cat("Average Expected Heterozygosity (He) for SNPs:", mean(He, na.rm=TRUE), "\n\n")

## ANSWER TO QUESTION: the theoretical maximum for expected heterozygosity is exactly 0.5. the formual for expected 
# heterozygosity yields 0.5 when the two allele frequencies are equal, the average He is 0.3115841 
