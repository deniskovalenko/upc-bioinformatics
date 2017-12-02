# 1) Load the data
library(genetics)
library(HardyWeinberg)
filename = "YRIChr1.rda"
load(filename)
#x <- X # for faster typing

# 2) How many individuals does the database contain? 
nrow(X)
# 2) What percentage of the variants is monomorphic?
total_vars <- sum(table(unlist(X)))
monomorphic <- table(unlist(X))[1] + table(unlist(X))[3]
monomorphic / total_vars * 100
# 2) Remove all monomorphic SNPs from the data bases. How many variants remain in the
# database? 
heteromorphic <- Filter(function(y)(length(unique(y))>1), X)
ncol(heteromorphic)
#Determine the genotype counts for these variants, and store them in matrix. 
geno <- data.frame(matrix(ncol = ncol(heteromorphic), nrow = 3))
colnames_geno <- colnames(heteromorphic)
colnames(geno) <- colnames_geno
rownames(geno) <- c("AA", "AB", "BB")

for(name in colnames(heteromorphic)){
  geno[1, name] <- dim(subset(X[name], X[name] == 0))[1]
  geno[2, name] <- dim(subset(X[name], X[name] == 1))[1]
  geno[3, name] <- dim(subset(X[name], X[name] == 2))[1]
}

# 2) Apply a chi-square test without continuity correction for 
# Hardy-Weinberg equilibrium to each SNP. 
# transform the geno matrix so that cols are AA, AB, BB and row variants (loci)
geno_transposed <- as.data.frame(t(geno))
# add chi column
geno_transposed$chi <- NA
# calc chi values
for(i in 1:dim(geno_transposed)[1]){
  num_vector <- as.numeric(geno_transposed[i, 1:3])
  names(num_vector) <- c("AA", "AB", "BB")
  geno_transposed$chi[i] <- HWChisq(num_vector,cc=0)[2]
}
# 2) How many SNPs are significant (use ?? = 0.05)?
nrow(subset(geno_transposed, geno_transposed$chi < 0.05))


# 3) How many markers of the remaining non-monomorphic markers would you expect to be out
# of equilibrium by the effect of chance alone?
nrow(geno_transposed) * 0.05


# 4) Apply an Exact test for Hardy-Weinberg equilibrium to each SNP. You can use function
# HWExactStats for fast computation. 
geno_transposed$HWExact <- HWExactStats(geno_transposed[,1:3])
# 4) How many SNPs are significant (use ?? = 0.05). Is the result
# consistent with the chi-square test?
nrow(subset(geno_transposed, geno_transposed$HWExact < 0.05))
# No, the results are a little different:
# chi significant: 160
# exact significant. 126


# 5) Apply a likelihood ratio test for Hardy-Weinberg equilibrium to each SNP, using the HWLratio function.
# add HWLratio column
geno_transposed$HWLratio <- NA
# calc chi values
for(i in 1:dim(geno_transposed)[1]){
  num_vector <- as.numeric(geno_transposed[i, 1:3])
  names(num_vector) <- c("AA", "AB", "BB")
  geno_transposed$HWLratio[i] <- HWLratio(num_vector)
}

# How many SNPs are significant (use ?? = 0.05). Is the result consistent with the chi-square
# test?
nrow(subset(geno_transposed, geno_transposed$HWLratio < 0.05))
# 99, result not consistent with the chi-square test of 160

# 6) Apply a permutation test for Hardy-Weinberg equilibrium to the first 10 SNPs, using the
# classical chi-square test (without continuity correction) as a test statistic. 
# add HWPerm column
geno_transposed$HWPerm <- NA
# calc HWPerm values
for(i in 1:10){
  num_vector <- as.numeric(geno_transposed[i, 1:3])
  names(num_vector) <- c("AA", "AB", "BB")
  print(num_vector)
  geno_transposed$HWPerm[i] <- HWPerm(num_vector, nperm = 1000)[2]
}

# 6) List the 10 p-values, together with the 10 p-values of the exact tests. Are the result consistent?
geno_transposed[1:10,]
# In this case the values are consistent

# 7)  Depict all SNPs simultaeneously in a ternary plot, and comment on your result (because many
# genotype counts repeat, you may use UniqueGenotypeCounts to speed up the computations)
SNPs_geno <-  geno_transposed[,1:3]
SNPs_geno <- UniqueGenotypeCounts(SNPs_geno)[1:3]
HWTernaryPlot(SNPs_geno)

# It seems as if most of the variants in the database follow the boundaries of Hardy-Weinberg equilibrium
# Also, we can notice that the BB-homozygotes don't dominate at all - they are a numerical minority in
# all genotype counts



# 8) Can you explain why half of the ternary diagram is empty?
A <- sum(geno_transposed[,1]) + 0.5*sum(geno_transposed[,2])
B <- sum(geno_transposed[,3]) + 0.5*sum(geno_transposed[,2])
A/B

# There are much more B-alleles in the dataset, 11.13 times more, to be exact. In no examples of the genotype counts
# do number of B-alleles surpass that of A-alleles











