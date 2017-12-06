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
monomorphic / total_vars * 100 # 96.2 %
# 2) Remove all monomorphic SNPs from the data bases. How many variants remain in the
# database? 
heteromorphic <- Filter(function(y)(length(unique(y))>1), X)
ncol(heteromorphic) # 3035
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
# add chi p value column
geno_transposed$chi_p <- NA

# calc chi values
for(i in 1:dim(geno_transposed)[1]){
  num_vector <- as.numeric(geno_transposed[i, 1:3])
  names(num_vector) <- c("AA", "AB", "BB")
  chiStatistics <- HWChisq(num_vector,cc=0)
  geno_transposed$chi_p[i] <- chiStatistics[2]
  geno_transposed$allelle_f[i] <- chiStatistics[4]
}
# 2) How many SNPs are significant (use ?? = 0.05)?
nrow(subset(geno_transposed, geno_transposed$chi_p < 0.05))


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
# exact significant. 126.
# But compared to 3035 non-monomorfic markers it's insignificant difference, so in general results are consistent.
#Chi test is closer to expected amount out of equilibrium


# 5) Apply a likelihood ratio test for Hardy-Weinberg equilibrium to each SNP, using the HWLratio function.
# add HWLratio column
geno_transposed$HWLratio <- NA
# calc HWLratio values
for(i in 1:dim(geno_transposed)[1]){
  num_vector <- as.numeric(geno_transposed[i, 1:3])
  names(num_vector) <- c("AA", "AB", "BB")
  geno_transposed$HWLratio[i] <- as.numeric(HWLratio(num_vector)[1])
}

# How many SNPs are significant (use ?? = 0.05). Is the result consistent with the chi-square test?
geno_transposed$HWLratio_p <- NA
for(i in 1:dim(geno_transposed)[1]){
  num_vector <- as.numeric(geno_transposed[i, 1:3])
  names(num_vector) <- c("AA", "AB", "BB")
  geno_transposed$HWLratio_p[i] <- as.numeric(HWLratio(num_vector)[3])
}
nrow(subset(geno_transposed, geno_transposed$HWLratio_p < 0.05))
# 145, result is almost consistent with the chi-square test of 160
# delete this useless column
geno_transposed$HWLratio_p <- NULL

# 6) Apply a permutation test for Hardy-Weinberg equilibrium to the first 10 SNPs, using the
# classical chi-square test (without continuity correction) as a test statistic. 
# add HWPerm column
geno_transposed$HWPerm_p <- NA
# calc HWPerm values
for(i in 1:10){
  num_vector <- as.numeric(geno_transposed[i, 1:3])
  names(num_vector) <- c("AA", "AB", "BB")
  geno_transposed$HWPerm_p[i] <- HWPerm(num_vector, nperm = 1000)[2]
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



# 9) Make a histogram of the p-values obtained in the chi-square test
# chi vector stores p-val from HWChisq
pvalues <- unlist(geno_transposed$chi_p)
hist(pvalues, breaks=20)

# 9) What distribution would you expect if HWE would hold for the data set? What distribution do you observe?
# if pval < 0.05 then it's unlikely that genotype differences are due to chance.
#Beta distribution could be possible one. Null distribution or uniform?

# 9) make a Q-Q plot of the p values obtained in the chi-square test against the quantiles of the distribution that you consider relevant. 
qqnorm(pvalues)
qqline(pvalues,  datax = FALSE, distribution = qlogis() ,probs = c(0.25, 0.75), qtype = 7)
# 9) What is your conclusion?


# 10) Imagine that for a particular marker the counts of the two homozygotes are accidentally interchanged. Would this affect the statistical tests for HWE?
test_geno <- geno_transposed[1:10,]
tmpAA <- test_geno$AA
test_geno$AA <- test_geno$BB
test_geno$BB <- tmpAA
test_geno$HWExact <- HWExactStats(test_geno[,1:3])
for(i in 1:10){
  num_vector <- as.numeric(test_geno[i, 1:3])
  names(num_vector) <- c("AA", "AB", "BB")
  test_geno$HWPerm[i] <- HWPerm(num_vector, nperm = 1000, verbose = FALSE)[2]
  test_geno$HWLratio[i] <- HWLratio(num_vector, verbose = FALSE)
  chiStats <- HWChisq(num_vector,cc=0, verbose = FALSE)
  test_geno$pval[i] = chiStats[2]
  test_geno$chisq[i] = chiStats[1]
  test_geno$D[i] = chiStats[3]
  test_geno$p[i] = chiStats[4]
  test_geno$f[i] = chiStats[5]
  test_geno$expected[i] = chiStats[6]
}
test_geno[1:3,]
geno_transposed[1:3,]
#D Doesn't seem like it affects something -> all statistics counted previously doesn't affect results, as well as all outputs from HWChisq(). It makes sense, as we have just replaced labels for A and B, in a nutshell.

# 11)  Compute the inbreeding coefficient ( ^f) for each SNP, and make a histogram of ^f. You can
#use function HWf for this purpose. 
# add inbreeding coefficient column
geno_transposed$inbrd <- NA
# calc inbreeding coefficient values
for(i in 1:dim(geno_transposed)[1]){
  num_vector <- as.numeric(geno_transposed[i, 1:3])
  names(num_vector) <- c("AA", "AB", "BB")
  geno_transposed$inbrd[i] <- HWf(num_vector)
}
# Give descriptive statistics (mean, standard deviation, etc) of ^f calculated over the set of SNPs. 
inbrd_mean <- mean(geno_transposed$inbrd) # mean -0.005368393
inbrd_sd   <- sd(geno_transposed$inbrd)   # sd 0.1372449
inbrd_min  <- min(geno_transposed$inbrd)  # -0.9067086
inbrd_max  <- max(geno_transposed$inbrd)  # 1
# What distribution do you think ^f follows? Use a probability plot to confirm your idea
# At first we thought that it follows normal distribution but it turns out that it is beta distribution:
# got help from https://www.statmethods.net/advgraphs/probability.html
x1 <- seq(-4,4,length=100)*inbrd_sd + inbrd_mean
hx1 <- dnorm(x1,inbrd_mean,inbrd_sd)
plot(x1, hx1)



# 12) Make a plot of the observed chi-square statistics agains the inbreeding coefficient ( ^f). 
#make chi value column
geno_transposed$chi_val <- NA
# calc chi values' p-values
for(i in 1:dim(geno_transposed)[1]){
  num_vector <- as.numeric(geno_transposed[i, 1:3])
  names(num_vector) <- c("AA", "AB", "BB")
  geno_transposed$chi_val[i] <- as.numeric(HWChisq(num_vector,cc=0)[1])
}

chi_mean <- mean(geno_transposed$chi_val) # mean 2.01789
chi_sd   <- sd(geno_transposed$chi_val)   # sd 11.49767
chi_min  <- min(geno_transposed$chi_val)  # min  4.070658e-05
chi_max  <- max(geno_transposed$chi_val)  # 107
# plot the distributions
x2 <- seq(-4,4,length=100)*chi_sd + chi_mean
hx2 <- dnorm(x2,chi_mean,chi_sd)
plot(x1, hx1) # normal distribution
plot(x2, hx2) # normal distribution
# What do you observe? Can you give an equation that relates the two statistics?
# both plots have exaclty the same shape but are different with regards to height and width


# 13) Make a chi-square probability plot of the observed chi-square statistics against their theoretical quantiles.

#Does the sample statistic follow a chi-square distribution?


# 14) Simulate SNPs under the assumption of Hardy-Weinberg equilibrium. 

#First, we need to get allele frequencies vector
allelef <- numeric(3035)
for(i in 1:dim(geno_transposed)[1]){
  aa <-geno_transposed[i, 1]
  ab <-geno_transposed[i,2]
  bb <-geno_transposed[i,3]
  allelef[i] = (aa+ab/2)/(aa+ab+bb)
}

# maybe use HwChisq allele frequency instead. Docs says it gives frequency of AA, but in my case it's 1-frequency(AA).. Which is correct? Let's use the library one
allele_lib <- unlist(geno_transposed$allelle_f)
simulated <- HWData(nm = 3035, n = 107, p=allele_lib)
simulated <- data.frame(simulated)
for(i in 1:dim(simulated)[1]){
  num_vector <- as.numeric(simulated[i, 1:3])
  names(num_vector) <- c("AA", "AB", "BB")
  chiStats <- HWChisq(num_vector,cc=0, verbose = FALSE)
  simulated$pval[i] = chiStats[2]
  simulated$chisq[i] = chiStats[1]
  simulated$D[i] = chiStats[3]
  simulated$p[i] = chiStats[4]
  simulated$f[i] = chiStats[5]
  simulated$expected[i] = chiStats[6]
}

simulated[1:10,]
qqplot(unlist(geno_transposed$HWExact), unlist(simulated$pval), ylab = "simulated", xlab="original")

# from the plot we can see that original and simulated dataset come from same distribution

# 15) Is genotyping error a problem?
# In general, the number of markers out of equilibrium is relatively small (126 or 160, compared to expected 150).
# But we've read that range 5-10% is being considered "aproaching significant", and there are 218 and 183 markers, 
# according to Chi test and exact test, which isn't significant increase. 
# Also, QQ plot for original and simlated data's statistics shows they belong to same distribtion.
# As relatively small amount of markers are out of equilibrium (compared to expected) we can assume it's due to genotyping error.
nrow(subset(geno_transposed, geno_transposed$chi < 0.1))
nrow(subset(geno_transposed, geno_transposed$HWExact < 0.1))
