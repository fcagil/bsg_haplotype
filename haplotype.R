rm(list=ls())
#install.packages("haplo.stats")
library(haplo.stats)
library(genetics)
APOE <- read.table("APOE.dat", header = TRUE)

#2
#How many individuals and how many SNPs are there in the database? What percentage of
#the data is missing?

#number of individuals
nrow <- nrow(APOE)
ncol <- ncol(APOE)-1

nmis <- function(x) {
  y <- sum(is.na(x))
  return(y)
}

nmiss <- apply(APOE,2,nmis)
#number of missings
sum(nmiss)


#3
#Assuming all SNPs are bi-allelic, how many haplotypes can theoretically be found for this
#data set?
2**162


#4
#Estimate haplotype frequencies using the haplo.stats package (set the minimum posterior
#probability to 0.001). How many haplotypes do you find? List the estimated probabilities in
#decreasing order. Which haplotype number is the most common?

alleles.matrix <- APOE
alleles.matrix <- alleles.matrix[,-1]
names <- colnames(alleles.matrix)
#alleles.matrix <- (gsub("/", "", alleles.matrix))
new.alleles.matrix <- matrix(0, nrow = nrow(alleles.matrix), ncol = ncol(alleles.matrix)*2)


j <- 1
for(i in 1:ncol){
  new.alleles.matrix[,j] <-   substr(alleles.matrix[,i],1,1)
  new.alleles.matrix[,j+1] <- substr(alleles.matrix[,i],3,3)
  j <- j+2
}

alleles.matrix[1:5,1:5]
new.alleles.matrix[1:5,1:5]
nrow(new.alleles.matrix)
ncol(new.alleles.matrix)

length(names)

table(new.alleles.matrix[,6])

nrow(new.alleles.matrix)

Haplo.Res <- haplo.em(new.alleles.matrix,locus.label=names,control=haplo.em.control(min.posterior=0.001))
Haplo.Res

attributes(Haplo.Res)
nrow(Haplo.Res$haplotype) #How many haplotypes do you find?

hprob <- Haplo.Res$hap.prob
hprob

df = data.frame(c(1:31),hprob)
df
df <- df[order(hprob, decreasing = TRUE),]
rownames(df) <- 1:nrow(df)
library(xtable)
print(xtable(df, digits=6), type = "latex")

#5
#Is the haplotypic constitution of any of the individuals in the database ambiguous or uncertain?
#For how many? What is the most likely haplotypic constitution of individual NA20763? (identify
#the constitution by the corresponding haplotype numbers).

hist(Haplo.Res$post, main="Histogram of posterior probabilities", xlab="Posterior probability", ylab="Counts", labels = TRUE, ylim=c(0,110))

sum(Haplo.Res$nreps>1)

table(Haplo.Res$nreps)

which(APOE$id=='NA20763')
Haplo.Res$indx.subj[Haplo.Res$indx.subj==which(APOE$id=='NA20763')] #meaning that we have 3 hp for it


NA20763_proba <- Haplo.Res$post[Haplo.Res$indx.subj==which(APOE$id=='NA20763')]
NA20763_hp <-Haplo.Res$hap1code[Haplo.Res$indx.subj==which(APOE$id=='NA20763')]
NA20763_df = data.frame(NA20763_proba,NA20763_hp)
print(xtable(NA20763_df, digits=6), type = "latex")

Haplo.Res$indx.subj
Haplo.Res$hap1code
 
#6
#Suppose we would delete polymorphism rs374311741 from the database prior to haplotype estimation. 
#Would this affect the results obtained?

APOE$rs374311741

#7
#Remove all genetic variants that have a minor allele frequency below 0.10 from the database,
#and re-run haplo.em. How does this affect the number of haplotypes?
alleles.matrix[1:5,1:5]
maf <- function(x){
  x <- genotype(x)
  out <- summary(x)
  af1 <- min(out$allele.freq[,2],na.rm=TRUE)
  af1[af1==1] <- 0 
  return(af1)
}
maf.per.snp <- apply(alleles.matrix,2,maf)
#find the position of maf<0.1
maf.per.snp
list <- which(maf.per.snp < 0.1)
list
#create a new data without those columns
alleles.matrix.maf <- alleles.matrix[,-list]


names <- colnames(alleles.matrix.maf)
new.alleles.matrix <- matrix(0, nrow = nrow(alleles.matrix.maf), ncol = ncol(alleles.matrix.maf)*2)

ncol(alleles.matrix.maf)
j <- 1
for(i in 1:ncol(alleles.matrix.maf)){
  new.alleles.matrix[,j] <-   substr(alleles.matrix.maf[,i],1,1)
  new.alleles.matrix[,j+1] <- substr(alleles.matrix.maf[,i],3,3)
  j <- j+2
}

alleles.matrix.maf[1:5,1:5]
new.alleles.matrix[1:5,1:5]
nrow(new.alleles.matrix)
ncol(new.alleles.matrix)


Haplo.Res <- haplo.em(new.alleles.matrix,locus.label=names,control=haplo.em.control(min.posterior=0.001))
Haplo.Res


nrow(Haplo.Res$haplotype)
Haplo.Res$haplotype



#8
#We could consider the newly created haplotypes in our last run of haplo.em as the alleles
#of a new superlocus. Which is, under the assumption of Hardy-Weinberg equilibrium, the most
#likely genotype at this new locus? What is the probability of this genotype? Which genotype is
#the second most likely, and what is its probability?
Haplo.Res
hprob <- Haplo.Res$hap.prob


df <- data.frame(c("A", "B", "C", "D", "E", "F", "G", "H"), hprob)
print(xtable(df, digits=6), type = "latex", include.rownames = FALSE)



m <- hprob%*%t(hprob)
m
# which(m == max(m), arr.ind = TRUE)


# Which genotype is the second most likely, and what is its probability?
m.sec <- hprob%*%t(hprob)
m.sec[8,8] <- 0
m.sec
which(m.sec == max(m.sec), arr.ind = TRUE)
max(m.sec)
# 
