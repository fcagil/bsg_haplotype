install.packages("haplo.stats")
library(haplo.stats)
library(genetics)
APOE <- read.table("APOE.dat", header = TRUE)

#2
#How many individuals and how many SNPs are there in the database? What percentage of
#the data is missing?

#number of individuals
nrow <- nrow(APOE)-1
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
#4
#Estimate haplotype frequencies using the haplo.stats package (set the minimum posterior
#probability to 0.001). How many haplotypes do you find? List the estimated probabilities in
#decreasing order. Which haplotype number is the most common?

#I just converted the data into the format that we have in the example script.
alleles.matrix <- t(APOE)
names <- alleles.matrix[1,]
alleles.matrix <- alleles.matrix[-1,]
alleles.matrix <- (gsub("/", "", alleles.matrix))
new.alleles.matrix <- matrix(0, nrow = nrow(alleles.matrix), ncol = ncol(alleles.matrix)*2)

j <- 1
for(i in 1:ncol(APOE)){
  new.alleles.matrix[,j] <-   substr(alleles.matrix[,i],1,1)
  new.alleles.matrix[,j+1] <- substr(alleles.matrix[,i],2,2)
  j <- j+2
}
rownames(alleles.matrix) <-colnames(APOE)


Haplo.Res <- haplo.em(new.alleles.matrix,locus.label=names,control=haplo.em.control(min.posterior=0.001))
Haplo.Res

##not sure if next part is true or not##

nHaploPossible <- 2^nrow
nHaploPossible

#5
#Is the haplotypic constitution of any of the individuals in the database ambiguous or uncertain?
#For how many? What is the most likely haplotypic constitution of individual NA20763? (identify
#the constitution by the corresponding haplotype numbers).

#6
#Suppose we would delete polymorphism rs374311741 from the database prior to haplotype estimation. 
#Would this affect the results obtained?