setwd("D:\\workspace\\bioinformatics\\BIOLF449\\tp4\\input\\data")
#read manifest file MANIFEST.txt
t <- read.table('MANIFEST.txt', header = T)

#Read annotations file
ensembl <-read.table('..\\mart_export.csv', header = T,sep=",")
allAnnotations <- data.frame( geneName = ensembl$Ensembl.Gene.ID,
                              start = ensembl$Gene.Start..bp., end = ensembl$Gene.End..bp.)

#this method returns a vector with the gene sizes for given gene names looked up in the 
#annotations file
#because in some cases the size of the same gene is variable, I was instructed to take
#the delta between max(end) and min(start)
findGeneSizeByGeneNameInAnnotations <- function(annotations, geneNames){
  communGenes <- annotations[annotations$geneName %in% geneNames, ]
  #communGenes contains duplicated data with sizes which don't always match
  #I was advised to take the union between min(start) and max(end)
  max(communGenes$end) - min(communGenes$start)
}

#generates fpkm array
fpkm <- function(df){
  librarySize <- sum(df$fragmentCount)/ 1e6
  regionSize <- df$size/ 1e3
  fragmentCount <- df$fragmentCount
  fragmentCount/(librarySize * regionSize)
}

#computes factor for generating tpm
calculateTFactor <- function(df){
  fragmentCount <- df$fragmentCount
  regionSize <- df$size/ 1e3
  
  sum(fragmentCount/regionSize) 
}

#generates tpm array
tpm <- function (df, factor) {
  regionSize <- df$size/ 1e3
  fragmentCount <- df$fragmentCount
  
  fragmentCount/(regionSize * factor)
}


"
createCountsData method creates a data frame with the following structure:
fragmentCount        geneName   size       fpkm          tpm
1          2016 ENSG00000000003  12882 3.15720689 2.450846e-05
2            23 ENSG00000000005  15083 0.03076351 2.388080e-07
3          1770 ENSG00000000419  23688 1.50744224 1.170183e-05
4           777 ENSG00000000457  44636 0.35118144 2.726117e-06
5           617 ENSG00000000460 192073 0.06480592 5.030690e-07
6           617 ENSG00000000938  23213 0.53622826 4.162580e-06

"
createCountsData <- function(input, allAnnotations){
  counts <- input
  #rename V2 to fragmentCount
  colnames(counts)[2] <- "fragmentCount"
  
  #add a colum with the genes names with substring after . removed
  counts$geneName <- sapply(strsplit(as.vector(counts$V1), "\\."), `[[`, 1)
  
  #remove original gene name column V1 as we will be using gene name
  counts <- counts[,-1]
  
  #select only data which starts with ENSG
  counts <- counts[grep('ENSG', counts$geneName),]
  
  #order counts by gene name in order to match the same elements in emsemble dataset
  counts <- counts[order(counts$geneName),]
  
  communGeneNames <-intersect(counts$geneName, allAnnotations$geneName)
  
  counts <- counts[counts$geneName %in% communGeneNames, ]

  #Add size column
  counts$size <-findGeneSizeByGeneNameInAnnotations(allAnnotations, communGeneNames)
  
  #Add fpkm column
  counts$fpkm <-fpkm(counts)
  
  #Add tpm column
  tFactor <- calculateTFactor(counts)
  counts$tpm <- tpm(counts, tFactor)
  
  counts
}

#Compute the FPKM and TPM for all samples
"
The method below will return the result in a data frame where I store fpkm results for 
sample i in fpkm<i> column, similarly for tpkm<i>
> head(allCounts)
       fpkm1       fpkm 2     fpkm 3       fpkm 4      fpkm 5       fpkm 6       tpm 7   counts$tpm counts$fpkm
1 3.15720689 2.450846e-05 0.59166346 5.359428e-06 2.612277020 2.069057e-05 1.216907255 8.028549e-06  0.59644494
2 0.03076351 2.388080e-07 0.02180334 1.974998e-07 0.008879913 7.033347e-08 0.003144717 2.074728e-08  0.01242459
3 1.50744224 1.170183e-05 1.85868452 1.683640e-05 1.530863700 1.212523e-05 1.239457289 8.177323e-06  1.15390247
4 0.351

"
calculateAllFPKM_TPM <- function(filenames) {
  sample1 <- read.table(as.character(filenames[1]), header = F)
  counts1 <- createCountsData(sample1, allAnnotations)
  #all counts are initialized with first sample
  allCounts <- data.frame(fpkm1 = counts1$fpkm, tpm1 = counts1$tpm)
  
  #read in a loop all remaining files
  for(i in 2:length(filenames)) {
    sample <- read.table(as.character(filenames[i]), header = F)
    
    counts <- createCountsData(sample, allAnnotations)
    
    allCounts <- cbind(allCounts, counts$fpkm)
    allCounts <- cbind(allCounts, counts$tpm)

    
    colnames(allCounts)[i] <- paste('fpkm', i)
    colnames(allCounts)[i + 1] <- paste('tpm', i + 1)
  }
  allCounts
}

allCounts <- calculateAllFPKM_TPM(t$filename)


#Draw a MA plot sample 1 vs sample 2 using the FPKM

sample1<- read.table(as.character(t$filename[1]), header = F)
sample2<- read.table(as.character(t$filename[2]), header = F)

counts1 <- createCountsData(sample1, allAnnotations)
counts2 <- createCountsData(sample2, allAnnotations)


meanExpression <- function( A, B) {
  1/2 * (A + B)
} 

log2FC <- function( A, B) {
  log2(A) - log2(B)
} 


plotMa <- function(input1, input2) {
  pdf("maplot.pdf")
  x <- meanExpression(input1, input2)
  y <- log2FC(input1, input2)
  plot(x, y, main = "Ma plot sample 1 vs sample 2", xlab = "mean expression", ylab = "log fold change")
  
  plot(x, y, main = "Ma plot sample 1 vs sample 2 logaritmic scale", xlab = "mean expression", ylab = "log fold change",  log="x")
  dev.off()
}


plotMa(counts1$fpkm, counts2$fpkm )


#Compute the FC, Log ratio and Pvalue by comparing the 10 first samples vs the 10 last (MANIFEST table order)

"
This function returns an array with the fpkm values calculated from given filename
"
fpkmFromSample <- function(filename){
  sample <- read.table(as.character(filename), header = F)
  counts <- createCountsData(sample, allAnnotations)
  counts$fpkm
}

"
population standard deviation = square root (1/n * sum (xi - xavg)^2 )
"
populationStandardDeviation <- function(v) {
  average <- sum(v)/length(v)
  sqrt(sum((v- average)^2)/(length(v)))
}


"
sample standard deviation = square root (1/n-1 * sum (xi - xavg)^2 )
"
sampleStandardDeviation <- function(v) {
  average <- sum(v)/length(v)
  sqrt(sum((v- average)^2)/(length(v) -1))
}


"
Puts together the FPKM values for each sample file as a separate column.
Sample files are looked up in the manifest file from start index to end index.

This method returns a data frame with endIndex-startIndex columns, each
column storing the fpkm values for the sample.
"
sampleGroupFPKM <- function(filenames, startIndex, endIndex) {
  #init result data frame with the first sample
  fpkmBySample <- data.frame(fpkmFromSample(filenames[startIndex]))
  #read in a loop all remaining files, for each add the column (using cbind)
  for(i in (startIndex + 1):endIndex) {
    fpkmBySample <- cbind(fpkmBySample, fpkmFromSample(filenames[i]))
  }
  fpkmBySample
}

"
Returns the mean , sample standard deviation and the sample size for given group
        mean         sd size
1 2.12038246 1.35165583    3
2 0.02048225 0.01100145    3
3 1.63233015 0.19637812    3
4 0.40758162 0.09295069    3
"
groupStatistics <- function(group) {
  sampleGroupMean <- apply(group, 1 , mean)
  sampleGroupStandardDeviation <- apply(group, 1 , sampleStandardDeviation)
  sampleGroupSize <-apply(group, 1 , length)
  sampleGroupStatistics <- data.frame(mean = sampleGroupMean, sd = sampleGroupStandardDeviation, size = sampleGroupSize)
  sampleGroupStatistics
}

"
Compute the mean, standard deviation and size for the 2 samples (first 10 files and last 10 files)
head(first_10_samples)
        mean         sd size
1 2.12038246 1.35165583    10
2 0.02048225 0.01100145    10
3 1.63233015 0.19637812    10
4 0.40758162 0.09295069    10
"
first_10_samples <- groupStatistics(sampleGroupFPKM(t$filename, 1, 10))
last_10_samples <- groupStatistics(sampleGroupFPKM(t$filename, length(t$filename) - 10, length(t$filename)))


FC <- function( a, b) {
  ifelse(a > b , a/b , - b/a)
}

"
Compute fold change and log ratio between the average values of first group sample (first 10 files)
and the second group sample (last 10 files)
"
samplesFC <- FC(first_10_samples$mean, last_10_samples$mean)
samplesLogRatio <- log2FC(first_10_samples$mean, last_10_samples$mean)

"
one sample t test formula
t = (sample mean - population mean)/(standard deviation / square root(sample size))
"
calculateTValueOneSample <- function (v, populationMean){
  sample.mean <- sum(v)/length(v)
  standard.deviation <- sampleStandardDeviation(v)
  (sample.mean - populationMean)/(standard.deviation/ sqrt(length(v)))
}


"
Wikipedia suggests using Welch's t-test for testing the hypothesis that 2 populations
have equal means : https://en.wikipedia.org/wiki/Welch%27s_t-test

t = (average1 - average2)/square root( sd1 ^2 / n1 + sd2^2 / n2)

"
calculateTValueWelch <- function (sample1, sample2) {
  (sample1$mean - sample2$mean)/sqrt( sample1$sd^2 / sample1$size + sample2$sd^2 / sample2$size  )
}

calculatePValue <- function (t.value, degreesOfFreedom) {
  2*pt(-abs(t.value), df=degreesOfFreedom )
}


tValues <- calculateTValueWelch(first_10_samples, last_10_samples)
degreesOfFreedom <- last_10_samples$size[1] -1 #sample size -1

pValues <- calculatePValue(tValues, degreesOfFreedom)


"
Install library to get q value
"
source("https://bioconductor.org/biocLite.R")
biocLite("qvalue")
browseVignettes("qvalue")
library(qvalue)

qobj <- qvalue(p = pValues)
qValues <- qobj$qvalues

# we aggregate together the pvalues we obtained for the second sample group with the logFC 
# and the q value calculated from p value
volcanoData <- data.frame( pValues = pValues, logFC = samplesLogRatio, qValues = qValues)

significant <- volcanoData[volcanoData$pValues < 0.05,]
nonSignificant<-volcanoData[volcanoData$pValues >= 0.05,]

pdf("volcano_plot.pdf")

plot(nonSignificant$logFC, -log10(nonSignificant$pValues), main = "FPKM first 10 samples vs last 10 samples (pvalue)", xlab = "log2 Fold Change", ylab = "-log10 pvalue" , col="black")
points(significant$logFC, -log10(significant$pValues), col="red")


plot(nonSignificant$logFC, -log10(nonSignificant$qValues), main = "FPKM first 10 samples vs last 10 samples(qvalue)", xlab = "log2 Fold Change", ylab = "-log10 qValues" , col="black")
points(significant$logFC, -log10(significant$qValues), col="red")
dev.off()

