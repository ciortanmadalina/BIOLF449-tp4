setwd("D:\\workspace\\bioinformatics\\BIOLF449\\tp4\\input\\data")
#read manifest file MANIFEST.txt
t <- read.table('MANIFEST.txt', header = T)

#Read annotations file
ensembl <-read.table('..\\mart_export.csv', header = T,sep=",")
allAnnotations <- data.frame( geneName = ensembl$Ensembl.Gene.ID, size = ensembl$Gene.End..bp - ensembl$Gene.Start..bp)

#this method generated the read size column
readSize <- function(v, commun){
  filteredIds <- v[v$geneName %in% commun, ]
  removeDuplicates <- filteredIds[match(unique(filteredIds$geneName), filteredIds$geneName),]
  #sort by gene name and return
  removeDuplicates[order(removeDuplicates$geneName), 'size']
}

#generates fpkm column
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

#generates tpm column
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
  #rename V2 to counts
  colnames(counts)[2] <- "fragmentCount"
  
  #add a colum with the genes names with substring after . removed
  counts$geneName <- sapply(strsplit(as.vector(counts$V1), "\\."), `[[`, 1)
  
  #remove original gene name column V1 as we will be using gene name
  counts <- counts[,-1]
  
  #select only data which starts with ENSG
  counts <- counts[grep('ENSG', counts$geneName),]
  
  #order counts by gene name in order to match the same elements in emsemble dataset
  counts <- counts[order(counts$geneName),]
  
  commun <-intersect(counts$geneName, allAnnotations$geneName)
  
  counts <- counts[counts$geneName %in% commun, ]

  #Add size column
  counts$size <-readSize(allAnnotations, commun)
  
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

aggregateSamples <- function(filenames, startIndex, endIndex) {
  sample1 <- read.table(as.character(filenames[startIndex]), header = F)
  counts1 <- createCountsData(sample1, allAnnotations)
  
  sumCounts <- counts1$fpkm
  
  #read in a loop all remaining files
  for(i in (startIndex + 1):endIndex) {
    sample <- read.table(as.character(filenames[i]), header = F)
    
    counts <- createCountsData(sample, allAnnotations)
    sumCounts <- sumCounts + counts$fpkm
  }
  sumCounts
}

sum_sample_first_10 <- aggregateSamples(t$filename, 1, 10)
sum_sample_last_10 <- aggregateSamples(t$filename, length(t$filename) - 10, length(t$filename))

FC <- function( a, b) {
  ifelse(a > b , a/b , - b/a)
}

samplesFC <- FC(sum_sample_first_10, sum_sample_last_10)
samplesLogRatio <- log2FC(sum_sample_first_10, sum_sample_last_10)

pvalue <- cor.test(sum_sample_first_10,sum_sample_last_10)


# Draw a volcano plot for the same groups (Pvalue using t-test)



