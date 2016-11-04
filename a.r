setwd("D:\\workspace\\bioinformatics\\BIOLF449\\tp4\\input")
ensembl <-read.table('mart_export.csv', header = T,sep=",")

#Compute the average number of transcripts per gene
frequency <-table(as.vector(ensembl$Ensembl.Gene.ID))

average <- function(v){
  sum(v)/length(v)
}

averageApproach1 <- average(frequency)
averageApproach2 <- mean(frequency)
all (averageApproach1 == averageApproach2)

#But if we consider that input data gives us duplicated values for the pair geneId - transcriptId 
#due to other factors, we only want to keep the unique combinations
#For instance:
#ENSG00000100023,ENST00000335025,22,21652270,21700015,1,54.20,NM_148176,PPIL2
#ENSG00000100023,ENST00000335025,22,21652270,21700015,1,54.20,NM_148175,PPIL2
#ENSG00000100023,ENST00000398831,22,21652270,21700015,1,54.20,,PPIL2
#should only count ENST00000335025 and ENST00000398831 for ENSG00000100023

noDuplicatesData <-unique(ensembl[c("Ensembl.Gene.ID", "Ensembl.Transcript.ID")])
frequencyNoDuplicates<-table(as.vector(noDuplicatesData$Ensembl.Gene.ID))
averageNoDuplicates <- mean(frequencyNoDuplicates)
#What is the maximum, minimum and average number of transcripts per gene?

min(frequency)
max(frequency)
mean(frequency)

min(frequencyNoDuplicates)
max(frequencyNoDuplicates)
mean(frequencyNoDuplicates)

#Plot the number of transcript per gene distribution (histogram)
frequencyHistogram = as.vector(frequencyNoDuplicates)

pdf("a.pdf")

hist(frequencyHistogram, col = 'red', nclass =200, main ="Histogram for number of transcripts per gene distribution",
     xlab = 'Number of transcripts', ylab = 'Frequency')

hist(frequencyHistogram, col = 'red', nclass =300, main ="Histogram for number of transcripts per gene distribution xlim 40",
     xlab = 'Number of transcripts', ylab = 'Frequency', xlim = c(0,40))

dev.off()


