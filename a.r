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

#What is the maximum, minimum and average number of transcripts per gene?

min(frequency)
max(frequency)
mean(frequency)


#Plot the number of transcript per gene distribution (histogram)
frequencyHistogram = as.vector(frequency)

pdf("a.pdf")

hist(frequencyHistogram, col = 'red', nclass =200, main ="Histogram for number of transcripts per gene distribution",
     xlab = 'Number of transcripts', ylab = 'Frequency')

hist(frequencyHistogram, col = 'red', nclass =300, main ="Histogram for number of transcripts per gene distribution xlim 40",
     xlab = 'Number of transcripts', ylab = 'Frequency', xlim = c(0,40))

dev.off()


