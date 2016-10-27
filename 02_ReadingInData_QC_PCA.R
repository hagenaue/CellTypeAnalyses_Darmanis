#These code documents contain all of the code documenting the cell type analyses that we applied to Darmanis' single-cell RNAseq data to validate our methodology
#Megan Hagenauer and Alek Pankonin
#October 27, 2016

#**************************************

#Alek's code from the initial joining and processing of the Darmanis human cortex single-cell RNAseq data downloaded from GEO:

SampleInfo<-read.csv("SampleInfo_forR.csv", header=T, colClasses ="character")

#Loading all of the separate files representing the RNAseq data for each cell into the R workspace and creating an object representing each of them:
Rename<-list.files()#must remove everything except for darmanis files for the code to run this step
for (i in 1:length(Rename)) assign(paste("",Rename[i], sep = ""), read.table(Rename[i]))
rm(Rename)
rm(i)

library(plyr)
ls()


Allobjects<-ls()

#Taking all of the objects representing the RNASeq data for individual cells and combining them into a single data matrix. Note that the only reason this code works is because: 1) there were only 2 other objects in the workspace 2) the RNAseq data for each cell is in an identical format - same length, same genes, same order. Otherwise, we would have needed to use a join function instead of cbind.


JoinedDarmanis<-cbind(GSM1657871_1772078217.C03.csv[,2],GSM1657872_1772078217.C04.csv[,2])
colnames(JoinedDarmanis)[c(1:2)]<-c("GSM1657871_1772078217.C03.csv", "GSM1657872_1772078217.C04.csv")

for(i in 3:length(Allobjects)){
  
JoinedDarmanis<-cbind(JoinedDarmanis, (get(Allobjects[i])[,2]))
colnames(JoinedDarmanis)[i]<-(Allobjects)[i]

}

head(JoinedDarmanis)
JoinedDarmanis2<-data.frame(GSM1657871_1772078217.C03.csv[,1], JoinedDarmanis)
colnames(JoinedDarmanis2)[1]<-"GeneSymbols"
write.csv(JoinedDarmanis2, "JoinedDarmanis.csv")

ColnamesForSignalData<-data.frame(colnames(JoinedDarmanis), stringsAsFactors=F)
str(ColnamesForSignalData)
colnames(ColnamesForSignalData)[1]<-"ColnamesForSignalData"
SampleInfo<-as.data.frame(SampleInfo, stringsAsFactors=F)
str(SampleInfo)
colnames(SampleInfo)[1]<-"ColnamesForSignalData"
SampleInfo_Reordered<-join(ColnamesForSignalData, SampleInfo, by="ColnamesForSignalData")
str(SampleInfo_Reordered)

#cleaning up the workspace:
rm(list = ls(pattern = "GSM")) 

#Taking a look at the Darmanis data:

#The last 3 rows of data are notes about the quality of the RNAseq, not actual measurements ("no feature"), ("ambiguous"), and ("alignment not unique")

#Some of this data isn't numeric. Why?
str(JoinedDarmanis[22086,])
#Ah - Alek has the first column being defined as the gene names (instead of using them as the column names)
#And what is variable "V2"? It looks like RNAseq data. 

is.numeric(JoinedDarmanis[22086,-c(1,2)])
#still having trouble, but no sign of non-numeric data in there so I think it is just unhappy because I'm treating a list as a matrix
is.numeric(as.matrix(JoinedDarmanis[22086,-c(1,2)]))
#That works.

GeneNamesForJoinedDarmanis<-JoinedDarmanis2[,1]
write.csv(GeneNamesForJoinedDarmanis, "GeneNamesForJoinedDarmanis.csv")

is.numeric(JoinedDarmanis)
JoinedDarmanis_AsNumMatrix<-JoinedDarmanis
#JoinedDarmanis_AsNumMatrix<-as.matrix(JoinedDarmanis)
#is.numeric(JoinedDarmanis_AsNumMatrix)

#Let's take a peek at their distribution:
png("Histogram_NumberOfReadsNoFeature.png")
hist(JoinedDarmanis_AsNumMatrix[22086,], col=4, breaks=40)
dev.off()
#Seems to be relatively normally distributed. There are some cell samples with a higher number of No feature reads than others (e.g., >2.5 E6), but I don't know enough about RNAseq to know if they are worth throwing out.

png("Histogram_NumberOfReadsAmbiguous.png")
hist(JoinedDarmanis_AsNumMatrix[22087,], col=2, breaks=40)
dev.off()
#Seems to be relatively normally distributed, but with a skew. There are some cell samples with a higher number of Ambiguous reads than others (e.g., >12500), but I don't know enough about RNAseq to know if they are worth throwing out.

png("Histogram_NumberOfReadsNotUniqueAlign.png")
hist(JoinedDarmanis_AsNumMatrix[22088,], col=3, breaks=40)
dev.off()
#Seems to be relatively normally distributed, but with a skew. There are some cell samples with a higher number of Non-Unique reads than others (e.g., >1,250,000), but I don't know enough about RNAseq to know if they are worth throwing out.

#Let's see if these types of measurements correlate:

png("Plot_NoFeatureVsAmbigAlignment.png")
plot(JoinedDarmanis_AsNumMatrix[22086,]~JoinedDarmanis_AsNumMatrix[22087,])
dev.off()

png("Plot_NoFeatureVsNotUniqueAlign.png")
plot(JoinedDarmanis_AsNumMatrix[22086,]~JoinedDarmanis_AsNumMatrix[22088,])
dev.off()

png("Plot_AmbigAlignmentVsNotUniqueAlign.png")
plot(JoinedDarmanis_AsNumMatrix[22087,]~JoinedDarmanis_AsNumMatrix[22088,])
dev.off()

#Yeah, so there may be some samples with generally lower quality of RNA, or perhaps a general range of quality. I really don't know enough about this to know if folks typically correct for that. Or perhaps certain cell types have more transcripts that aren't well mapped?

#I wonder if the principal components of variation in the data correlate with those variables?  

#I should probably log transform the reads first before messing around with PCA:
#Let's check out the average distribution of the signal first:
temp<-apply(JoinedDarmanis_AsNumMatrix[-c(22086:22088),], 1, mean)
png("Histogram_AverageReadsPerProbe.png")
hist(temp, breaks=1000)
dev.off()
temp<-apply(JoinedDarmanis_AsNumMatrix[-c(22086:22088),], 1, max)
png("Histogram_MaxReadsPerProbe.png")
hist(temp, breaks=1000)
dev.off()
#Wow, that is super skewed with a huge number of probes with an average, max, or median of 0. Yeah, let's log transform that data.

#Oh wait - log transformation in RNAseq data is awkward, because there are 0's, which convert into -Inf. It seems that folks typically log transform shifted data instead (data+1)
JoinedDarmanis_AsNumMatrix_Log2<-log2((JoinedDarmanis_AsNumMatrix[-c(22086:22088),]+1))
row.names(JoinedDarmanis_AsNumMatrix_Log2)<-GeneNamesForJoinedDarmanis[-c(22086:22088)]
write.csv(JoinedDarmanis_AsNumMatrix_Log2, "JoinedDarmanis_AsNumMatrix_Log2.csv")

temp<-apply(JoinedDarmanis_AsNumMatrix_Log2, 1, mean)
png("Histogram_AverageReadsPerProbeLog2.png")
hist(temp, breaks=1000)
dev.off()

temp<-apply(JoinedDarmanis_AsNumMatrix_Log2, 1, max)
png("Histogram_MaxReadsPerProbeLog2.png")
hist(temp, breaks=1000)
dev.off()

#That's still super skewed, but not as bad as before. I'm guessing that it is difficult to get reads >0 for many transcripts in single cell data. That means that any analysis based on sample-sample correlations is going to be largely driven by a few highly-expressed transcripts.
#Apparently there are other transformation methods out there that work better for RNAseq - VST is implemented by DESeq. VOOM applies a tranformation to the read counts that supposedly then makes the RNAseq data compatible with analyses included in the limma package. https://seqqc.wordpress.com/2015/02/16/should-you-transform-rna-seq-data-log-vst-voom/


#Time to recycle some code:
library(gdata)
library(fields)
library(stats)
library(car)
library(affy)
library(preprocessCore)
library(multtest)

#The fact that this data (unlike our microarray data) isn't quantile normalized may make the results of PCA a little funky. Let's start with sample-sample correlations:

#Visualize the sample-sample correlations using a heatmap:
png("09 Sample Sample Correlations Heatmap.png")
image(cor(JoinedDarmanis_AsNumMatrix_Log2), main="Visualizing the correlations between entire samples (by index#)", xlab="Red=Less correlated, Light yellow=Highly correlated")
dev.off()
#Note that the heatmap can be tailored to focus on a certain level of correlation by using the command zlim=c(lower limit, upper limit)

#Visualize the sample-sample correlations using a boxplot:
png("09 Boxplot Sample Sample Correlations.png", width=2000, height=300)
boxplot(data.frame(cor(JoinedDarmanis_AsNumMatrix_Log2)), cex=0.25, las=3, par(cex.axis=0.75), main="Boxplot of sample-sample correlations", xlab="Subject", ylab="Sample-Sample Correlations")
Median10thQuantile<-median(apply((cor(JoinedDarmanis_AsNumMatrix_Log2)), 1, quantile, 0.1))
MedianQuantile<-median(apply((cor(JoinedDarmanis_AsNumMatrix_Log2)), 1, quantile, 0.5))
abline(a=Median10thQuantile, b=0, col=2)
abline(a=MedianQuantile, b=0, col=3)
mtext(paste("Median Sample-Sample Correlation=", round(MedianQuantile, digits=3), sep=" ")) 
dev.off()

#Yes, there are whole groups of samples that very much do not correlate with each other. 

#It might be worthwile to perhaps try running things with all the genes with extremely low reads thrown out:

temp<-apply(JoinedDarmanis_AsNumMatrix_Log2, 1, max)
sum(temp==1)
#[1] 1620
#That's actually not horrible - so there are only a small number of genes that had 1 read or less to begin with.
sum(temp<3)
#[1] 3442
#And less than 4000 with 7 read or fewer to begin with maximum.  
sum(temp<3)/length(temp)
#[1] 0.1558524
#...or 15%.


#Run principal components analysis (PCA) to determine which major gradients of sample-sample correlations exist in the data (i.e. who is similar to whom):
pcaNormFiltered<-prcomp(t(JoinedDarmanis_AsNumMatrix_Log2))
tmp<-pcaNormFiltered$x[,1:4]
write.table(tmp, "PCA_1_4.txt", sep="\t")

PC1<-pcaNormFiltered$x[,1]
PC2<-pcaNormFiltered$x[,2]

PC3<-pcaNormFiltered$x[,3]
PC4<-pcaNormFiltered$x[,4]

#Output a scree plot for the PCA:
png("09 PCA Scree Plot1.png")
plot(summary(pcaNormFiltered)$importance[2,]~(c(1:ncol(JoinedDarmanis_AsNumMatrix_Log2))), main="Variance Explained by Each Principal Component", xlab="PC #", ylab="Proportion of Variance Explained", col=2)
dev.off()

png("09 PCA Scree Plot2.png")
plot(summary(pcaNormFiltered)$importance[3,]~(c(1:ncol(JoinedDarmanis_AsNumMatrix_Log2))), main="Variance Explained by Each Principal Component", xlab="PC #", ylab="Cumulative Proportion of Variance Explained", col=3)
dev.off()

#Output a scatterplot illustrating the relationship between Principal components 1 & 2 (PC1 & PC2):
png("09 PC1 vs PC2.png")
plot(PC1~PC2, main="Principal Components Analysis of Normalized Filtered Data", col=as.factor(SampleInfo_Reordered[,3]))
dev.off()

#Output a scatterplot illustrating the relationship between Principal components 1 & 2 (PC1 & PC2):
png("09 PC2 vs PC3.png")
plot(PC2~PC3, main="Principal Components Analysis of Normalized Filtered Data", col=as.factor(SampleInfo_Reordered[,3]))
dev.off()

#Output a scatterplot illustrating the relationship between Principal components 3 & 4 (PC3 & PC4):
png("09 PC3 vs PC4.png")
plot(PC3~PC4, main="Principal Components Analysis of Normalized Filtered Data", col=as.factor(SampleInfo_Reordered[,3]))
dev.off()

#Nice clustering - just like they mentioned in the paper.

#I wonder how the top PCs correlate with the read stats: ("no feature"), ("ambiguous"), and ("alignment not unique")

png("PC1vsNoFeature.png")
plot(PC1~JoinedDarmanis_AsNumMatrix[22086,])
dev.off()
cor(PC1,JoinedDarmanis_AsNumMatrix[22086,])
#[1] -0.5423068
#That's a pretty strong correlation

png("PC2vsNoFeature.png")
plot(PC2~JoinedDarmanis_AsNumMatrix[22086,])
dev.off()

png("PC1vsAmbigAlign.png")
plot(PC1~JoinedDarmanis_AsNumMatrix[22087,])
dev.off()
cor(PC1,JoinedDarmanis_AsNumMatrix[22087,])
#[1] -0.6105314
#That's a pretty strong correlation

png("PC2vsAmbigAlign.png")
plot(PC2~JoinedDarmanis_AsNumMatrix[22087,])
dev.off()

png("PC1vsNotUniqueAlign.png")
plot(PC1~JoinedDarmanis_AsNumMatrix[22088,])
dev.off()
cor(PC1,JoinedDarmanis_AsNumMatrix[22088,])
#[1] -0.3295242
#A much weaker correlation

png("PC2vsNotUniqueAlign.png")
plot(PC2~JoinedDarmanis_AsNumMatrix[22088,])
dev.off

temp<-scale(JoinedDarmanis_AsNumMatrix[22086,])+scale(JoinedDarmanis_AsNumMatrix[22087,])

png("PC1vsAVEscaledNoFeatureAmbigAlign.png")
plot(PC1~temp[,1])
dev.off()
cor(cbind(PC1,temp[,1]))
#cor=-0.6727594
summary.lm(lm(PC1~JoinedDarmanis_AsNumMatrix[22086,]+JoinedDarmanis_AsNumMatrix[22087,]))
# Coefficients:
#                                       Estimate Std. Error t value Pr(>|t|)    
# (Intercept)                          1.312e+02  7.435e+00  17.653  < 2e-16 ***
# JoinedDarmanis_AsNumMatrix[22086, ] -4.529e-05  5.373e-06  -8.429 4.57e-16 ***
# JoinedDarmanis_AsNumMatrix[22087, ] -9.859e-03  8.412e-04 -11.721  < 2e-16 ***
# ---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 67.78 on 458 degrees of freedom
# Multiple R-squared:  0.457,	Adjusted R-squared:  0.4546 
# F-statistic: 192.7 on 2 and 458 DF,  p-value: < 2.2e-16

#So even though they correlate, both are highly related to PC1 in their own way. Interesting. 
#later I'll check to see whether that is just artifact or reflects a particular cell type.

#Also, note that some of the colnames are still a little funky - I may need to go back and re-do Alek's compilation of the individual files.
#Done - I went back and re-did it.

