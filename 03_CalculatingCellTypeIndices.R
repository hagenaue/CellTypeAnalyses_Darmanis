#These code documents contain all of the code documenting the cell type analyses that we applied to Darmanis' single-cell RNAseq data to validate our methodology
#Megan Hagenauer and Alek Pankonin
#October 27, 2016

#**************************************


#cell type analysis: stolen from the Allen Brain Atlas code


#Reading in the database of cell type specific gene expression:

CellTypeSpecificGenes_Master3<-read.csv("CellTypeSpecificGenes_Master3.csv", header=T)

colnames(CellTypeSpecificGenes_Master3)

 [1] "Umbrella.Cell.Type"    "Specific.Cell.Type"    "Brain.Region"          "Gene.Symbol..Human."  
 [5] "Gene.Symbol..Mouse."   "Species"               "Age"                   "Statistical.Criterion"
 [9] "Specificity"           "Comparison"            "Platform"              "Citation"             
[13] "Tag"                   "CellType_Primary" 

table(CellTypeSpecificGenes_Master3$CellType_Primary)

#Note: I had to reverse the order of the next 3 lines to properly remove NAs - I should change that before attempting to release the code in a way that other people will use.

colnames(CellTypeSpecificGenes_Master3)[4]<-"GeneSymbol_Human"

sum(is.na(CellTypeSpecificGenes_Master3$GeneSymbol_Human))
#[1] 364

CellTypeSpecificGenes_Master3NoNA<-CellTypeSpecificGenes_Master3[is.na(CellTypeSpecificGenes_Master3$GeneSymbol_Human)==F,]

CellTypeSpecificGenes_Master3NoNA[,4]<-as.character(CellTypeSpecificGenes_Master3NoNA[,4])
str(CellTypeSpecificGenes_Master3NoNA)

#Note to self - if I try to make a generalizable version of this code, I will need to change it so that it references the rodent gene symbols instead of the human gene symbols.

########################

#Getting the Darmanis data ready for the cell type analysis:

length(unique(GeneNamesForJoinedDarmanis))
#All of the gene names are unique, so I think Alek's code for averaging by gene is unnecessary.

#Instead, I'm just going to z-score my log-transformed data.
#I wonder where the NAs are coming from.
#Ah -after messing around, I believe either NaNs or Inf are produced whenever a gene doesn't have any variability associated with it
#So let's remove the completely invariable data first:

JoinedDarmanis_StDev<-apply(JoinedDarmanis_AsNumMatrix_Log2, 1, sd) 
sum(JoinedDarmanis_StDev==0)
#[1] 458

JoinedDarmanis_AsNumMatrix_Log2_NoSD0<-JoinedDarmanis_AsNumMatrix_Log2[JoinedDarmanis_StDev>0,]
temp<-GeneNamesForJoinedDarmanis[-c(22086:22088)]
GeneNamesForJoinedDarmanis_NoSD0<-temp[JoinedDarmanis_StDev>0]


ZscoreDarmanis<-t(scale(t(JoinedDarmanis_AsNumMatrix_Log2_NoSD0), center=T, scale=T))#Zscores the data 
write.csv(ZscoreDarmanis, "ZscoreDarmanis.csv")
sum(is.na(ZscoreDarmanis))
#looks good now

#ZscoreDarmanis<-as.data.frame(JoinedDarmanis_AsNumMatrix_Log2)


#######################

#Pulling out the cell type specific gene expression from the Darmanis data:
temp<-data.frame(row.names(ZscoreDarmanis), ZscoreDarmanis, stringsAsFactors=F)

colnames(temp)[1]<-"GeneSymbol_Human"

sum(is.na(temp[,1]))
#[1] 0


sum(temp[,1] %in% CellTypeSpecificGenes_Master3[,4])
#[1] 2374

sum(CellTypeSpecificGenes_Master3[,4]  %in%  temp[,1])
# [1] 2882

#Note: NAs were causing a serious problem with this join function.  Fixed now. :)
library(plyr)
ZscoreDarmanis_Expression_CellType<-join(CellTypeSpecificGenes_Master3, temp, by="GeneSymbol_Human", type="inner")
dim(ZscoreDarmanis_Expression_CellType)
# [1] 2882  480
#It is making all possible combinations - some of the cell type specific genes are found in more than one index.

write.csv(ZscoreDarmanis_Expression_CellType, "ZscoreDarmanis_Expression_CellType.csv")

###############################################

AVE_Expression_CellType_Primary_bySample<-matrix(NA, nrow=length(names(table(ZscoreDarmanis_Expression_CellType$CellType_Primary))), ncol=(ncol(ZscoreDarmanis_Expression_CellType)-14))

row.names(AVE_Expression_CellType_Primary_bySample)<-names(table(ZscoreDarmanis_Expression_CellType$CellType_Primary))
colnames(AVE_Expression_CellType_Primary_bySample)<-colnames(temp)[-1]


for(i in c(15:ncol(ZscoreDarmanis_Expression_CellType))){
AVE_Expression_CellType_Primary_bySample[,(i-14)]<-tapply(ZscoreDarmanis_Expression_CellType[,i], ZscoreDarmanis_Expression_CellType$CellType_Primary, function(y) mean(y, na.rm=T))
}

head(AVE_Expression_CellType_Primary_bySample)

png("CorrMatrixCellTypeVsCellType_HeatMap.png")
heatmap(cor(t(AVE_Expression_CellType_Primary_bySample[,is.na(AVE_Expression_CellType_Primary_bySample[1,])==F])), margins = c(15, 15), cex.lab=0.5)
dev.off()

CorrelationMatrixCellTypeVsCellType<-cor(t(AVE_Expression_CellType_Primary_bySample[,is.na(AVE_Expression_CellType_Primary_bySample[1,])==F]))

write.csv(CorrelationMatrixCellTypeVsCellType, "CorrelationMatrixCellTypeVsCellType.csv")
#Huh - all correlations are positive. Perhaps because some samples simply have less reads or more artifacts?  

#let's see some examples:

png("AstrocyteByEndothelial.png")
plot(AVE_Expression_CellType_Primary_bySample[1,]~AVE_Expression_CellType_Primary_bySample[2,])
dev.off()
#Interesting - 3 distinct populations, and one outlier that looks like both. 

AVE_Expression_CellType_Tag_bySample<-matrix(NA, nrow=length(names(table(ZscoreDarmanis_Expression_CellType$Tag))), ncol=ncol(ZscoreDarmanis_Expression_CellType)-14)
row.names(AVE_Expression_CellType_Tag_bySample)<-names(table(ZscoreDarmanis_Expression_CellType$Tag))
colnames(AVE_Expression_CellType_Tag_bySample)<-colnames(temp)[-1]

for(i in c(15:ncol(ZscoreDarmanis_Expression_CellType))){
AVE_Expression_CellType_Tag_bySample[,(i-14)]<-tapply(ZscoreDarmanis_Expression_CellType[,i], ZscoreDarmanis_Expression_CellType$Tag, mean)
}

head(AVE_Expression_CellType_Tag_bySample)

png("CorrMatrixCellIndexVsCellIndex_HeatMap.png", width=1000, height=1000)
heatmap(cor(t(AVE_Expression_CellType_Tag_bySample[,is.na(AVE_Expression_CellType_Tag_bySample[1,])==F])), cex.lab=0.3, margins = c(20, 20))
dev.off()

CorrelationMatrixCellIndexVsCellIndex<-cor(t(AVE_Expression_CellType_Tag_bySample[,is.na(AVE_Expression_CellType_Tag_bySample[1,])==F]))

write.csv(CorrelationMatrixCellIndexVsCellIndex, "CorrelationMatrixCellIndexVsCellIndex.csv")

################################################

#Alright, so part of the trouble here is that there hasn't been any removal of overlapping probes yet, and we aren't averaging by tag. Let's go ahead and do that.


#Making a storage matrix to store information about overlap between primary indices:
CellTypeSpecificGenes_Master3_Overlap<-matrix(0, length(table(ZscoreDarmanis_Expression_CellType$CellType_Primary)), length(table(ZscoreDarmanis_Expression_CellType$CellType_Primary)) )

colnames(CellTypeSpecificGenes_Master3_Overlap)<-names(table(ZscoreDarmanis_Expression_CellType$CellType_Primary))
row.names(CellTypeSpecificGenes_Master3_Overlap)<-names(table(ZscoreDarmanis_Expression_CellType$CellType_Primary))

#Quantifying overlap between primary cell type indices:
for(i in 1: length(table(ZscoreDarmanis_Expression_CellType$CellType_Primary))){
  for(j in 1: length(table(ZscoreDarmanis_Expression_CellType$CellType_Primary))){
    
    CellTypeSpecificGenes_Master3_Overlap[i,j]<-sum(ZscoreDarmanis_Expression_CellType[ZscoreDarmanis_Expression_CellType$CellType_Primary==names(table(ZscoreDarmanis_Expression_CellType$CellType_Primary)[i]), 4]%in%ZscoreDarmanis_Expression_CellType[ZscoreDarmanis_Expression_CellType$CellType_Primary==names(table(ZscoreDarmanis_Expression_CellType$CellType_Primary)[j]), 4])/length(ZscoreDarmanis_Expression_CellType[ZscoreDarmanis_Expression_CellType$CellType_Primary==names(table(ZscoreDarmanis_Expression_CellType$CellType_Primary)[i]), 4])
    
  }
}

write.csv(CellTypeSpecificGenes_Master3_Overlap, "CellTypeSpecificGenes_Master3_Overlap.csv")



#What happens if we eliminate overlap between primary categories and then make master indices:

dim(ZscoreDarmanis_Expression_CellType)
# [1] 2882  480

#Making an empty first row for the storage matrix:
ZscoreDarmanis_Expression_CellType_NoPrimaryOverlap<-matrix(0, 1, (length(ZscoreDarmanis_Expression_CellType[1,])))
colnames(ZscoreDarmanis_Expression_CellType_NoPrimaryOverlap)<-colnames(ZscoreDarmanis_Expression_CellType)

for(i in 1: length(table(ZscoreDarmanis_Expression_CellType$CellType_Primary))){
  
  #Choosing all data for a particular primary cell type:
  TempCurrentIndexAllInfo<-ZscoreDarmanis_Expression_CellType[ZscoreDarmanis_Expression_CellType$CellType_Primary==names(table(ZscoreDarmanis_Expression_CellType$CellType_Primary)[i]), ] 
  
  #All of the gene symbols within the current primary cell type:
  TempCurrentIndex<-ZscoreDarmanis_Expression_CellType[ZscoreDarmanis_Expression_CellType$CellType_Primary==names(table(ZscoreDarmanis_Expression_CellType$CellType_Primary)[i]), 4] 
  
  #All of the gene symbols within all other primary cell types:
  TempAllOtherIndices<-ZscoreDarmanis_Expression_CellType[ZscoreDarmanis_Expression_CellType$CellType_Primary%in%names(table(ZscoreDarmanis_Expression_CellType$CellType_Primary)[-i]), 4]
  
  #Grabs only rows of data with gene symbols not found in other primary cell type indices:
  ZscoreDarmanis_Expression_CellType_NoPrimaryOverlap<-rbind(ZscoreDarmanis_Expression_CellType_NoPrimaryOverlap, TempCurrentIndexAllInfo[(TempCurrentIndex%in%TempAllOtherIndices)==F,])
  
}

dim(ZscoreDarmanis_Expression_CellType_NoPrimaryOverlap)
# [1] 2464  480

#removing that one dummy row:
ZscoreDarmanis_Expression_CellType_NoPrimaryOverlap<-ZscoreDarmanis_Expression_CellType_NoPrimaryOverlap[-1,]

dim(ZscoreDarmanis_Expression_CellType_NoPrimaryOverlap)
# [1] 2463  480

write.csv(ZscoreDarmanis_Expression_CellType_NoPrimaryOverlap, "ZscoreDarmanis_Expression_CellType_NoPrimaryOverlap.csv")


CellTypeSpecificGenes_Master3_PrimaryTable_NoPrimaryOverlap<-table(ZscoreDarmanis_Expression_CellType_NoPrimaryOverlap$CellType_Primary)

write.csv(CellTypeSpecificGenes_Master3_PrimaryTable_NoPrimaryOverlap, "CellTypeSpecificGenes_Master3_PrimaryTable_NoPrimaryOverlap.csv")

ZscoreDarmanis_Expression_CellType_NoPrimaryOverlap_Mean<-matrix(0, nrow=length(table(ZscoreDarmanis_Expression_CellType_NoPrimaryOverlap$CellType_Primary)), ncol=(length(ZscoreDarmanis_Expression_CellType_NoPrimaryOverlap[1,])-14))

row.names(ZscoreDarmanis_Expression_CellType_NoPrimaryOverlap_Mean)<-names(table(ZscoreDarmanis_Expression_CellType_NoPrimaryOverlap$CellType_Primary))

colnames(ZscoreDarmanis_Expression_CellType_NoPrimaryOverlap_Mean)<-colnames(ZscoreDarmanis_Expression_CellType_NoPrimaryOverlap)[-c(1:14)]

#Old version of code:
# for(i in c(15:length(ZscoreDarmanis_Expression_CellType_NoPrimaryOverlap[1,]))){
# ZscoreDarmanis_Expression_CellType_NoPrimaryOverlap_Mean[,i-14]<-tapply(ZscoreDarmanis_Expression_CellType_NoPrimaryOverlap[,i], ZscoreDarmanis_Expression_CellType_NoPrimaryOverlap$CellType_Primary, mean)
# }


#I went back and changed this so that it averaged by tag first, then by primary cell category

temp<-data.frame(ZscoreDarmanis_Expression_CellType_NoPrimaryOverlap$Tag, ZscoreDarmanis_Expression_CellType_NoPrimaryOverlap$CellType_Primary) 
dim(temp)
# [1] 2463    2

CellTypePrimaryVsTag<-unique(temp)
dim(CellTypePrimaryVsTag)
# [1] 38  2

colnames(CellTypePrimaryVsTag)<-c("Tag","CellType_Primary")
head(CellTypePrimaryVsTag)
#                                    Tag CellType_Primary
# 1     Astrocyte_All_Darmanis_PNAS_2015        Astrocyte
# 21     Astrocyte_All_Cahoy_JNeuro_2008        Astrocyte
# 72     Astrocyte_All_Zhang_JNeuro_2014        Astrocyte
# 102      Astrocyte_All_Doyle_Cell_2008        Astrocyte
# 118  Astrocyte_All_Zeisel_Science_2015        Astrocyte
# 309 Endothelial_All_Darmanis_PNAS_2015      Endothelial

ZscoreDarmanis_Expression_CellType_NoPrimaryOverlap_MeanTag<-matrix(0, nrow=length(table(ZscoreDarmanis_Expression_CellType_NoPrimaryOverlap$Tag)), ncol=(length(ZscoreDarmanis_Expression_CellType_NoPrimaryOverlap[1,])-14))

row.names(ZscoreDarmanis_Expression_CellType_NoPrimaryOverlap_MeanTag)<-names(table(ZscoreDarmanis_Expression_CellType_NoPrimaryOverlap$Tag))

colnames(ZscoreDarmanis_Expression_CellType_NoPrimaryOverlap_MeanTag)<-colnames(ZscoreDarmanis_Expression_CellType_NoPrimaryOverlap)[-c(1:14)]

for(i in c(15:length(ZscoreDarmanis_Expression_CellType_NoPrimaryOverlap[1,]))){
  ZscoreDarmanis_Expression_CellType_NoPrimaryOverlap_MeanTag[,i-14]<-tapply(ZscoreDarmanis_Expression_CellType_NoPrimaryOverlap[,i], ZscoreDarmanis_Expression_CellType_NoPrimaryOverlap$Tag, mean)
}

head(ZscoreDarmanis_Expression_CellType_NoPrimaryOverlap_MeanTag)

write.csv(ZscoreDarmanis_Expression_CellType_NoPrimaryOverlap_MeanTag, "ZscoreDarmanis_Expression_CellType_NoPrimaryOverlap_MeanTag.csv")


#Making histograms for each cell type tag:

for(i in 1:nrow(ZscoreDarmanis_Expression_CellType_NoPrimaryOverlap_MeanTag)){
  png(paste("Histogram_", row.names(ZscoreDarmanis_Expression_CellType_NoPrimaryOverlap_MeanTag)[i], ".png", sep=""))
  hist(ZscoreDarmanis_Expression_CellType_NoPrimaryOverlap_MeanTag[i,], main=row.names(ZscoreDarmanis_Expression_CellType_NoPrimaryOverlap_MeanTag)[i], breaks=40, col=i)
  dev.off()
}

#Interesting - it is easy to say what is *not* the cell type of interest, but the values for what could be the cell type of interest range greatly. I'm guessing that this is partially a property of the skewed variability and signal values in the data itself, but I'm not sure. I wonder if it correlates at all with the read qc stats for the samples.



png("Heatmap_CellType_NoPrimaryOverlap_MeanTag.png", height=1000, width=1000)
heatmap(cor(t(ZscoreDarmanis_Expression_CellType_NoPrimaryOverlap_MeanTag[,is.na(ZscoreDarmanis_Expression_CellType_NoPrimaryOverlap_MeanTag[1,])==F])), cex.lab=0.3, margins = c(20, 20))
dev.off()

CellType_NoPrimaryOverlap_MeanTag_CorrMatrix<-cor(t(ZscoreDarmanis_Expression_CellType_NoPrimaryOverlap_MeanTag[,is.na(ZscoreDarmanis_Expression_CellType_NoPrimaryOverlap_MeanTag[1,])==F]))

head(CellType_NoPrimaryOverlap_MeanTag_CorrMatrix)

write.csv(CellType_NoPrimaryOverlap_MeanTag_CorrMatrix, "CellType_NoPrimaryOverlap_MeanTag_CorrMatrix.csv")

#Interestingly, in this dataset, Ziesel's predictions are not as good as the other predictions - perhaps because Ziesel included a larger number of genes (less specific?)
#Mural cells don't seem to be clearly distinguished from other cell types.
#Doyle's predictions still don't match anyone else for the most part, even for cells like oligodendrocytes.
#OPCs still don't seem specific.
#And predictions made using the Darmanis dataset are the most clear. Go Figure.


temp2<-tapply(ZscoreDarmanis_Expression_CellType_NoPrimaryOverlap[,i], ZscoreDarmanis_Expression_CellType_NoPrimaryOverlap$Tag, mean)
Tag<-names(temp2)

CellTypePrimaryVsTag2<-join(as.data.frame(Tag), as.data.frame(CellTypePrimaryVsTag), by="Tag")


for(i in c(1:length(ZscoreDarmanis_Expression_CellType_NoPrimaryOverlap_MeanTag[1,]))){
  ZscoreDarmanis_Expression_CellType_NoPrimaryOverlap_Mean[,i]<-tapply(ZscoreDarmanis_Expression_CellType_NoPrimaryOverlap_MeanTag[,i], CellTypePrimaryVsTag2[,2], mean)
}

head(ZscoreDarmanis_Expression_CellType_NoPrimaryOverlap_Mean)


write.csv(ZscoreDarmanis_Expression_CellType_NoPrimaryOverlap_Mean, "ZscoreDarmanis_Expression_CellType_NoPrimaryOverlap_Mean.csv")

#Making histograms for each primary cell type:

for(i in 1:nrow(ZscoreDarmanis_Expression_CellType_NoPrimaryOverlap_Mean)){
  png(paste("Histogram_", row.names(ZscoreDarmanis_Expression_CellType_NoPrimaryOverlap_Mean)[i], ".png", sep=""))
  hist(ZscoreDarmanis_Expression_CellType_NoPrimaryOverlap_Mean[i,], main=row.names(ZscoreDarmanis_Expression_CellType_NoPrimaryOverlap_Mean)[i], breaks=40, col=i)
  dev.off()
}

#Interesting - it is easy to say what is *not* the cell type of interest, but the values for what could be the cell type of interest range greatly. I'm guessing that this is partially a property of the skewed variability and signal values in the data itself, but I'm not sure. I wonder if it correlates at all with the read qc stats for the samples.


is.numeric(ZscoreDarmanis_Expression_CellType_NoPrimaryOverlap_Mean)

png("Heatmap_CorMatrixPrimaryCellsNoOverlap.png", height=1000, width=1000)
heatmap(cor(t(ZscoreDarmanis_Expression_CellType_NoPrimaryOverlap_Mean[,is.na(ZscoreDarmanis_Expression_CellType_NoPrimaryOverlap_Mean[1,])==F])), cex.lab=0.3, margins = c(20, 20))
dev.off()

CellTypeIndexNoNA3_NormBest_NoPrimaryOverlap_CorrMatrix<-cor(t(ZscoreDarmanis_Expression_CellType_NoPrimaryOverlap_Mean[,is.na(ZscoreDarmanis_Expression_CellType_NoPrimaryOverlap_Mean[1,])==F]))

head(CellTypeIndexNoNA3_NormBest_NoPrimaryOverlap_CorrMatrix)

write.csv(CellTypeIndexNoNA3_NormBest_NoPrimaryOverlap_CorrMatrix, "CellTypeIndexNoNA3_NormBest_NoPrimaryOverlap_CorrMatrix.csv")

##################################################

