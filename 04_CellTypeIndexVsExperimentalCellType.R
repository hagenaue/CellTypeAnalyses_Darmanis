#These code documents contain all of the code documenting the cell type analyses that we applied to Darmanis' single-cell RNAseq data to validate our methodology
#Megan Hagenauer and Alek Pankonin
#October 27, 2016

#**************************************


#Comparing our cell type indices to the cell type identity provided in the paper (GEO website)
#Note that they determined cell identity by data clustering (not morphology)


ExperimentalCellType<-SampleInfo_Reordered[,3]

#Plotting all the cell type tags for each cell type:
ZscoreDarmanis_Expression_CellType_NoPrimaryOverlap_MeanTag_byActualCellType<-matrix(0, nrow=nrow(ZscoreDarmanis_Expression_CellType_NoPrimaryOverlap_MeanTag), ncol=length(unique(ExperimentalCellType)))
row.names(ZscoreDarmanis_Expression_CellType_NoPrimaryOverlap_MeanTag_byActualCellType)<-row.names(ZscoreDarmanis_Expression_CellType_NoPrimaryOverlap_MeanTag)
colnames(ZscoreDarmanis_Expression_CellType_NoPrimaryOverlap_MeanTag_byActualCellType)<-names(table(ExperimentalCellType))

for(i in 1: nrow(ZscoreDarmanis_Expression_CellType_NoPrimaryOverlap_MeanTag)){
  ZscoreDarmanis_Expression_CellType_NoPrimaryOverlap_MeanTag_byActualCellType[i,]<-tapply(ZscoreDarmanis_Expression_CellType_NoPrimaryOverlap_MeanTag[i,], ExperimentalCellType, mean)
}

write.csv(ZscoreDarmanis_Expression_CellType_NoPrimaryOverlap_MeanTag_byActualCellType, "ZscoreDarmanis_Expression_CellType_NoPrimaryOverlap_MeanTag_byActualCellType.csv")


ColorsForCellTypes<-as.character(CellTypePrimaryVsTag[,2])
ColorsForCellTypes[ColorsForCellTypes=="Astrocyte"]<-"lavender"
ColorsForCellTypes[ColorsForCellTypes=="Endothelial"]<-"orange"
ColorsForCellTypes[ColorsForCellTypes=="Microglia"]<-"darkolivegreen3"
ColorsForCellTypes[ColorsForCellTypes=="Mural"]<-"yellow"
ColorsForCellTypes[ColorsForCellTypes=="Neuron_All"]<-"darkviolet"
ColorsForCellTypes[ColorsForCellTypes=="Neuron_Interneuron"]<-"blue"
ColorsForCellTypes[ColorsForCellTypes=="Neuron_Projection"]<-"red"
ColorsForCellTypes[ColorsForCellTypes=="Oligodendrocyte"]<-"lightpink"
ColorsForCellTypes[ColorsForCellTypes=="Oligodendrocyte_Immature"]<-"ivory2"
ColorsForCellTypes[ColorsForCellTypes=="RBC"]<-"indianred4"

#There is a mysterious legend that keeps popping up... that is just the *original* margin parameters. Um. Yeah. Weird, right?

for(i in 1:length(unique(ExperimentalCellType))){
  png(paste("CellTypeTagVs_", colnames(ZscoreDarmanis_Expression_CellType_NoPrimaryOverlap_MeanTag_byActualCellType)[i], ".png", sep=""), height=800, width=800)
  barplot(height=ZscoreDarmanis_Expression_CellType_NoPrimaryOverlap_MeanTag_byActualCellType[,i], width=0.3, space = NULL, main=colnames(ZscoreDarmanis_Expression_CellType_NoPrimaryOverlap_MeanTag_byActualCellType)[i], ylab="Average Cell Type Index (Mean Z-score for All Genes)", xlab=NULL, names.arg=row.names(ZscoreDarmanis_Expression_CellType_NoPrimaryOverlap_MeanTag_byActualCellType), las=2, par(mar=c(25,6,4,2)+0.1), cex.lab=1.3, cex.main=1.3, font.axis=2, font.lab=2, font.main=2,col=ColorsForCellTypes, args.legend = c(x="topright", legend=""))
  dev.off()
}


#Note that any references to the cell type specific gene lists originally derived from Darmanis are self-referential.


str(ZscoreDarmanis_Expression_CellType_NoPrimaryOverlap_MeanTag)
table(ExperimentalCellType)

CorrelationCoefficients_CellTypeIndexVsCellType<-matrix(0, 38, 6)
colnames(CorrelationCoefficients_CellTypeIndexVsCellType)<-names(table(ExperimentalCellType[which(ExperimentalCellType%in%c("fetal_quiescent","fetal_replicating","hybrid")==F)]))
row.names(CorrelationCoefficients_CellTypeIndexVsCellType)<-row.names(ZscoreDarmanis_Expression_CellType_NoPrimaryOverlap_MeanTag)

ZscoreDarmanis_Expression_CellType_NoPrimaryOverlap_MeanTag_noFetal<-ZscoreDarmanis_Expression_CellType_NoPrimaryOverlap_MeanTag[,which(ExperimentalCellType%in%c("fetal_quiescent","fetal_replicating","hybrid")==F)]
str(ZscoreDarmanis_Expression_CellType_NoPrimaryOverlap_MeanTag_noFetal)

ExperimentalCellType_NoFetal<-ExperimentalCellType[which(ExperimentalCellType%in%c("fetal_quiescent","fetal_replicating","hybrid")==F)]
str(ExperimentalCellType_NoFetal)

for(i in c(1:38)){
  for(j in c(1:6)){
    temp<-summary.lm(lm(ZscoreDarmanis_Expression_CellType_NoPrimaryOverlap_MeanTag_noFetal[i,]~(ExperimentalCellType_NoFetal==names(table(ExperimentalCellType_NoFetal))[j])))
    CorrelationCoefficients_CellTypeIndexVsCellType[i,j]<-temp$r.squared
  }
}

write.csv(CorrelationCoefficients_CellTypeIndexVsCellType, "CorrelationCoefficients_CellTypeIndexVsCellType.csv")
