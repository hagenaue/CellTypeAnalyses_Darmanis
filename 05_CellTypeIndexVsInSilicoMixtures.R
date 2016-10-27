#These code documents contain all of the code documenting the cell type analyses that we applied to Darmanis' single-cell RNAseq data to validate our methodology
#Megan Hagenauer and Alek Pankonin
#October 27, 2016

#**************************************

#determining how tightly we can track in silico ratios:

AstrocyteMixtures<-matrix(0, length(names(table(ZscoreDarmanis_Expression_CellType_NoPrimaryOverlap$Tag))), 30)
row.names(AstrocyteMixtures)<-names(table(ZscoreDarmanis_Expression_CellType_NoPrimaryOverlap$Tag))
AstrocyteMixtureRatios<-c(rep(15, 5), rep(20, 5), rep(25, 5), rep(30,5), rep(35,5), rep(40,5))

for(i in c(15, 20, 25, 30, 35, 40)){
  
  temp3<-matrix(0, length(names(table(ZscoreDarmanis_Expression_CellType_NoPrimaryOverlap$Tag))), 5)  
  
  for(j in 1:5){
    temp<-cbind(ZscoreDarmanis_Expression_CellType_NoPrimaryOverlap[, sample(which(ExperimentalCellType=="astrocytes"), i, replace=T)+14],ZscoreDarmanis_Expression_CellType_NoPrimaryOverlap[, sample(which(ExperimentalCellType!="astrocytes"&ExperimentalCellType!="fetal_quiescent"&ExperimentalCellType!="fetal_replicating"), (100-i), replace=T)+14])
    temp2<-apply(temp, 1, mean)
    temp3[,j]<-tapply(temp2, ZscoreDarmanis_Expression_CellType_NoPrimaryOverlap$Tag, mean)
  }
  
  AstrocyteMixtures[,c(i:(i+4))-14]<-temp3
}

head(AstrocyteMixtures)

png("AstrocyteMixtureRatiosVsIndex_NoDarmanis.png")
plot(apply(AstrocyteMixtures[c(1, 3:5),], 2, mean)~AstrocyteMixtureRatios, pch=18, ylab="AstrocyteIndex")
BestFitLine<-lm(apply(AstrocyteMixtures[c(1, 3:5),], 2, mean)~AstrocyteMixtureRatios)
abline(BestFitLine)
abline(a=(BestFitLine$coefficients[1]+summary.lm(BestFitLine)[[6]]), b=BestFitLine$coefficients[2], col="grey")
abline(a=(BestFitLine$coefficients[1]-summary.lm(BestFitLine)[[6]]), b=BestFitLine$coefficients[2], col="grey")

abline(a=BestFitLine$coefficients[1]+(BestFitLine$coefficients[2]*30), b=0, col="red")
abline(a=BestFitLine$coefficients[1]+(BestFitLine$coefficients[2]*27.5), b=0, col="green")
abline(a=BestFitLine$coefficients[1]+(BestFitLine$coefficients[2]*25), b=0, col="red")
dev.off()


OligodendrocyteMixtures<-matrix(0, length(names(table(ZscoreDarmanis_Expression_CellType_NoPrimaryOverlap$Tag))), 30)
row.names(OligodendrocyteMixtures)<-names(table(ZscoreDarmanis_Expression_CellType_NoPrimaryOverlap$Tag))
OligodendrocyteMixtureRatios<-c(rep(15, 5), rep(20, 5), rep(25, 5), rep(30,5), rep(35,5), rep(40,5))

for(i in c(15, 20, 25, 30, 35, 40)){
  
  temp3<-matrix(0, length(names(table(ZscoreDarmanis_Expression_CellType_NoPrimaryOverlap$Tag))), 5)  
  
  for(j in 1:5){
    temp<-cbind(ZscoreDarmanis_Expression_CellType_NoPrimaryOverlap[, sample(which(ExperimentalCellType=="oligodendrocytes"), i, replace=T)+14],ZscoreDarmanis_Expression_CellType_NoPrimaryOverlap[, sample(which(ExperimentalCellType!="oligodendrocytes"&ExperimentalCellType!="fetal_quiescent"&ExperimentalCellType!="fetal_replicating"), (100-i), replace=T)+14])
    temp2<-apply(temp, 1, mean)
    temp3[,j]<-tapply(temp2, ZscoreDarmanis_Expression_CellType_NoPrimaryOverlap$Tag, mean)
  }
  
  OligodendrocyteMixtures[,c(i:(i+4))-14]<-temp3
}

head(OligodendrocyteMixtures)

png("OligodendrocyteMixtureRatiosVsIndex_NoDarmanis.png")
plot(apply(OligodendrocyteMixtures[c(29:31, 33:34),], 2, mean)~OligodendrocyteMixtureRatios, pch=18, ylab="OligodendrocyteIndex")
BestFitLine<-lm(apply(OligodendrocyteMixtures[c(29:31, 33:34),], 2, mean)~OligodendrocyteMixtureRatios)
abline(BestFitLine)
abline(a=(BestFitLine$coefficients[1]+summary.lm(BestFitLine)[[6]]), b=BestFitLine$coefficients[2], col="grey")
abline(a=(BestFitLine$coefficients[1]-summary.lm(BestFitLine)[[6]]), b=BestFitLine$coefficients[2], col="grey")

abline(a=BestFitLine$coefficients[1]+(BestFitLine$coefficients[2]*30), b=0, col="red")
abline(a=BestFitLine$coefficients[1]+(BestFitLine$coefficients[2]*27.5), b=0, col="green")
abline(a=BestFitLine$coefficients[1]+(BestFitLine$coefficients[2]*25), b=0, col="red")
dev.off()

MicrogliaMixtures<-matrix(0, length(names(table(ZscoreDarmanis_Expression_CellType_NoPrimaryOverlap$Tag))), 30)
row.names(MicrogliaMixtures)<-names(table(ZscoreDarmanis_Expression_CellType_NoPrimaryOverlap$Tag))
MicrogliaMixtureRatios<-c(rep(5, 5), rep(10, 5), rep(15, 5), rep(20,5), rep(25,5), rep(30,5))

for(i in c(5, 10, 15, 20, 25, 30)){
  
  temp3<-matrix(0, length(names(table(ZscoreDarmanis_Expression_CellType_NoPrimaryOverlap$Tag))), 5)  
  
  for(j in 1:5){
    temp<-cbind(ZscoreDarmanis_Expression_CellType_NoPrimaryOverlap[, sample(which(ExperimentalCellType=="microglia"), i, replace=T)+14],ZscoreDarmanis_Expression_CellType_NoPrimaryOverlap[, sample(which(ExperimentalCellType!="microglia"&ExperimentalCellType!="fetal_quiescent"&ExperimentalCellType!="fetal_replicating"), (100-i), replace=T)+14])
    temp2<-apply(temp, 1, mean)
    temp3[,j]<-tapply(temp2, ZscoreDarmanis_Expression_CellType_NoPrimaryOverlap$Tag, mean)
  }
  
  MicrogliaMixtures[,c(i:(i+4))-4]<-temp3
}

head(MicrogliaMixtures)



png("MicrogliaMixtureRatiosVsIndex_NoDarmanis.png")
plot(apply(MicrogliaMixtures[c(11:12),], 2, mean)~MicrogliaMixtureRatios, pch=18, ylab="MicrogliaIndex")
BestFitLine<-lm(apply(MicrogliaMixtures[c(11:12),], 2, mean)~MicrogliaMixtureRatios)
abline(BestFitLine)
abline(a=(BestFitLine$coefficients[1]+summary.lm(BestFitLine)[[6]]), b=BestFitLine$coefficients[2], col="grey")
abline(a=(BestFitLine$coefficients[1]-summary.lm(BestFitLine)[[6]]), b=BestFitLine$coefficients[2], col="grey")

abline(a=BestFitLine$coefficients[1]+(BestFitLine$coefficients[2]*15), b=0, col="red")
abline(a=BestFitLine$coefficients[1]+(BestFitLine$coefficients[2]*17.5), b=0, col="green")
abline(a=BestFitLine$coefficients[1]+(BestFitLine$coefficients[2]*20), b=0, col="red")
dev.off()

EndothelialMixtures<-matrix(0, length(names(table(ZscoreDarmanis_Expression_CellType_NoPrimaryOverlap$Tag))), 30)
row.names(EndothelialMixtures)<-names(table(ZscoreDarmanis_Expression_CellType_NoPrimaryOverlap$Tag))
EndothelialMixtureRatios<-c(rep(1, 5), rep(6, 5), rep(11, 5), rep(16,5), rep(21,5), rep(26,5))

for(i in c(1, 6, 11, 16, 21, 26)){
  
  temp3<-matrix(0, length(names(table(ZscoreDarmanis_Expression_CellType_NoPrimaryOverlap$Tag))), 5)  
  
  for(j in 1:5){
    temp<-cbind(ZscoreDarmanis_Expression_CellType_NoPrimaryOverlap[, sample(which(ExperimentalCellType=="endothelial"), i, replace=T)+14],ZscoreDarmanis_Expression_CellType_NoPrimaryOverlap[, sample(which(ExperimentalCellType!="endothelial"&ExperimentalCellType!="fetal_quiescent"&ExperimentalCellType!="fetal_replicating"), (100-i), replace=T)+14])
    temp2<-apply(temp, 1, mean)
    temp3[,j]<-tapply(temp2, ZscoreDarmanis_Expression_CellType_NoPrimaryOverlap$Tag, mean)
  }
  
  EndothelialMixtures[,c(i:(i+4))]<-temp3
}

head(EndothelialMixtures)


png("EndothelialMixtureRatiosVsIndex_NoDarmanis.png")
plot(apply(EndothelialMixtures[c(6,8:9),], 2, mean)~EndothelialMixtureRatios, pch=18, ylab="EndothelialIndex")
BestFitLine<-lm(apply(EndothelialMixtures[c(6,8:9),], 2, mean)~EndothelialMixtureRatios)
abline(BestFitLine)
abline(a=(BestFitLine$coefficients[1]+summary.lm(BestFitLine)[[6]]), b=BestFitLine$coefficients[2], col="grey")
abline(a=(BestFitLine$coefficients[1]-summary.lm(BestFitLine)[[6]]), b=BestFitLine$coefficients[2], col="grey")

abline(a=BestFitLine$coefficients[1]+(BestFitLine$coefficients[2]*5), b=0, col="red")
abline(a=BestFitLine$coefficients[1]+(BestFitLine$coefficients[2]*7.5), b=0, col="green")
abline(a=BestFitLine$coefficients[1]+(BestFitLine$coefficients[2]*10), b=0, col="red")
dev.off()

NeuronMixtures<-matrix(0, length(names(table(ZscoreDarmanis_Expression_CellType_NoPrimaryOverlap$Tag))), 30)
row.names(NeuronMixtures)<-names(table(ZscoreDarmanis_Expression_CellType_NoPrimaryOverlap$Tag))
NeuronMixtureRatios<-c(rep(25, 5), rep(30, 5), rep(35, 5), rep(40,5), rep(45,5), rep(50,5))

for(i in c(25, 30, 35, 40, 45, 50)){
  
  temp3<-matrix(0, length(names(table(ZscoreDarmanis_Expression_CellType_NoPrimaryOverlap$Tag))), 5)  
  
  for(j in 1:5){
    temp<-cbind(ZscoreDarmanis_Expression_CellType_NoPrimaryOverlap[, sample(which(ExperimentalCellType=="neurons"), i, replace=T)+14],ZscoreDarmanis_Expression_CellType_NoPrimaryOverlap[, sample(which(ExperimentalCellType!="neurons"&ExperimentalCellType!="fetal_quiescent"&ExperimentalCellType!="fetal_replicating"), (100-i), replace=T)+14])
    temp2<-apply(temp, 1, mean)
    temp3[,j]<-tapply(temp2, ZscoreDarmanis_Expression_CellType_NoPrimaryOverlap$Tag, mean)
  }
  
  NeuronMixtures[,c(i:(i+4))-24]<-temp3
}

head(NeuronMixtures)

png("NeuronMixtureRatiosVsIndex_NoDarmanis.png")
plot(apply(NeuronMixtures[c(16,18),], 2, mean)~NeuronMixtureRatios, pch=18, ylab="NeuronIndex")
BestFitLine<-lm(apply(NeuronMixtures[c(16,18),], 2, mean)~NeuronMixtureRatios)
abline(BestFitLine)
abline(a=(BestFitLine$coefficients[1]+summary.lm(BestFitLine)[[6]]), b=BestFitLine$coefficients[2], col="grey")
abline(a=(BestFitLine$coefficients[1]-summary.lm(BestFitLine)[[6]]), b=BestFitLine$coefficients[2], col="grey")

abline(a=BestFitLine$coefficients[1]+(BestFitLine$coefficients[2]*35), b=0, col="red")
abline(a=BestFitLine$coefficients[1]+(BestFitLine$coefficients[2]*37.5), b=0, col="green")
abline(a=BestFitLine$coefficients[1]+(BestFitLine$coefficients[2]*40), b=0, col="red")
dev.off()

