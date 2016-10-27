#These code documents contain all of the code documenting the cell type analyses that we applied to Darmanis' single-cell RNAseq data to validate our methodology
#Megan Hagenauer and Alek Pankonin
#October 27, 2016

#**************************************

##########Comparing results to PSEA:####################

source("https://bioconductor.org/biocLite.R")
biocLite("PSEA")

library("PSEA")
marker

#here is the code for the marker function that averages expression over the cell type specific marker probesets:
#expr is the expression dataset
#id is the list of cell type specific marker probesets

function (expr, id, sampleSubset = NULL, targetMean = 1) 
{
  if (class(expr) == "ExpressionSet") 
    expr <- exprs(expr)
  if (!is.matrix(expr)) {
    warning("expr is not a matrix. Coercing to matrix.")
    expr <- as.matrix(expr)
  }
  if (any(sapply(id, is.null))) 
    stop("An element of id is null.")
  if (is.null(sampleSubset)) 
    sampleSubset <- 1:ncol(expr)
  if (any(is.na(sampleSubset))) 
    stop("NA in sampleSubset.")
  mrkrs <- sapply(id, function(x) {
    
    #This code is just pulling out the expression data for the cell type specific markers:
    xind <- match(x, rownames(expr))
    if (any(is.na(xind))) 
      stop("Unmatched marker.")
    
    #This code is normalizing data: "The reference signal is scaled to have an average value determined by targetMean (1 by default)". "The expression value of each transcript is normalized before averaging over multiple transcripts so that they have same weight in the final reference signal (i.e to avoid that highly expressed transcripts dominate the reference signal)." Note - I believe they don't automatically correct for the number of transcripts representing a single gene. "If the expression matrix contains multiple measures of the same transcript (for instance several probes of a microarray mea- suring the same transcript) they can be averaged before taking the average over different transcripts."
    
    else apply(targetMean * expr[xind, , drop = FALSE]#The targetMean is multiplied by the expression for the markers
               /apply(expr[xind, sampleSubset, drop = FALSE], 1, mean, na.rm = TRUE), #Divided by the mean for each row for the expression of the markers
               2, mean, na.rm = TRUE)#So basically for each column, the expression values for each row are divided by the mean for the row, and then averaged. That is important to note because that means all of the values are essentially a fold change.
  })
  refSignal <- apply(t(mrkrs), 2, mean)
  return(refSignal)
}
<environment: namespace:PSEA>

  #Let's test out the normalization aspect of this:
  expr<-rbind(c(1.5, 2, 0, -1, 2.5), c(2.5, 3, 1, 0, 3.5), c(3.5, 4, 2, 1, 4.5), c(1, 1, 1, 1, 1))

#       [,1] [,2] [,3] [,4] [,5]
# [1,]  1.5    2    0   -1  2.5
# [2,]  2.5    3    1    0  3.5
# [3,]  3.5    4    2    1  2.5
# [4,]  1.0    1    1    1  1.0

  xind<-c(1:4)
  sampleSubset<-c(1:5)
  targetMean<-1
  
  temp<-apply(expr[xind, sampleSubset, drop = FALSE], 1, mean, na.rm = TRUE)
  #[1] 1 2 3 1  #Yep - that calculated the mean for each row.
  
  targetMean * expr[xind, , drop = FALSE]/temp
  
  #           [,1]     [,2]      [,3]       [,4] [,5]
  # [1,] 1.500000 2.000000 0.0000000 -1.0000000 2.50
  # [2,] 1.250000 1.500000 0.5000000  0.0000000 1.75
  # [3,] 1.166667 1.333333 0.6666667  0.3333333 1.50
  # [4,] 1.000000 1.000000 1.0000000  1.0000000 1.00
  
  #Does that fit what we would expect?
  #The top row should be entirely divided by 1 - yes
  #The second row should be halved - yes
  #The third row should be divided by 3 - yes
  #The fourth rown is still all ones - yes
  
  #So yeah, they are normalizing their data as a ratio of the original mean.  This could make results where the average value for a gene is close to 0 super funky - like in Darmanis' RNAseq data.
  
  apply(targetMean * expr[xind, , drop = FALSE]/temp, 2, mean, na.rm = TRUE)
  #[1] 1.22916667 1.45833333 0.54166667 0.08333333 1.68750000
  
  #And this is just averaging: 
  mean(c(1.5, 1.25, 1.1666, 1))
  #[1] 1.22915

  #Cool - I think I understand their method now. 
#So maybe it is worth running a head to head comparison - especially with the Darmanis data.
 
  #First I need to make an id list:
 # id List of strings. Names of the transcripts to use to generate the reference sig- nal. Names correspond to row names in expr if it is a matrix or row names in exprs(expr) if it is an ExpressionSet). 
  #To make things fair, I'm going to grab from the cell type specific gene list that has already been filtered to exclude overlap.
  #Note: unlike ours, their function doesn't include normalizing and averaging by gene first, and although their paper states that they filter the cell type specific markers based on co-expression they do not actually include that filtering step in their package, so really the only thing being explicitly compared here is the scaling step.
  
  
  
 colnames(ZscoreDarmanis_Expression_CellType_NoPrimaryOverlap)
  
 #GeneSymbol Human: ZscoreDarmanis_Expression_CellType_NoPrimaryOverlap[,4]
 #PrimaryCell Type: ZscoreDarmanis_Expression_CellType_NoPrimaryOverlap[,14]
 
 #I'm going to start without filtering for SD0 genes, then try a different version with filtering:
 #JoinedDarmanis_AsNumMatrix_Log2_NoSD0<-JoinedDarmanis_AsNumMatrix_Log2[JoinedDarmanis_StDev>0,]
 #temp<-GeneNamesForJoinedDarmanis[-c(22086:22088)]
 #GeneNamesForJoinedDarmanis_NoSD0<-temp[JoinedDarmanis_StDev>0]
 
 row.names(JoinedDarmanis_AsNumMatrix_Log2)<-GeneNamesForJoinedDarmanis[-c(22086:22088)]
 row.names(JoinedDarmanis_AsNumMatrix_Log2_NoSD0)<-GeneNamesForJoinedDarmanis_NoSD0
 
   
#For astrocytes:
 astrocyte_genes<-ZscoreDarmanis_Expression_CellType_NoPrimaryOverlap[ZscoreDarmanis_Expression_CellType_NoPrimaryOverlap[,14]=="Astrocyte",4]
 
astrocyte_reference<-marker(expr=JoinedDarmanis_AsNumMatrix_Log2, id=astrocyte_genes)

png("PSEA_AstrocyteVsCellType.png", width=1000, height=400) 
boxplot(astrocyte_reference~ExperimentalCellType, ylab="Astrocyte Reference Signal", col="lavender")
dev.off()
 
astrocyte_reference<-marker(expr=JoinedDarmanis_AsNumMatrix_Log2_NoSD0, id=astrocyte_genes)

png("PSEA_AstrocyteVsCellType_noSD0.png", width=1000, height=400) 
boxplot(astrocyte_reference~ExperimentalCellType, ylab="Astrocyte Reference Signal", col="lavender")
dev.off()

summary.lm(lm(astrocyte_reference~ExperimentalCellType))
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept)                            2.16385    0.08788  24.622  < 2e-16 ***
#   ExperimentalCellTypeendothelial       -1.18268    0.17795  -6.646 8.58e-11 ***
#   ExperimentalCellTypefetal_quiescent   -1.74148    0.10989 -15.847  < 2e-16 ***
#   ExperimentalCellTypefetal_replicating -1.38754    0.16394  -8.464 3.55e-16 ***
#   ExperimentalCellTypehybrid            -0.77211    0.13466  -5.734 1.79e-08 ***
#   ExperimentalCellTypemicroglia         -1.84364    0.19404  -9.501  < 2e-16 ***
#   ExperimentalCellTypeneurons           -1.07571    0.10680 -10.072  < 2e-16 ***
#   ExperimentalCellTypeoligodendrocytes  -1.55634    0.14143 -11.005  < 2e-16 ***
#   ExperimentalCellTypeOPC               -1.49400    0.18527  -8.064 6.54e-15 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 0.692 on 457 degrees of freedom
# Multiple R-squared:  0.3992,	Adjusted R-squared:  0.3886 
# F-statistic: 37.95 on 8 and 457 DF,  p-value: < 2.2e-16

row.names(ZscoreDarmanis_Expression_CellType_NoPrimaryOverlap_Mean)

png("Us_AstrocyteVsCellType_noSD0.png", width=1000, height=400) 
boxplot(ZscoreDarmanis_Expression_CellType_NoPrimaryOverlap_Mean[1,]~ExperimentalCellType, ylab="Astrocyte Index", col="lavender")
dev.off()
summary.lm(lm(ZscoreDarmanis_Expression_CellType_NoPrimaryOverlap_Mean[1,]~ExperimentalCellType))

# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept)                            0.84067    0.03145   26.73   <2e-16 ***
#   ExperimentalCellTypeendothelial       -1.01138    0.06368  -15.88   <2e-16 ***
#   ExperimentalCellTypefetal_quiescent   -1.12953    0.03932  -28.72   <2e-16 ***
#   ExperimentalCellTypefetal_replicating -0.94691    0.05867  -16.14   <2e-16 ***
#   ExperimentalCellTypehybrid            -0.64423    0.04819  -13.37   <2e-16 ***
#   ExperimentalCellTypemicroglia         -1.12173    0.06943  -16.16   <2e-16 ***
#   ExperimentalCellTypeneurons           -0.90227    0.03822  -23.61   <2e-16 ***
#   ExperimentalCellTypeoligodendrocytes  -1.04467    0.05061  -20.64   <2e-16 ***
#   ExperimentalCellTypeOPC               -0.99906    0.06630  -15.07   <2e-16 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 0.2476 on 457 degrees of freedom
# Multiple R-squared:  0.6797,	Adjusted R-squared:  0.6741 
# F-statistic: 121.2 on 8 and 457 DF,  p-value: < 2.2e-16

#So our method basically doubles the R-squared. To properly compare them, we should probably remove Darmanis though, and it might be good to know how much of this is just because we averaged by tag, diminishing the influence of Zeisel.

colnames(ZscoreDarmanis_Expression_CellType_NoPrimaryOverlap)
table(ZscoreDarmanis_Expression_CellType_NoPrimaryOverlap[,12])

#Cutting out the darmanis data:
astrocyte_genes<-ZscoreDarmanis_Expression_CellType_NoPrimaryOverlap[ZscoreDarmanis_Expression_CellType_NoPrimaryOverlap[,14]=="Astrocyte" & ZscoreDarmanis_Expression_CellType_NoPrimaryOverlap[,12]!="Darmanis_PNAS_2015",4]

astrocyte_reference<-marker(expr=JoinedDarmanis_AsNumMatrix_Log2, id=astrocyte_genes)

png("PSEA_AstrocyteVsCellType_NoDarmanis.png", width=1000, height=400) 
boxplot(astrocyte_reference~ExperimentalCellType, ylab="Astrocyte Reference Signal", col="lavender")
dev.off()

summary.lm(lm(astrocyte_reference~ExperimentalCellType))

# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept)                            1.98560    0.08984  22.102  < 2e-16 ***
#   ExperimentalCellTypeendothelial       -0.96444    0.18191  -5.302 1.79e-07 ***
#   ExperimentalCellTypefetal_quiescent   -1.53850    0.11234 -13.695  < 2e-16 ***
#   ExperimentalCellTypefetal_replicating -1.19230    0.16759  -7.114 4.37e-12 ***
#   ExperimentalCellTypehybrid            -0.60799    0.13766  -4.417 1.25e-05 ***
#   ExperimentalCellTypemicroglia         -1.65073    0.19836  -8.322 1.01e-15 ***
#   ExperimentalCellTypeneurons           -0.85153    0.10918  -7.799 4.26e-14 ***
#   ExperimentalCellTypeoligodendrocytes  -1.35235    0.14457  -9.354  < 2e-16 ***
#   ExperimentalCellTypeOPC               -1.28549    0.18940  -6.787 3.56e-11 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 0.7074 on 457 degrees of freedom
# Multiple R-squared:  0.3399,	Adjusted R-squared:  0.3284 
# F-statistic: 29.42 on 8 and 457 DF,  p-value: < 2.2e-16

#The r-squared values are equivalent if you leave in the sd0 genes or take them out.

astrocyte_reference<-marker(expr=JoinedDarmanis_AsNumMatrix_Log2_NoSD0, id=astrocyte_genes)

png("PSEA_AstrocyteVsCellType_noSD0_NoDarmanis.png", width=1000, height=400) 
boxplot(astrocyte_reference~ExperimentalCellType, ylab="Astrocyte Reference Signal", col="lavender")
dev.off()

summary.lm(lm(astrocyte_reference~ExperimentalCellType))

# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept)                            1.98560    0.08984  22.102  < 2e-16 ***
#   ExperimentalCellTypeendothelial       -0.96444    0.18191  -5.302 1.79e-07 ***
#   ExperimentalCellTypefetal_quiescent   -1.53850    0.11234 -13.695  < 2e-16 ***
#   ExperimentalCellTypefetal_replicating -1.19230    0.16759  -7.114 4.37e-12 ***
#   ExperimentalCellTypehybrid            -0.60799    0.13766  -4.417 1.25e-05 ***
#   ExperimentalCellTypemicroglia         -1.65073    0.19836  -8.322 1.01e-15 ***
#   ExperimentalCellTypeneurons           -0.85153    0.10918  -7.799 4.26e-14 ***
#   ExperimentalCellTypeoligodendrocytes  -1.35235    0.14457  -9.354  < 2e-16 ***
#   ExperimentalCellTypeOPC               -1.28549    0.18940  -6.787 3.56e-11 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 0.7074 on 457 degrees of freedom
# Multiple R-squared:  0.3399,	Adjusted R-squared:  0.3284 
# F-statistic: 29.42 on 8 and 457 DF,  p-value: < 2.2e-16

#What if we get rid of the hybrid and fetal cells (since they are ambiguous anyway):
summary.lm(lm(astrocyte_reference[ExperimentalCellType!="fetal_quiescent" & ExperimentalCellType!="fetal_replicating"&ExperimentalCellType!="hybrid"]~ExperimentalCellType[ExperimentalCellType!="fetal_quiescent" & ExperimentalCellType!="fetal_replicating"&ExperimentalCellType!="hybrid"]))
# Residual standard error: 0.7865 on 279 degrees of freedom
# Multiple R-squared:  0.2833,	Adjusted R-squared:  0.2705 
# F-statistic: 22.06 on 5 and 279 DF,  p-value: < 2.2e-16
#Huh - I guess that makes sense. 


#Making a version of our analysis that resembles theirs:
AstrocyteIndex_NoDarmanis_NoTagAvg<-apply(ZscoreDarmanis_Expression_CellType_NoPrimaryOverlap[ZscoreDarmanis_Expression_CellType_NoPrimaryOverlap[,14]=="Astrocyte"& ZscoreDarmanis_Expression_CellType_NoPrimaryOverlap[,12]!="Darmanis_PNAS_2015",c(15:ncol(ZscoreDarmanis_Expression_CellType_NoPrimaryOverlap))], 2, mean)

png("Us_AstrocyteVsCellType_noSD0_NoDarmanis_noTagAvg.png", width=1000, height=400) 
boxplot(AstrocyteIndex_NoDarmanis_NoTagAvg~ExperimentalCellType, ylab="Astrocyte Index", col="lavender")
dev.off()

summary.lm(lm(AstrocyteIndex_NoDarmanis_NoTagAvg~ExperimentalCellType))
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept)                            0.43802    0.02617  16.735  < 2e-16 ***
#   ExperimentalCellTypeendothelial       -0.55742    0.05300 -10.518  < 2e-16 ***
#   ExperimentalCellTypefetal_quiescent   -0.65972    0.03273 -20.157  < 2e-16 ***
#   ExperimentalCellTypefetal_replicating -0.53039    0.04883 -10.863  < 2e-16 ***
#   ExperimentalCellTypehybrid            -0.27360    0.04010  -6.822 2.86e-11 ***
#   ExperimentalCellTypemicroglia         -0.71832    0.05779 -12.430  < 2e-16 ***
#   ExperimentalCellTypeneurons           -0.39793    0.03181 -12.510  < 2e-16 ***
#   ExperimentalCellTypeoligodendrocytes  -0.55830    0.04212 -13.255  < 2e-16 ***
#   ExperimentalCellTypeOPC               -0.53108    0.05518  -9.625  < 2e-16 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 0.2061 on 457 degrees of freedom
# Multiple R-squared:  0.5224,	Adjusted R-squared:  0.514 
# F-statistic: 62.48 on 8 and 457 DF,  p-value: < 2.2e-16

#Our method still does way better.

#What happens then if we average by tag first?

ZscoreDarmanis_Expression_CellType_NoPrimaryOverlap_MeanTag_NoDarmanis<-matrix(0, length(names(table(ZscoreDarmanis_Expression_CellType_NoPrimaryOverlap[ZscoreDarmanis_Expression_CellType_NoPrimaryOverlap[,12]!="Darmanis_PNAS_2015",13]))), ncol(ZscoreDarmanis_Expression_CellType_NoPrimaryOverlap)-14)

row.names(ZscoreDarmanis_Expression_CellType_NoPrimaryOverlap_MeanTag_NoDarmanis)<-names(table(ZscoreDarmanis_Expression_CellType_NoPrimaryOverlap[ZscoreDarmanis_Expression_CellType_NoPrimaryOverlap[,12]!="Darmanis_PNAS_2015",13]))
colnames(ZscoreDarmanis_Expression_CellType_NoPrimaryOverlap_MeanTag_NoDarmanis)<-colnames(ZscoreDarmanis_Expression_CellType_NoPrimaryOverlap)[c(15:ncol(ZscoreDarmanis_Expression_CellType_NoPrimaryOverlap))]

for(i in 15:ncol(ZscoreDarmanis_Expression_CellType_NoPrimaryOverlap)){
  
ZscoreDarmanis_Expression_CellType_NoPrimaryOverlap_MeanTag_NoDarmanis[,i-14]<-tapply(ZscoreDarmanis_Expression_CellType_NoPrimaryOverlap[ZscoreDarmanis_Expression_CellType_NoPrimaryOverlap[,12]!="Darmanis_PNAS_2015",i], ZscoreDarmanis_Expression_CellType_NoPrimaryOverlap[ZscoreDarmanis_Expression_CellType_NoPrimaryOverlap[,12]!="Darmanis_PNAS_2015",13] , mean)

}

TagsNoDarmanis<-as.matrix(row.names(ZscoreDarmanis_Expression_CellType_NoPrimaryOverlap_MeanTag_NoDarmanis))

colnames(TagsNoDarmanis)<-"Tag"

CellTypePrimaryVsTag_NoDarmanis<-join(as.data.frame(TagsNoDarmanis), as.data.frame(CellTypePrimaryVsTag), by="Tag", type="left")

CellTypeIndices_NoDarmanis<-matrix(0, 10, ncol(ZscoreDarmanis_Expression_CellType_NoPrimaryOverlap_MeanTag_NoDarmanis))
row.names(CellTypeIndices_NoDarmanis)<-names(table(CellTypePrimaryVsTag_NoDarmanis[,2]))
colnames(CellTypeIndices_NoDarmanis)<-colnames(ZscoreDarmanis_Expression_CellType_NoPrimaryOverlap_MeanTag_NoDarmanis)

for(i in 1:ncol(ZscoreDarmanis_Expression_CellType_NoPrimaryOverlap_MeanTag_NoDarmanis)){
  
  CellTypeIndices_NoDarmanis[,i]<-tapply(ZscoreDarmanis_Expression_CellType_NoPrimaryOverlap_MeanTag_NoDarmanis[,i], CellTypePrimaryVsTag_NoDarmanis[,2], mean)
  
}

png("Us_AstrocyteVsCellType_noSD0_NoDarmanis_TagAvg.png", width=1000, height=400) 
boxplot(CellTypeIndices_NoDarmanis[1,]~ExperimentalCellType, ylab="Astrocyte Index", col="lavender")
dev.off()

summary.lm(lm(CellTypeIndices_NoDarmanis[1,]~ExperimentalCellType))
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept)                            0.57208    0.02860  20.001   <2e-16 ***
#   ExperimentalCellTypeendothelial       -0.70784    0.05791 -12.222   <2e-16 ***
#   ExperimentalCellTypefetal_quiescent   -0.81133    0.03577 -22.685   <2e-16 ***
#   ExperimentalCellTypefetal_replicating -0.64690    0.05336 -12.124   <2e-16 ***
#   ExperimentalCellTypehybrid            -0.40692    0.04383  -9.285   <2e-16 ***
#   ExperimentalCellTypemicroglia         -0.80703    0.06315 -12.779   <2e-16 ***
#   ExperimentalCellTypeneurons           -0.57574    0.03476 -16.563   <2e-16 ***
#   ExperimentalCellTypeoligodendrocytes  -0.72930    0.04603 -15.845   <2e-16 ***
#   ExperimentalCellTypeOPC               -0.67194    0.06030 -11.144   <2e-16 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 0.2252 on 457 degrees of freedom
# Multiple R-squared:  0.5672,	Adjusted R-squared:  0.5596 
# F-statistic: 74.87 on 8 and 457 DF,  p-value: < 2.2e-16

#So the best fit uses our averaging technique, and averages by cell type tag first.


#Oligodendrocytes:
#Cutting out the darmanis data:
oligodendrocyte_genes<-ZscoreDarmanis_Expression_CellType_NoPrimaryOverlap[ZscoreDarmanis_Expression_CellType_NoPrimaryOverlap[,14]=="Oligodendrocyte" & ZscoreDarmanis_Expression_CellType_NoPrimaryOverlap[,12]!="Darmanis_PNAS_2015",4]

oligodendrocyte_reference<-marker(expr=JoinedDarmanis_AsNumMatrix_Log2_NoSD0, id=oligodendrocyte_genes)

png("PSEA_OligodendrocyteVsCellType_noSD0_NoDarmanis.png", width=1000, height=400) 
boxplot(oligodendrocyte_reference~ExperimentalCellType, ylab="Oligodendrocyte Reference Signal", col="lightpink")
dev.off()

summary.lm(lm(oligodendrocyte_reference~ExperimentalCellType))
# Residual standard error: 0.5473 on 457 degrees of freedom
# Multiple R-squared:  0.3769,	Adjusted R-squared:  0.366 
# F-statistic: 34.55 on 8 and 457 DF,  p-value: < 2.2e-16


#Making a version of our analysis that resembles theirs:
OligodendrocyteIndex_NoDarmanis_NoTagAvg<-apply(ZscoreDarmanis_Expression_CellType_NoPrimaryOverlap[ZscoreDarmanis_Expression_CellType_NoPrimaryOverlap[,14]=="Oligodendrocyte"& ZscoreDarmanis_Expression_CellType_NoPrimaryOverlap[,12]!="Darmanis_PNAS_2015",c(15:ncol(ZscoreDarmanis_Expression_CellType_NoPrimaryOverlap))], 2, mean)

png("Us_OligodendrocyteVsCellType_noSD0_NoDarmanis_noTagAvg.png", width=1000, height=400) 
boxplot(OligodendrocyteIndex_NoDarmanis_NoTagAvg~ExperimentalCellType, ylab="Oligodendrocyte Index", col="lightpink")
dev.off()

summary.lm(lm(OligodendrocyteIndex_NoDarmanis_NoTagAvg~ExperimentalCellType))

# Residual standard error: 0.193 on 457 degrees of freedom
# Multiple R-squared:  0.4544,	Adjusted R-squared:  0.4448 
# F-statistic: 47.57 on 8 and 457 DF,  p-value: < 2.2e-16

png("Us_OligodendrocyteVsCellType_noSD0_NoDarmanis_TagAvg.png", width=1000, height=400) 
boxplot(CellTypeIndices_NoDarmanis[8,]~ExperimentalCellType, ylab="Oligodendrocyte Index", col="lightpink")
dev.off()

summary.lm(lm(CellTypeIndices_NoDarmanis[8,]~ExperimentalCellType))
# Residual standard error: 0.1907 on 457 degrees of freedom
# Multiple R-squared:  0.5006,	Adjusted R-squared:  0.4919 
# F-statistic: 57.27 on 8 and 457 DF,  p-value: < 2.2e-16


#Microglia:
#Cutting out the darmanis data:
Microglia_genes<-ZscoreDarmanis_Expression_CellType_NoPrimaryOverlap[ZscoreDarmanis_Expression_CellType_NoPrimaryOverlap[,14]=="Microglia" & ZscoreDarmanis_Expression_CellType_NoPrimaryOverlap[,12]!="Darmanis_PNAS_2015",4]

Microglia_reference<-marker(expr=JoinedDarmanis_AsNumMatrix_Log2_NoSD0, id=Microglia_genes)

png("PSEA_MicrogliaVsCellType_noSD0_NoDarmanis.png", width=1000, height=400) 
boxplot(Microglia_reference~ExperimentalCellType, ylab="Microglia Reference Signal", col="darkolivegreen3")
dev.off()

summary.lm(lm(Microglia_reference~ExperimentalCellType))
# Residual standard error: 1.131 on 457 degrees of freedom
# Multiple R-squared:  0.3608,	Adjusted R-squared:  0.3497 
# F-statistic: 32.25 on 8 and 457 DF,  p-value: < 2.2e-16


#Making a version of our analysis that resembles theirs:
MicrogliaIndex_NoDarmanis_NoTagAvg<-apply(ZscoreDarmanis_Expression_CellType_NoPrimaryOverlap[ZscoreDarmanis_Expression_CellType_NoPrimaryOverlap[,14]=="Microglia"& ZscoreDarmanis_Expression_CellType_NoPrimaryOverlap[,12]!="Darmanis_PNAS_2015",c(15:ncol(ZscoreDarmanis_Expression_CellType_NoPrimaryOverlap))], 2, mean)

png("Us_MicrogliaVsCellType_noSD0_NoDarmanis_noTagAvg.png", width=1000, height=400) 
boxplot(MicrogliaIndex_NoDarmanis_NoTagAvg~ExperimentalCellType, ylab="Microglia Index", col="darkolivegreen3")
dev.off()

summary.lm(lm(MicrogliaIndex_NoDarmanis_NoTagAvg~ExperimentalCellType))

# Residual standard error: 0.1738 on 457 degrees of freedom
# Multiple R-squared:  0.4234,	Adjusted R-squared:  0.4133 
# F-statistic: 41.95 on 8 and 457 DF,  p-value: < 2.2e-16

png("Us_MicrogliaVsCellType_noSD0_NoDarmanis_TagAvg.png", width=1000, height=400) 
boxplot(CellTypeIndices_NoDarmanis[3,]~ExperimentalCellType, ylab="Microglia Index", col="darkolivegreen3")
dev.off()

summary.lm(lm(CellTypeIndices_NoDarmanis[3,]~ExperimentalCellType))
# Residual standard error: 0.2062 on 457 degrees of freedom
# Multiple R-squared:  0.506,	Adjusted R-squared:  0.4973 
# F-statistic: 58.51 on 8 and 457 DF,  p-value: < 2.2e-16

#Endothelial:
#Cutting out the darmanis data:
Endothelial_genes<-ZscoreDarmanis_Expression_CellType_NoPrimaryOverlap[ZscoreDarmanis_Expression_CellType_NoPrimaryOverlap[,14]=="Endothelial" & ZscoreDarmanis_Expression_CellType_NoPrimaryOverlap[,12]!="Darmanis_PNAS_2015",4]

Endothelial_reference<-marker(expr=JoinedDarmanis_AsNumMatrix_Log2_NoSD0, id=Endothelial_genes)

png("PSEA_EndothelialVsCellType_noSD0_NoDarmanis.png", width=1000, height=400) 
boxplot(Endothelial_reference~ExperimentalCellType, ylab="Endothelial Reference Signal", col="orange")
dev.off()

summary.lm(lm(Endothelial_reference~ExperimentalCellType))
# Residual standard error: 0.8934 on 457 degrees of freedom
# Multiple R-squared:  0.3048,	Adjusted R-squared:  0.2927 
# F-statistic: 25.05 on 8 and 457 DF,  p-value: < 2.2e-16


#Making a version of our analysis that resembles theirs:
EndothelialIndex_NoDarmanis_NoTagAvg<-apply(ZscoreDarmanis_Expression_CellType_NoPrimaryOverlap[ZscoreDarmanis_Expression_CellType_NoPrimaryOverlap[,14]=="Endothelial"& ZscoreDarmanis_Expression_CellType_NoPrimaryOverlap[,12]!="Darmanis_PNAS_2015",c(15:ncol(ZscoreDarmanis_Expression_CellType_NoPrimaryOverlap))], 2, mean)

png("Us_EndothelialVsCellType_noSD0_NoDarmanis_noTagAvg.png", width=1000, height=400) 
boxplot(EndothelialIndex_NoDarmanis_NoTagAvg~ExperimentalCellType, ylab="Endothelial Index", col="orange")
dev.off()

summary.lm(lm(EndothelialIndex_NoDarmanis_NoTagAvg~ExperimentalCellType))

# Residual standard error: 0.1865 on 457 degrees of freedom
# Multiple R-squared:  0.2831,	Adjusted R-squared:  0.2706 
# F-statistic: 22.56 on 8 and 457 DF,  p-value: < 2.2e-16

png("Us_EndothelialVsCellType_noSD0_NoDarmanis_TagAvg.png", width=1000, height=400) 
boxplot(CellTypeIndices_NoDarmanis[2,]~ExperimentalCellType, ylab="Endothelial Index", col="orange")
dev.off()

summary.lm(lm(CellTypeIndices_NoDarmanis[2,]~ExperimentalCellType))
# Residual standard error: 0.219 on 457 degrees of freedom
# Multiple R-squared:  0.3284,	Adjusted R-squared:  0.3166 
# F-statistic: 27.93 on 8 and 457 DF,  p-value: < 2.2e-16

#Neuron_All:
#Cutting out the darmanis data:
Neuron_All_genes<-ZscoreDarmanis_Expression_CellType_NoPrimaryOverlap[ZscoreDarmanis_Expression_CellType_NoPrimaryOverlap[,14]=="Neuron_All" & ZscoreDarmanis_Expression_CellType_NoPrimaryOverlap[,12]!="Darmanis_PNAS_2015",4]

Neuron_All_reference<-marker(expr=JoinedDarmanis_AsNumMatrix_Log2_NoSD0, id=Neuron_All_genes)

png("PSEA_Neuron_AllVsCellType_noSD0_NoDarmanis.png", width=1000, height=400) 
boxplot(Neuron_All_reference~ExperimentalCellType, ylab="Neuron_All Reference Signal", col="darkviolet")
dev.off()

summary.lm(lm(Neuron_All_reference~ExperimentalCellType))
# Residual standard error: 0.7891 on 457 degrees of freedom
# Multiple R-squared:  0.3339,	Adjusted R-squared:  0.3222 
# F-statistic: 28.63 on 8 and 457 DF,  p-value: < 2.2e-16


#Making a version of our analysis that resembles theirs:
Neuron_AllIndex_NoDarmanis_NoTagAvg<-apply(ZscoreDarmanis_Expression_CellType_NoPrimaryOverlap[ZscoreDarmanis_Expression_CellType_NoPrimaryOverlap[,14]=="Neuron_All"& ZscoreDarmanis_Expression_CellType_NoPrimaryOverlap[,12]!="Darmanis_PNAS_2015",c(15:ncol(ZscoreDarmanis_Expression_CellType_NoPrimaryOverlap))], 2, mean)

png("Us_Neuron_AllVsCellType_noSD0_NoDarmanis_noTagAvg.png", width=1000, height=400) 
boxplot(Neuron_AllIndex_NoDarmanis_NoTagAvg~ExperimentalCellType, ylab="Neuron_All Index", col="darkviolet")
dev.off()

summary.lm(lm(Neuron_AllIndex_NoDarmanis_NoTagAvg~ExperimentalCellType))

# Residual standard error: 0.2539 on 457 degrees of freedom
# Multiple R-squared:  0.5653,	Adjusted R-squared:  0.5577 
# F-statistic: 74.29 on 8 and 457 DF,  p-value: < 2.2e-16


png("Us_Neuron_AllVsCellType_noSD0_NoDarmanis_TagAvg.png", width=1000, height=400) 
boxplot(CellTypeIndices_NoDarmanis[5,]~ExperimentalCellType, ylab="Neuron_All Index", col="darkviolet")
dev.off()

summary.lm(lm(CellTypeIndices_NoDarmanis[5,]~ExperimentalCellType))
# Residual standard error: 0.251 on 457 degrees of freedom
# Multiple R-squared:  0.5269,	Adjusted R-squared:  0.5186 
# F-statistic: 63.61 on 8 and 457 DF,  p-value: < 2.2e-16


###########################
