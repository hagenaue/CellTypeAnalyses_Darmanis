#These code documents contain all of the code documenting the cell type analyses that we applied to Darmanis' single-cell RNAseq data to validate our methodology
#Megan Hagenauer and Alek Pankonin
#October 27, 2016

#**************************************
#R: 

### R version 3.3.0 (2016-05-03) -- "Supposedly Educational"
### Copyright (C) 2016 The R Foundation for Statistical Computing
### Platform: x86_64-apple-darwin13.4.0 (64-bit)

#I also tend to use R-studio as a GUI:
### Version 0.99.896 – © 2009-2016 RStudio, Inc.
### Mozilla/5.0 (Macintosh; Intel Mac OS X 10_10_5) AppleWebKit/601.5.17 (KHTML, like Gecko)

#****************************************************

#(Potentially) Relevant code libraries:
##I apologize for the fact that I'm not entirely sure which of these I actually used - I just tend to load them all because I use them regularly
##Each of these packages requires installation before they can be loaded.

library(gdata)
library("fields")
library(stats)
library(car)
library(affy)
library(preprocessCore)
library(multtest)
library(plyr)
library("car")
library("drc")
library("gdata")
library("gtools")
library("lattice")
library("lmtest")
library("magic")
library("MASS")
library("plotrix")
library("proto")
library("sandwich")
library("spam")

#************************