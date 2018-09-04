# IVIVE workflow to be used in ICE #################################################################################################
# File name: ivive.R; 
# Original author: Xiaoqing Chang, Date: November 1st, 2017 #
# Code revision number and date: Version_1.0 2018-01-17 ###
# The workflow uses the one-compartment rat PK model that incorporates Monte Carlo simulation to simulate the population variance. The workflow is to predict the daily equivalent administered dose (EAD, mg/kg/day)
# that would lead to steady state blood concentration equivalent to the bioactive concentration from in vitro assays and compared to the predicted lowest effective levels (LELs) of in vivo assays, which is user provided.
# To execute the IVIVE workflow, the required inputs are chemical name, CASRN, in vitro active concentrations, fraction of chemical unbound to protein (fu), intrinsic clearance, and renal clearance.
# The outputs of the workflow are the EAD estimates from in vitro assays with or without in vivo fu adjustment. One example is provided, focusing on assays measuring estrogenic activity.
# In this example, the activity concentration at cut-off (ACC) from estrogen receptor (ER) pathway assays were applied to the workflow to obtain daily EADs that would result in a blood concentration equivalent to the ACC at steady state.
# To validate the results, the EADs were compared to the lowest effective levels in rat uterotrophic assays.#
# R packages required for executing code are stringr, readxl, dplyr, tidyr, scales, grid, ggplot2 #
# R function required for executing code is theme_complete_bw.R #
#wd <- "C:\\working_directory"
#setwd(wd)
# f0 <- '/Input data for Css_rat32_0.25kg (5 chem).xlsx'
# f1 <- '/In vitro assay data (5 chem).xlsx'
# nsamples <- 100 # user-provided, any number between 10 - 10,000
# quantile  <- 1 # user-provided, chosen from 1-5
####input
library("stringr")
rm(list=ls())
date.string <-Sys.Date()
args=(commandArgs(TRUE))
f0 <- args[1]
f1 <- args[2]
f2 <- args[3]
f3 <- args[4]
nsamples <- as.double(args[5]) # user-provided, any number between 10 - 10,000
# #f1 <- '2017447750_input1.txt'
#
# #f0 <- '2017447750_input.txt'
# f0 <- 'Input data for Css_rat32_0.25kg (5 chem).txt'
# f1 <- 'In vitro assay data (5 chem).txt'
# f3 <- 'Outf3l.txt'
# f2 <- 'Outf2l.txt'
# nsamples <- 1000

 #wth      <- 0.2  # the width of the boxplot,user-provided
 #size     <- 5
steadyState <- function(A=f0,inVitro=f1, CssOut=f3){
  #date.string <-Sys.Date()
  library(stringr)
  library(dplyr)
  A <- read.table(A, sep = "\t", header = TRUE)
  inVitro <- read.table(inVitro, sep = "\t", header = TRUE)
  options(stringsAsFactors = FALSE)
  #all.perc <- c("50perc", "95perc")
  #quantile  <- 1 # user-provided, chosen from 1-5 ############ Quantiles
  #perc     <- all.perc[quantile]
  Dose         <- rep(1/24, nrow(A))  # daily dose of 1 mg/kg/day converted to a unit of mg/kg/h
  Fub          <- A$Fub  # unitless
  CLrenalMean  <- A$renal.clearance  # L/h
  CLintrinMean <- A$intrinsic.clearance  # L/h
  MW           <- A$`MW`  # molecular weight
  Css          <- NULL  # Steady state blood concentration
  CssAll       <- NULL  # CssAll is accumulated Css value obtained from each round of calculation
  BWAll        <- NULL  # BW is body weight, BWAll is accumulated BW value obtained from each round of calculation
  QlAll        <- NULL  # Ql is liver blood flow, QlAll is accumulated Ql value obtained from each round of calculation
  CLintrinAll  <- NULL  # CLintrin is intrinsic clearance (L/h), CLintrinAll is accumulated CLintrin value obtained from each round of calculation
  CLrenalAll   <- NULL  # CLrenal is renal clearance (L/h), CLrenalAll is accumulated CLrenal value obtained from each round of calculation
#  std.CssAll   <-  NULL  # std.Css is steady state blood concentration (uM) at 1 mg/kg/day, std.CssAll accumulated std.Css value obtained from each round of calculation
  for (j in 1:nsamples) {           # The P-PK modeling and Monte Carlo simulation
    BW <- rnorm(1, 0.25, 0.25*0.2)  # assume a normal distribution, mean 0.25, std 0.25*0.2
    Ql <- rnorm(1,13750, 13750*0.2)*60/10^6  # assume normal distribution, mean 13750 (uL/min), Sd = 13750*0.2, also change the unit from ul/min to L/h
    CLrenal <- rnorm(length(CLrenalMean), CLrenalMean, CLrenalMean*0.2)  # assume the coefficeient of variance (cv) is 0.2 for renal clearance, cv=sd/mean
    CLintrin <- rnorm(length(CLintrinMean), CLintrinMean, CLintrinMean*0.2)  # assume the coefficeient of variance (cv) is 0.2 for intrinsic clearance, cv=sd/mean
    Css <- Dose*BW/(CLrenal+(Ql*Fub*CLintrin)/(Ql+Fub*CLintrin))  # unit is mg/L
    Css <- Css*1000/MW  # convert unit from mg/L to uM
    CssAll <- cbind(CssAll,Css)
    BWAll[j] <- BW
    QlAll[j] <- Ql
    CLrenalAll <- cbind(CLrenalAll,CLrenal)
    CLintrinAll <- cbind(CLintrinAll,CLintrin)
  }
  CssAll <- t(apply(CssAll, 1, quantile, probs= c(.5, .95))) ##get the 50th and 95th percentile for Css at 1mg/kg/day
  CssAll <- as.data.frame(CssAll)
  CssAlls <- cbind(A, CssAll)
  write.table(CssAlls,file=CssOut, sep="\t",row.names=F)
  CssAlls
}
invitroacc <- function(CssAlls=CssAlls, inVitro=f1, EADOut=f2){
  #date.string <-Sys.Date()
  library(stringr)
  library(dplyr)
  #library(tidyr)
  #library(scales)
  options(stringsAsFactors = FALSE)
  #all.perc <- c("50perc", "95perc")
  # quantile  <- 1
  # perc     <- all.perc[quantile]
  inVitro <- read.table( inVitro, sep = "\t", header = TRUE)
  inVitro.m <- merge(CssAlls, inVitro, by= c("CASRN"), all = TRUE) ## merge() two datasets based on c("column Name")
  # std.Css.5perc <- inVitro.m$`5%`
  # std.Css.25perc <- inVitro.m$`25%`
  std.Css.50perc <- inVitro.m$`50%`
  # std.Css.75perc <- inVitro.m$`75%`
  std.Css.95perc <- inVitro.m$`95%`
  # std.CssAll <- cbind(std.Css.50perc, std.Css.95perc)
  # std.Css <- std.CssAll[,quantile]
  #std.Css <- std.Css.50perc[,quantile]
  colnames(inVitro.m)[2] <- "Chemical"
  i<-ncol(CssAlls)+2
  # To calculate EAD.total and EAD.adj corresponding to in vitro concentation of each assay
  if(nrow(inVitro.m) ==1){
    EAD.total50 <- t(sapply(inVitro.m[,i:ncol(inVitro.m)], function(x) x/std.Css.50perc))
    EAD.total95 <- t(sapply(inVitro.m[,i:ncol(inVitro.m)], function(x) x/std.Css.95perc))
    colnames(EAD.total50) <- paste("EAD.total50", colnames(EAD.total50), sep = "_")
    colnames(EAD.total95) <- paste("EAD.total95", colnames(EAD.total95), sep = "_")
  } else {
    EAD.total50 <- sapply(inVitro.m[,i:ncol(inVitro.m)], function(x) x/std.Css.50perc)
    EAD.total95 <- sapply(inVitro.m[,i:ncol(inVitro.m)], function(x) x/std.Css.95perc)
    colnames(EAD.total50) <- paste("EAD.total50", colnames(EAD.total50), sep = "_")
    colnames(EAD.total95) <- paste("EAD.total95", colnames(EAD.total95), sep = "_")
  }


  #calculate the EAD ajusted by fraction unbound
  Fub2 <- inVitro.m$Fub
  if(nrow(inVitro.m) ==1){
    EAD.adj50 <- t(sapply(inVitro.m[,i:ncol(inVitro.m)], function(x) x/std.Css.50perc/Fub2))
    EAD.adj95 <- t(sapply(inVitro.m[,i:ncol(inVitro.m)], function(x) x/std.Css.95perc/Fub2))
    colnames(EAD.adj50) <- paste("EAD.adj50", colnames(EAD.adj50), sep = "_")
    colnames(EAD.adj95) <- paste("EAD.adj95", colnames(EAD.adj95), sep = "_")
  } else {
    EAD.adj50 <- sapply(inVitro.m[,i:ncol(inVitro.m)], function(x) x/std.Css.50perc/Fub2)
    EAD.adj95 <- sapply(inVitro.m[,i:ncol(inVitro.m)], function(x) x/std.Css.95perc/Fub2)
    colnames(EAD.adj50) <- paste("EAD.adj50", colnames(EAD.adj50), sep = "_")
    colnames(EAD.adj95) <- paste("EAD.adj95", colnames(EAD.adj95), sep = "_")
  }
  # Output <- cbind(inVitro.m, EAD.total, EAD.adj)
  Output <- cbind(inVitro.m, EAD.total50, EAD.total95, EAD.adj50, EAD.adj95)
  write.table(Output,file=EADOut, sep="\t",row.names=F,col.names=T)
  return(Output)
}

CssAlls<-steadyState(f0, f1, f3)

Eads<-invitroacc(CssAlls, f1, f2)

# A <- read.table(f0, sep = "\t", header = TRUE)
# inVitro <- read.table(f1, sep = "\t", header = TRUE)
# CMP1pk<-invitroacc(CssAlls=steadyState(A,inVitro), inVitro = inVitro)

