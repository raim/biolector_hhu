####Set working directory, files and IDs####

setwd("C:/Users/Katrin/biolector_hhu/data/20191121_BIQ980_biolector_pre2")
cfile <- read.csv("C:/Users/Katrin/biolector_hhu/data/calibrations/calibrations.csv")
parameters <- read.csv("C:/Users/Katrin/biolector_hhu/data/Parameters.csv", sep = ";")

expID = "21_BIQ980_2019_REDUCTION-1"
calID <- "20210315_calib_split_pre"


####preparation####

if (TRUE==TRUE){
  library(platexpress)
  
  dfile <- paste(expID, ".csv", sep="")
  lfile <- paste(expID, "layout.csv", sep = "_")
  map_overview <- read.csv(lfile, sep = ";")
  
  if (length(grep("Glc:",map_overview[1,2])>0 ) & length(grep("Ace:",map_overview[1,2])>0 )){
    id <- "GlcAce"
  } else if (length(grep("Ace:",map_overview[1,2])>0 )){
    id <- "Ace"
  } else if (length(grep("Glc:",map_overview[1,2])>0 )){
    id <- "Glc"
  } else if (length(grep("aTc:",map_overview[1,2])>0 )){
    id <- "aTc"
  }
}


####analysis####

if (id == "GlcAce"){
  #1: read plate map
  map_glc <- readPlateMap(lfile, sep = ";", fsep = ";", asep = ":",
                          fields = c("strain", "Glc", "Ace", "aTc"),
                          afields = c("Glc", "Ace", "aTc"))
  
  map_ace <- readPlateMap(lfile, sep = ";", fsep = ";", asep = ":",
                          fields = c("strain", "Glc", "Ace", "aTc"),
                          afields = c("Ace", "Glc", "aTc"))
  #1.2: initial carbon amount per well
  carbon_mol <- NULL
  for (i in 1:length(map_glc$well)){
    carboni <- (map_glc$amount[i]/100*1000)/parameters$Number[which(parameters$Parameter=="Glucose MW")]*
      parameters$Number[which(parameters$Parameter=="Glucose NC")]+
      (map_glc$Ace.amount[i]/100*1000)/parameters$Number[which(parameters$Parameter=="Glucose MW")]*
      parameters$Number[which(parameters$Parameter=="Glucose NC")]
    carbon_mol <- append(carbon_mol, carboni)
  }
  names(carbon_mol) <- map_glc$well
  map_carbon <- cbind(map_glc, carbon = "carbon", carbon.amount=round(carbon_mol, 3))
  map_carbon <- amountColors(map_carbon, amount = "carbon.amount", substance = "carbon", color = "carbon.color")
  
  
  #2: read plate data
  dat <- readPlateData(dfile,
                       type="BioLectorPro", time.conversion = 1/3600)
  
  #3: check for blanks
  if (any(map_carbon$blank)){
    #3: read experiment
    data <- readExperiment(dfile,
                           type = "BioLectorPro",
                           time.conversion = 1/3600,
                           layout = lfile, 
                           fields = c("strain","Glc","Ace", "aTc"), 
                           afields = c("Glc", "Ace", "aTc"),
                           group1 = c("strain"),
                           group2 = c("strain","amount", "Ace.amount"),
                           group2.color="color",
                           blank.id = "blank",
                           blank.data = c("Biomass","Riboflavine",
                                          "NADH - NADPH"))
    #4: rename data
    data <- prettyData(data, yids=c(scatter="Biomass", "NAD(P)H"="NADH - NADPH",
                                    DO="DO(Pst3)", pH="pH(HP8)",
                                    "riboflavin"="Riboflavine"))
    #5: correct pH trend
    data <- correctBlanks(data, yids="pH", mbins=100)
    
  } else {
    data <- readExperiment(dfile,
                           type = "BioLectorPro",
                           time.conversion = 1/3600,
                           layout = lfile, 
                           fields = c("strain","Glc","Ace", "aTc"), 
                           afields = c("Glc", "Ace", "aTc"),
                           group1 = c("strain"),
                           group2 = c("strain","amount","Ace.amount"),
                           group2.color="color")
    #4: rename data
    data <- prettyData(data, yids=c(scatter="Biomass", "NAD(P)H"="NADH - NADPH",
                                    DO="DO(Pst3)", pH="pH(HP8)",
                                    "riboflavin"="Riboflavine"))
  }
  data$pH$data <- data$pH$data + 7.2
  
  data$groups$group2.color <- 
    groupColors(map_carbon,
                group=getGroups(map_carbon,
                                by=c("strain","amount","Ace.amount")),
                color="carbon.color")
  
  #7: calibrate OD to biomass [g/L]
  # y = m*x+b ; m -> slope, b -> intercept
  calscatter <- which(cfile$ID==calID)[1]
  calriboflavin <- which(cfile$ID==calID)[2]
  scatter_m <- cfile$slope[calscatter]
  scatter_b <- cfile$intercept[calscatter]
  ribo_m <- cfile$slope[calriboflavin]
  ribo_b <- cfile$intercept[calriboflavin]
  
  data <- addData(data=data, dat=scatter_m*getData(data, "scatter")+scatter_b,
                  ID="biomass_scatter", replace=TRUE)
  data <- addData(data=data, dat=ribo_m*getData(data, "riboflavin")+ribo_b,
                  ID="biomass_riboflavin", replace=TRUE)
  data <- addData(data, ID="biomass_carbon", dat = getData(data, "biomass_scatter")*0.474/12)
  data <- addData(data, ID="biomass_carbon_cmmol", dat=getData(data, "biomass_carbon")*1000)
  #calculate yield
  carbon_mol <- carbon_mol[colnames(data$biomass_carbon$data)] #reorder vector
  yield <- NULL
  for (i in 1:length(data$wells$plate)){
    yi <- data$biomass_carbon$data[,i]/carbon_mol[[i]]
    yield <- cbind(yield, yi)
  }
  colnames(yield) <- data$wells$plate
  data <- addData(data, ID="yield", dat=yield)
  
  #8: dpseg
  dpp <- platexpress::dpseg_plate(data, yid="biomass_scatter", P=0.001)
  
  #add dpseg piecewise linear model to platexpress object
  data <- addModel(fit=dpp, data=data, ID="dpseg", replace=TRUE)
  data <- addModel(fit=dpp, data=data, ID="rates",add.slopes=TRUE, replace=TRUE)
  
  #9: plots
  pdf(paste0(expID,"_report.pdf"))
  viewMap(map_glc, title = "plate layout: glucose")
  viewMap(map_ace, title = "plate layout: acetate")
  viewMap(map_carbon, text = "carbon.amount", color = "carbon.color", title = "plate layout: carbon [cmol]")
  viewPlate(dat, legpos = "bottomright", add.legend = TRUE)
  viewPlate(data, legpos = "bottomright", add.legend = TRUE)
  par(mfrow=c(1,1))
  viewGroups(data, stats="range",yids="scatter",
             g2.legend=FALSE, g1.legend = FALSE,lwd.orig=0,
             xlab = "Time [h]", ylab = "scatter [a.u.]",
             mai = c(.5,.5,0,0))
  viewGroups(data, stats="range",yids="riboflavin",
             g2.legend=FALSE, g1.legend = FALSE,lwd.orig=0,
             xlab = "Time [h]", ylab = "riboflavine fluorescence [a.u.]",
             mai = c(.5,.5,0,0))
  viewGroups(data, stats="range",yids="DO",
             g2.legend=FALSE, g1.legend = FALSE,lwd.orig=0,
             xlab = "Time [h]", ylab = "Diss. O2 [%]",
             mai = c(.5,.5,0,0))
  viewGroups(data, stats="range",yids="pH",
             g2.legend=FALSE, g1.legend = FALSE,lwd.orig=0,
             xlab = "Time [h]", ylab = "pH",
             mai = c(.5,.5,0,0)) 
  viewGroups(data, stats="range",yids="NAD(P)H",
             g2.legend=FALSE, g1.legend = FALSE,lwd.orig=0,
             xlab = "Time [h]", ylab = "NAD(P)H fluorescence [a.u.]",
             mai = c(.5,.5,0,0))
  viewGroups(data, stats="range",yids="biomass_scatter",
             g2.legend=FALSE, g1.legend = FALSE,lwd.orig=0,
             xlab = "Time [h]", ylab = "biomass [g/L]",
             mai = c(.5,.5,0,0))
  viewGroups(data,stats="range",yids="dpseg",log="y",
             lwd.orig=0,g2.legend=FALSE, g1.legend = FALSE,
             xlab = "Time [h]", ylab = "dpseg",
             mai = c(.5,.5,0,0))
  viewGroups(data,stats="blub",yids="rates",log="",
             lwd.orig=0,g2.legend=FALSE, g1.legend = FALSE,
             xlab = "Time [h]", ylab =expression("growth rate, h"^-1),
             mai = c(.5,.5,0,0))
  viewGroups(data, yids = "biomass_carbon_cmmol", stats = "range",
             g2.legend = FALSE, g1.legend = FALSE, lwd.orig = 0,
             xlab = "Time [h]", ylab = "carbon amount of biomass [cmmol]",
             mai = c(.5,.5,0,0))
  viewGroups(data, yids = "biomass_carbon_cmmol", stats = "range", log = "y",
             g2.legend = FALSE, g1.legend = FALSE, lwd.orig = 0,
             xlab = "Time [h]", ylab = "carbon amount of biomass [cmmol]",
             mai = c(.5,.5,0,0))
  viewGroups(data, yids="yield", stats = "mean",
             g2.legend=FALSE, g1.legend = FALSE,lwd.orig=0,
             xlab = "Time [h]", ylab = "yield",
             mai = c(.5,.5,0,0))
  dev.off()
  
} else if ( id=="aTc" ){
  #1: read plate map
  map <- readPlateMap(lfile, sep = ";", fsep = ";", asep = ":",
                      fields = c("strain", "medium", "aTc"),
                      afields = c("aTc"))
  
  #2: read plate data
  dat <- readPlateData(dfile,
                       type="BioLectorPro", time.conversion = 1/3600, dec=",")
  
  #3: check for blanks
  if (any(map_carbon$blank)){
    #3: read experiment
    data <- readExperiment(dfile,
                           type = "BioLectorPro",
                           time.conversion = 1/3600,
                           dec = ",",
                           layout = lfile, 
                           fields = c("strain","medium","aTc"), 
                           afields = c("aTc"),
                           group1 = c("strain"),
                           group2 = c("strain","amount"),
                           group2.color="color",
                           blank.id = "blank",
                           blank.data = c("Biomass","Riboflavine",
                                          "NADH - NADPH"))
    #4: rename data
    data <- prettyData(data, yids=c(scatter="Biomass", "NAD(P)H"="NADH - NADPH",
                                    DO="DO(Pst3)", pH="pH(HP8)",
                                    "riboflavin"="Riboflavine"))
    #5: correct pH trend
    data <- correctBlanks(data, yids="pH", mbins=100)
    
  } else {
    data <- readExperiment(dfile,
                           type = "BioLectorPro",
                           time.conversion = 1/3600,
                           dec = ",",
                           layout = lfile, 
                           fields = c("strain","medium","aTc"), 
                           afields = c("aTc"),
                           group1 = c("strain"),
                           group2 = c("strain","amount"),
                           group2.color="color")
    #4: rename data
    data <- prettyData(data, yids=c(scatter="Biomass", "NAD(P)H"="NADH - NADPH",
                                    DO="DO(Pst3)", pH="pH(HP8)",
                                    "riboflavin"="Riboflavine"))
  }
  
  ## correct aTc concentration
  data$layout$amount <- 150*data$layout$amount/max(data$layout$amount)
  ##correct pH
  data$pH$data <- data$pH$data + 7.2
  
  #7: calibrate OD to biomass
  ## y = m*x+b ; m -> slope, b -> intercept
  calscatter <- which(cfile$ID==calID)[1]
  calriboflavin <- which(cfile$ID==calID)[2]
  scatter_m <- cfile$slope[calscatter]
  scatter_b <- cfile$intercept[calscatter]
  ribo_m <- cfile$slope[calriboflavin]
  ribo_b <- cfile$intercept[calriboflavin]
  
  data <- addData(data=data, dat=scatter_m*getData(data, "scatter")+scatter_b,
                  ID="biomass_scatter", replace=TRUE)
  data <- addData(data=data, dat=ribo_m*getData(data, "riboflavin")+ribo_b,
                  ID="biomass_riboflavin", replace=TRUE)
  
  #8: dpseg
  dpp <- platexpress::dpseg_plate(data, yid="biomass_scatter", P=0.001)
  
  #add dpseg piecewise linear model to platexpress object
  data <- addModel(fit=dpp, data=data, ID="dpseg", replace=TRUE)
  data <- addModel(fit=dpp, data=data, ID="rates",add.slopes=TRUE, replace=TRUE)
  
  ##prepare for dose-response
  res <- data$layout
  tmp <- getData(data, ID="biomass_scatter", xrng = 50)
  res <- cbind(res, biomass=tmp[match(as.character(res$well),
                                      names(tmp))])
  
  #9: plots
  pdf(paste0(expID,"_report.pdf"))
  viewMap(map)
  viewPlate(dat, legpos = "bottomright", add.legend = TRUE)
  viewPlate(data, legpos = "bottomright", add.legend = TRUE)
  par(mfrow=c(1,1))
  viewGroups(data, stats="range",yids="scatter",
             g2.legend=FALSE, g1.legend = FALSE,lwd.orig=0,
             xlab = "Time [h]", ylab = "scatter [a.u.]",
             mai = c(.5,.5,0,0))
  viewGroups(data, stats="range",yids="riboflavin",
             g2.legend=FALSE, g1.legend = FALSE,lwd.orig=0,
             xlab = "Time [h]", ylab = "riboflavine fluorescence [a.u.]",
             mai = c(.5,.5,0,0))
  viewGroups(data, stats="range",yids="DO",
             g2.legend=FALSE, g1.legend = FALSE,lwd.orig=0,
             xlab = "Time [h]", ylab = "Diss. O2 [%]",
             mai = c(.5,.5,0,0))
  viewGroups(data, stats="range",yids="pH",
             g2.legend=FALSE, g1.legend = FALSE,lwd.orig=0,
             xlab = "Time [h]", ylab = "pH",
             mai = c(.5,.5,0,0)) 
  viewGroups(data, stats="range",yids="NAD(P)H",
             g2.legend=FALSE, g1.legend = FALSE,lwd.orig=0,
             xlab = "Time [h]", ylab = "NAD(P)H fluorescence [a.u.]",
             mai = c(.5,.5,0,0))
  viewGroups(data, stats="range",yids="biomass_scatter",
             g2.legend=FALSE, g1.legend = FALSE,lwd.orig=0,
             xlab = "Time [h]", ylab = "biomass [g/L]",
             mai = c(.5,.5,0,0))
  viewGroups(data,stats="range",yids="dpseg",log="y",
             lwd.orig=0,g2.legend=FALSE, g1.legend = FALSE,
             xlab = "Time [h]", ylab = "dpseg",
             mai = c(.5,.5,0,0))
  viewGroups(data,stats="blub",yids="rates",log="",
             lwd.orig=0,g2.legend=FALSE, g1.legend = FALSE,
             xlab = "Time [h]", ylab = expression("growth rate, h"^-1),
             mai = c(.5,.5,0,0))
  doseResponse(res, val="biomass", amount="amount",
               wells=data$groups$group1[[1]],col="color",
               ylab="biomass [g/L]", xlab="aTc [ng/mL]") 
  dev.off()
}
