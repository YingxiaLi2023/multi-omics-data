########################################################

# Set the working directory to the directory 'multi-omics-data' 
# of the electronic appendix (outcomment the following line
# and replace 'pathtomulti-omics-data' by the path to 'multi-omics-data'
# on your computer):

## setwd("pathtomulti-omics-data/multi-omics-data/Data")

########################################################


# OpenML dataset ids to querry
library(OpenML)
library(mlr)
library(ParamHelpers)
library(farff)



# --> This file can be downloaded from the Github page in
# the subfolder "data":

# Redownload as survival information was not included the first time
load("./datset_ids.RData")



#26-3-5-1(BRCA)
nams <- c( "BLCA", "COAD", "ESCA", "HNSC", 
           "KIRC", "KIRP", "LAML","LIHC", "LGG",
          "LUAD", "LUSC",  "OV", "PAAD", "SARC", 
          "SKCM", "STAD", "UCEC") #, "PCPG", "PRAD", "TGCT","THCA", "THYM", "BRCA",


nam <- "PRAD"
for(nam in nams){
  # download dataset
  dat_part1 <- getOMLDataSet(datset_ids[[nam]][[1]])
  dat_part2 <- getOMLDataSet(datset_ids[[nam]][[2]])
  
  dat <- cbind.data.frame(dat_part1, dat_part2)
  
  
  blocknames <- c("clinical", "cnv", "mirna", "mutation", "rna")
  
  blockinds <- lapply(paste0("_", blocknames), function(x) grep(x, names(dat)))
  # --> blockinds is list of length 5, where the first list element contains the indices
  # of the clinical data, the second list element that of the cnv data
  # and so on (see "blocknames" above).
  
  
  clindata <- dat[,blockinds[[1]]]
  cnvdata <- dat[,blockinds[[2]]]
  mirnadata <- dat[,blockinds[[3]]]
  mutationdata <- dat[,blockinds[[4]]]
  rnadata <- dat[,blockinds[[5]]]
  surdata <- dat[,1:3]
  
  #head(names(dat))
  # --> "bcr_patient_barcode", "time", und "status" do not belong to the covariates
  # and have to be removed.
  
  save(clindata, cnvdata, mirnadata,mutationdata, rnadata,surdata,
       file = paste(nam,".RData", sep = ""))
}

#load("E:/18_omics_datas/UCEC.RData")

#########download "BRCA", due to "BRCA" has 3 parts

nam <- "BRCA"


dat_part1 <- getOMLDataSet(datset_ids[[nam]][[1]])
dat_part2 <- getOMLDataSet(datset_ids[[nam]][[2]])

dat <- cbind.data.frame(dat_part1, dat_part2)

dat_part3 <- getOMLDataSet(datset_ids[[nam]][[3]])
dat <- cbind.data.frame(dat, dat_part3)


blocknames <- c("clinical", "cnv", "mirna", "mutation", "rna")

blockinds <- lapply(paste0("_", blocknames), function(x) grep(x, names(dat)))
# --> blockinds is list of length 5, where the first list element contains the indices
# of the clinical data, the second list element that of the cnv data
# and so on (see "blocknames" above).


clindata <- dat[,blockinds[[1]]]
cnvdata <- dat[,blockinds[[2]]]
mirnadata <- dat[,blockinds[[3]]]
mutationdata <- dat[,blockinds[[4]]]
rnadata <- dat[,blockinds[[5]]]
surdata <- dat[,1:3]
#head(names(dat))
# --> "bcr_patient_barcode", "time", und "status" do not belong to the covariates
# and have to be removed.

save(clindata, cnvdata, mirnadata,mutationdata, rnadata, surdata,
     file = paste(nam,".RData", sep = ""))
