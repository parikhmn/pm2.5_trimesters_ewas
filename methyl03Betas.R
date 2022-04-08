# Methylation Data Analysis on β values
# Milan Parikh

knitr::opts_knit$set(root.dir = "../")
setwd("/path/to/outputFolder1/")

##### Import Libraries and Data #####
suppressMessages(library(minfi))

load("WB.noob2.RData")
load("betas.rcp_norm.RData")

setwd("/path/to/outputFolder3/")

suppressPackageStartupMessages({
  library(CpGassoc) # for running association analysis between methylation levels values and phenotype of interest
  library(data.table) # for fast aggregation of large data 
  library(qqman) # for visualization of data
  library(IlluminaHumanMethylationEPICanno.ilm10b2.hg19) # for annotation for Illumina's EPIC methylation arrays
  library(DMRcate) # for regional analysis
  library(MASS) # for basic statistics
  library(sandwich) # for linear regression (robust sandwich variance estimator)
  library(lmtest) # for testing Linear Regression Models
  library(stringi) # string manipulation
})

##### Cleaning up data #####
# phenotype data
pData(WB.noob2)$Array <- pData(WB.noob2)$Sentrix_Position
pData(WB.noob2)$chip <- pData(WB.noob2)$Sentrix_ID
pheno <- data.frame(pData(WB.noob2))
pheno$race <- factor(pheno$race)
pheno$chip <- factor(pheno$chip)

# chromosome annotations
IlluminaAnnot <- as.data.frame(getAnnotation(WB.noob2))
dim(IlluminaAnnot)

# remove unnecessary data
rm(WB.noob2)

# clean up CpG loci
dim(betas.rcp_norm)
betas.clean <- rmSNPandCH(betas.rcp_norm,  mafcut = 0.05, and = TRUE, rmcrosshyb = TRUE, rmXY= TRUE)
nCpG <- dim(betas.clean)[1]
nCpG

rm(betas.rcp_norm)

# prep
suppressMessages(library(sva))
suppressMessages(library(lumi))


##### Generate Cell Proportions #####
library(EpiDISH)
load(file = "../centEpiFibIC.m.rda")

out.l <- epidish(beta.m = betas.clean, ref.m = centEpiFibIC.m, method = "RPC")
ef <- out.l$estF

png("cell proportions.png", width = 10, height = 10, units = 'in', res =600)
boxplot(ef)
dev.off()

pheno$cell.prob.Epi<-ef[,1]
pheno$cell.prob.Fib<-ef[,2]
pheno$cell.prob.IC<-ef[,3]


##### EWAS without SVA #####

# create covariate matrices for future models
covA <- model.matrix(~c_sex+race+matAge+study+Array+0, data=pheno)
covAP <- model.matrix(~c_sex+race+matAge+study+Array+cell.prob.Epi+0, data=pheno)

predictorList <- c("tri1pm", "tri2pm", "tri12pm")
for (i in c(1,3)) {
  varName <- paste0("results.covA_", predictorList[i])
  assign(varName, cpg.assoc(betas.clean, pheno[,predictorList[i]], covariates = covA, logit.transform = FALSE))
  save(list = varName, file=paste0("results.norm.covA.betavalue_", predictorList[i], ".RData"))
  rm(list = varName)
  
  varName <- paste0("results.covAP_", predictorList[i])
  assign(varName, cpg.assoc(betas.clean, pheno[,predictorList[i]], covariates = covAP, logit.transform = FALSE))
  save(list = varName, file=paste0("results.norm.covAP.betavalue_", predictorList[i], ".RData"))
  rm(list = varName)
}

##### EWAS with SVA #####
modA <- model.matrix(~c_sex+race+matAge+study+Array, data=pheno)
modAP <- model.matrix(~c_sex+race+matAge+study+Array+cell.prob.Epi, data=pheno)

for (i in c(2)) {
  assign(paste0("modA_", predictorList[i]), model.matrix(~c_sex+race+matAge+study+get(predictorList[i])+Array, data=pheno))
  assign(paste0("svaA_", predictorList[i]), sva(betas.clean, get(paste0("modA_", predictorList[i])), modA))
  assign(paste0("covSA_", predictorList[i]), model.matrix(~c_sex+race+matAge+study+get(paste0("svaA_", predictorList[i]))$sv+Array+0, data=pheno))
  varName <- paste0("results.covSA_", predictorList[i])
  assign(varName, cpg.assoc(betas.clean, pheno[,predictorList[i]], covariates = get(paste0("covSA_", predictorList[i])), logit.transform = FALSE))
  save(list = varName, file=paste0("results.norm.covSA.betavalue_", predictorList[i], ".RData"))
  rm(list = varName)
  
  assign(paste0("modAP_", predictorList[i]), model.matrix(~c_sex+race+matAge+study+get(predictorList[i])+Array+cell.prob.Epi, data=pheno))
  assign(paste0("svaAP_", predictorList[i]), sva(betas.clean, get(paste0("modAP_", predictorList[i])), modAP))
  assign(paste0("covSAP_", predictorList[i]), model.matrix(~c_sex+race+matAge+study+get(paste0("svaAP_", predictorList[i]))$sv+Array+cell.prob.Epi+0, data=pheno))
  varName <- paste0("results.covSAP_", predictorList[i])
  assign(varName, cpg.assoc(betas.clean, pheno[,predictorList[i]], covariates = get(paste0("covSAP_", predictorList[i])), logit.transform = FALSE))
  save(list = varName, file=paste0("results.norm.covSAP.betavalue_", predictorList[i], ".RData"))
  rm(list = varName)
}

