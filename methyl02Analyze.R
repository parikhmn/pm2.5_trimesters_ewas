# Methylation Data Initial Analysis
# Milan Parikh

knitr::opts_knit$set(root.dir = "../")
setwd("/path/to/outputFolder1/")

##### Import Libraries and Data #####
suppressMessages(library(minfi))

load("WB.noob2.RData")
load("betas.rcp_norm.RData")

setwd("/path/to/outputFolder2/")

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

# convert beta values to M values
m.clean <- beta2m(betas.clean)


##### Generate Cell Proportions #####
library(EpiDISH)
load(file = "../centEpiFibIC.m.rda")

out.l <- epidish(beta.m = betas.clean, ref.m = centEpiFibIC.m, method = "RPC")
ef <- out.l$estF

png("cell proportions.png", width = 10, height = 10, units = 'in', res = 600)
boxplot(ef)
dev.off()

pheno$cell.prob.Epi<-ef[,1]
pheno$cell.prob.Fib<-ef[,2]
pheno$cell.prob.IC<-ef[,3]


##### EWAS without SVA #####

# create covariate matrices for future models
covA <- model.matrix(~c_sex+race+matAge+study+Array+0, data=pheno)
covAP <- model.matrix(~c_sex+race+matAge+study+Array+cell.prob.Epi+0, data=pheno)

predictorList1 <- c("tri1pm", "tri12pm")
for (i in 1:length(predictorList)) {
  varName <- paste0("results.covAP_", predictorList[i])
  assign(varName, cpg.assoc(betas.clean, pheno[,predictorList[i]], covariates = covAP, logit.transform = TRUE))
  png(paste0(predictorList[i], "_c_sex+race+matAge+study+Array+cell.prob.Epi_normDNAm_Mvalue.png"), width = 10, height = 10, units = 'in', res =600)
  plot(get(varName), main=paste0("QQ plot for association between methylation(Mvalue) and ", predictorList[i], "\nc_sex+race+matAge+study+Array+cell.prob.Epi"), gcdisplay=TRUE)
  dev.off()
  save(list = varName, file=paste0("results.norm.covAP.mvalue_", predictorList[i], ".RData"))
  rm(list = varName)
}

##### EWAS with SVA #####
modA <- model.matrix(~c_sex+race+matAge+study+Array, data=pheno)
modAP <- model.matrix(~c_sex+race+matAge+study+Array+cell.prob.Epi, data=pheno)

predictorList2 <- c("tri2pm")
for (i in 1:length(predictorList2)) {
  assign(paste0("modAP_", predictorList2[i]), model.matrix(~c_sex+race+matAge+study+get(predictorList2[i])+Array+cell.prob.Epi, data=pheno))
  assign(paste0("svaAP_", predictorList2[i]), sva(m.clean, get(paste0("modAP_", predictorList2[i])), modAP))
  assign(paste0("covSAP_", predictorList2[i]), model.matrix(~c_sex+race+matAge+study+get(paste0("svaAP_", predictorList2[i]))$sv+Array+cell.prob.Epi+0, data=pheno))
  varName <- paste0("results.covSAP_", predictorList2[i])
  assign(varName, cpg.assoc(betas.clean, pheno[,predictorList2[i]], covariates = get(paste0("covSAP_", predictorList2[i])), logit.transform = TRUE))
  png(paste0(predictorList2[i], "_c_sex+race+matAge+study+sva+Array+cell.prob.Epi_normDNAm_Mvalue.png"), width = 10, height = 10, units = 'in', res =600)
  plot(get(varName), main=paste0("QQ plot for association between methylation(Mvalue) and ", predictorList2[i], "\nc_sex+race+matAge+study+sva+Array+cell.prob.Epi"), gcdisplay=TRUE)
  dev.off()
  save(list = varName, file=paste0("results.norm.covSAP.mvalue_", predictorList2[i], ".RData"))
  rm(list = varName)
}


##### Plot Stuff #####
file_names <- as.list(dir(pattern="*.RData"))
lapply(file_names, load, .GlobalEnv)

selectedData <- list(tri1pm_c_sex_race_matAge_study_Array_cell.prob.Epi_normDNAm_Mvalue = results.covAP_tri1pm,
                     tri2pm_c_sex_race_matAge_study_sva_Array_cell.prob.Epi_normDNAm_Mvalue = results.covSAP_tri2pm,
                     tri12pm_c_sex_race_matAge_study_Array_cell.prob.Epi_normDNAm_Mvalue = results.covAP_tri12pm)

figureTitles <- c("Manhattan Plot: Associations of DNAm vs\nTri 1 PM2.5 + Child Sex + Mat Race + Mat Age + Study Batch + Sample Row + Cell Comp",
                  "Manhattan Plot: Associations of DNAm vs\nTri 2 PM2.5 + Child Sex + Mat Race + Mat Age + Study Batch + Sample Row + Cell Comp + SVA",
                  "Manhattan Plot: Associations of DNAm vs\nCombined PM2.5 + Child Sex + Mat Race + Mat Age + Study Batch + Sample Row + Cell Comp")

plotFigures <- function(dataset, name, title) {
  
  IlluminaAnnot.use<-data.frame(cbind(IlluminaAnnot$Name, IlluminaAnnot$chr, IlluminaAnnot$pos,IlluminaAnnot$UCSC_RefGene_Name))
  colnames(IlluminaAnnot.use)<-c('CPG.Labels', 'chr', 'pos', 'UCSC_RefGene_Name')
  
  dataset$coefficients$CPG.Labels<-rownames(dataset$coefficients)
  
  
  dataset.use<-merge(merge(dataset$results, IlluminaAnnot.use, by="CPG.Labels", all.x=TRUE),
                     dataset$coefficients, by="CPG.Labels")
  
  
  datamanhat <- data.frame(CpG=dataset.use$CPG.Labels,Chr=as.character(dataset.use$chr),
                           Mapinfo=dataset.use$pos, UCSC_RefGene_Name=dataset.use$UCSC_RefGene_Name,
                           Pval=dataset.use$P.value, Eff.Size = dataset.use$effect.size, Std.Error = dataset.use$std.error,
                           FDR= dataset.use$FDR, gc.p.value=dataset.use$gc.p.value, Holm.sig=dataset.use$Holm.sig)
  
  write.csv(datamanhat, paste(name, '.csv', sep=""))
  
  #'## Manhattan plot for cell-type adjusted EWAS
  #' Reformat the variable Chr (so we can simplify and use a numeric x-axis)
  datamanhat$Chr <- as.numeric(sub("chr","",datamanhat$Chr))
  datamanhat$Mapinfo <- as.numeric(datamanhat$Mapinfo)
  
  datamanhat<-datamanhat[!is.na(datamanhat$Chr),]
  
  plot_title2 <- paste0("Manhattan Plot ", name)
  #' the function manhattan needs data.frame including CpG, Chr, MapInfo and Pvalues
  png(paste0(plot_title2, ".png"), res = 1200, width = 10, height = 10, units = "in")  #often 8X8 or 12*8 look good also
  manhattan(datamanhat, "Chr", "Mapinfo", "Pval", "CpG",
            genomewideline = -log10(0.05/(nCpG)), suggestiveline = FALSE,
            main = title, ylim=c(0,8))
  dev.off()
}

mapply(plotFigures, dataset = selectedData, name = names(selectedData), title = figureTitles)
