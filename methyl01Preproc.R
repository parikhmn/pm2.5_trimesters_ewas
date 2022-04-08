# Methylation Data Preprocessing QC
# Milan Parikh

knitr::opts_knit$set(root.dir = "../")

# import libraries
suppressPackageStartupMessages({
  library(minfi)
  library(shinyMethyl)
  library(pryr)
  library(matrixStats)
  library(limma)
  library(reshape,scales)
  require(sva)
  library(IlluminaHumanMethylationEPICmanifest)
  library(IlluminaHumanMethylationEPICanno.ilm10b2.hg19)
})

# path to methylation idat files
idatPath <- paste0("/path/to/idats/")

# import phenotypic data with all PM2.5 and covariate data listed for each sample ID
targets.pheno <- read.csv("/path/to/sample_sheet.csv", strip.white=T, stringsAsFactors=F)

# read in IDAT images for sample IDs in targets.pheno
targets.pheno$Basename <- paste0(targets.pheno$Sentrix_ID,"_",targets.pheno$Sentrix_Position)
WB.i <- read.metharray.exp(base=idatPath, targets=targets.pheno, verbose=T, force=T, extended=T)

# change working directory to simplify output
setwd("/path/to/outputFolder1/")


####### Sex analysis #######
# check on sex matching
Gbeta <- mapToGenome(WB.i) # convert to MethylSet
sexPreds1 <- getSex(Gbeta)
mismatchTable1 <- table(pData(WB.i)$c_sex, sexPreds1$predictedSex, useNA='ifany')
mismatchFreq <- as.data.frame(mismatchTable1)
mismatchFreqSorted <- sort_df(mismatchFreq, vars = "Freq")
compSexTable <- data.frame(cbind(pData(WB.i)$c_sex, sexPreds1$predictedSex))
colnames(compSexTable) <- c("Actual", "Predicted")
compSexTable$Mismatch <- apply(compSexTable, 1,
                               function(rowCells) (rowCells[1] == mismatchFreqSorted[1,1] & rowCells[2] == mismatchFreqSorted[1,2]) | (rowCells[1] == mismatchFreqSorted[2,1] & rowCells[2] == mismatchFreqSorted[2,2]))
WB.ii<-WB.i[,!compSexTable$Mismatch] # remove mismatched subjects

# visualize sex intensities
Gbeta = addSex(Gbeta)
png("male_female_XYintensity.png", width = 10, height = 10, units = 'in', res =600)  #often 8X8 or 12*8 look good also
plotSex(Gbeta, id=pData(WB.i)$SEX)
dev.off()

# recheck sex matching
Gbeta2 <- mapToGenome(WB.ii)
sexPreds2 <- getSex(Gbeta2)
table(pData(WB.ii)$c_sex, sexPreds2$predictedSex, useNA='ifany')

save(Gbeta2, file="Gbeta2.RData")

####### Bead counting #######
suppressMessages(library(wateRmelon))
bc <- beadcount(WB.ii) # count number of beads identifying each probe (remove < 3)

beadc <- function(x){
  sum(is.na(x))
}

bab <- apply(bc, 1, beadc)
badbead <- which(bab > (ncol(bc)*0.05)) # exclude if larger than 5% of samples
length(badbead)
WB <- WB.ii[-badbead,] # removing bad beads

dim(WB.ii)
dim(WB)

# Clean up workspace
rm(WB.i, Gbeta2, Gbeta, sexPreds1, sexPreds2, mismatchFreq, mismatchFreqSorted, compSexTable); gc()


####### Normalization #######
WB.noob <- preprocessNoob(WB)

minfiQC(WB.noob, fixOutliers = TRUE, verbose = FALSE)
qc <- getQC(WB.noob)
head(qc)
WB.noob <- addQC(WB.noob, qc=qc) # store qc measures back in dataset

png("Plot3.png", width = 10, height = 10, units = 'in', res =600)
plotQC(qc, badSampleCutoff = 10.5)
dev.off()

# Distribution of beta-values: before and after normalization
png("Plot4.png", width = 10, height = 10, units = 'in', res =300)
densityPlot(WB, main = "density plots before and after preprocessing", pal="#440154FF", ylim=c(0,4.5))
densityPlot(WB.noob, add = F, pal = "#FDE725FF")
legend("topleft", c("Noob","Raw"), 
       lty=c(1,1), title="Normalization", 
       bty='n', cex=1.3, col=c("#FDE725FF","#440154FF"))
dev.off()


####### Detection P Values #######
# finds P values for difference between signal and noise
library(ewastools)
detP_oob <- ewastools::detectionP.minfi(WB)
dim(detP_oob)

# select probes whose signal was not significantly different from noise and remove
detP_oob[detP_oob > 0.01] <- NA
detP_oob <- na.omit(detP_oob)
dim(detP_oob)

# filter out bad p-value probes from dataset
WB.noob2 <- WB.noob[rownames(getAnnotation(WB.noob)) %in% rownames(detP_oob),]
dim(WB.noob2)


####### RCP Probe Type Adjustment ########
library(ENmix)
betas.rcp_norm <- rcp(WB.noob2)
dim(betas.rcp_norm)


####### Save data for next step #######
save(betas.rcp_norm,file="betas.rcp_norm.RData") # normalized data to use in methyl02
save(WB.noob2, file="WB.noob2.RData") # unadjusted data

