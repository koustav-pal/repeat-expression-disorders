
#################################
# setting the working directory #
#################################

setwd("~/Google Drive/EH") # the path will need to be changed!
#setwd("C:/Users/IDDO/Google Drive/EH")
myfile <- "All_samples_enrolment.txt"

#############################################
# reading the text file into an R dataframe #
#############################################

raw_counts <- read.table(myfile, header=TRUE, 
                         row.names=1, skip=5) 

metadata <- read.table(myfile, header=FALSE,
                       nrows=5, row.names=1, stringsAsFactors=FALSE)

colnames(metadata) <- colnames(raw_counts)

###################################
# filter out low coverage samples #
###################################

filt <- colSums(raw_counts) >= 40000 # in this case i filtered out samples below 100000 miRNA counts
raw_counts <- raw_counts[,filt]
metadata <- metadata[,filt]

#################################
# filter out low coverage miRNA #
#################################

# option 2: keep only miRNA with at least 5 counts (sum across all samples)

keep <- (rowSums(raw_counts)>=25300)>0
filtered_counts2 <- raw_counts[keep,]

########################################
# define metadata variables as factors #
########################################


batch <- metadata[2,]
batch <- factor(batch, levels=unique(as.character(batch))) 
sex <- metadata[3,]
sex <- factor(sex, levels=c("M","F"))
Onset <- metadata[4,]
Onset <- factor(Onset, levels=c("Bulbar","Non_bulbar"))
Riluzole <- metadata [5,]
Riluzole <- factor (Riluzole, levels=c("Yes","No"))

plot(batch, sex, xlab="batch")
plot(batch, log10(colSums(raw_counts)), xlab="batch",
     ylab="log10 total miRNA count", ylim=c(3,7), col=rainbow(20))
abline(a=4, b=0)
plot(Age, log10(colSums(raw_counts)), xlab="Age",
     ylab="log10 total miRNA count", ylim=c(3,7), notch=FALSE)

#################################################
# define DESeq parameters countData and colData #
#################################################

library (DESeq2)
countData <- as.matrix(round(filtered_counts2)) 
colData <- data.frame (batch, sex, Onset, Riluzole)

#############################################
# build the DESeq dataset with model design #
#############################################

dds <- DESeqDataSetFromMatrix(countData = countData,
                              colData = colData,
                              design = ~ batch + sex + Onset + Riluzole) 

####################
# PCA using DESeq2 #
####################

vst <- varianceStabilizingTransformation(dds)
vst.matrix <- assay(vst) 

# PCA plots
plotPCA(vst, intgroup=c("batch"))

data <- plotPCA(vst, intgroup=c("Age"), returnData=TRUE)
percentVar <- round(100 * attr(data, "percentVar"))
data$sums <- colSums(revertdata)

library(ggplot2)
pd <- position_dodge(.5) 
ggplot(data, aes(PC1, PC2, color=Batch, shape=Rx)) +
  geom_point() +
  #  geom_text(aes(PC1-2.5, PC2-0), size=3) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance"))

###################################################
# ComBat - correcting the counts for batch effect #
####################################################

library(sva)
library(edgeR)

cpmdata = cpm(filtered_counts2, log=TRUE, normalized.lib.sizes=TRUE, prior.count=0.25)

modcombat = model.matrix(~ sex + Onset + Riluzole, data=data.frame(t(cpmdata)))

combat_edata = ComBat(dat=cpmdata, batch=batch, par.prior=FALSE, prior.plots=TRUE)

revertdata = t(((2^t(combat_edata))*colSums(filtered_counts2)/1000000)-0.25)
revertdata[revertdata<0] <- 0

par(pty="s")
plot (as.matrix(filtered_counts2), revertdata, log="xy",
      xlab="raw_counts", ylab="batch-corrected counts", col="green", cex=.5)
abline(a=0,b=1, col="blue")

write.csv(revertdata, file="Batch_corrected_counts.csv")

