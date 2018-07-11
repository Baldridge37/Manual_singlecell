# Manual PCA pseudotime analysis

```
setwd("C:/Documents/normalise_workflow/")

install.packages("devtools") # If not already installed
install.packages("C:/Documents/pkgmaker_0.22.zip", lib="C:/Rpackages") 
install.packages("pkgmaker")


#Will probably need to install stringi, data.table and httpuv
library(pkgmaker)
library(devtools)
library(ggplot2)
library(reshape2)
library(scater)
library(SingleCellExperiment)
#install_github("kieranrcampbell/phenopath")


#Phenodata includes; sample name, cell number, plate, well,	sample_number,	sample,	sample_id,
#Cell,	Batch,	read_count,	feature_number,	Cell_type
phn<- read.table("Combined_annotation.txt", header=TRUE, sep="\t")
counts<-read.table("new_counts.txt", header = TRUE)
tpms<-read.table("protoplast_tpm.txt", header = TRUE)

#May not need feature data
feat<-read.table("C:/Documents/SingleCell/gene_feature.txt", header=TRUE, sep="\t")
fd <- new("AnnotatedDataFrame", data = feat)

#Create Single cell experiment with counts data then add tpms.
sce <- SingleCellExperiment(assays = list(counts = as.matrix(counts)), colData = phn)

tpm(sce) <- tpms
exprs(sce) <- log2(tpm(sce) + 1)

sce <- calculateQCMetrics(sce)

plotPhenoData(sce, aes(x = total_counts, y = total_features))

#### Extra QC
total_C <- sce@colData$total_counts
feature_c <- sce@colData$total_features
mean(total_C)
mean(feature_c)

###Filter cells to keep
sce$to_keep <- sce$total_counts > 2.5e5 & sce$total_features > 1.5e3 & sce$cell_number == 1
PhenoPlot<-plotPhenoData(sce, aes(x = total_counts, y = total_features, colour = to_keep)) +
  labs(title ="Read counts and genes detected for all sequenced wells" ,
       subtitle = "QC pass: total_features > 1500 and total_counts > 25000",
       x="Read counts", y="Genes detected")+
  scale_colour_discrete(name = "Passed\nfiltering")

plotPCA(sce, colour_by=total_features)

ggsave(filename="Readvsfeat_scatter.tiff", plot=PhenoPlot, dpi=600, height=12, width=15, units="cm")

sce_qc <- sce[, sce$to_keep]
keep_feature <- rowSums(exprs(sce_qc) > 0) > 0
sce_qc <- sce_qc[keep_feature,]

plotQC(sce_qc, type = 'find', var = 'total_features', ntop = 2e3)

###normalise for total features. Could also do read count
m <- model.matrix(~ sce_qc$total_features)
sce_qc <- normaliseExprs(sce_qc, design = m)
exprs(sce_qc) <- norm_exprs(sce_qc)

sce <- sce_qc
rm(sce_qc)
print(sce)

plotPhenoData(sce, aes(x = total_counts, y = total_features))

plotPCA(sce, ncomponents=4, colour_by = "AT5G07320") #A9

####Now cluster and analyse with SC3
biocLite("SC3")
library(SC3)

rowData(sce)$feature_symbol <- rownames(sce)
sce <- sc3(sce, ks = 2:5, biology = TRUE, n_cores = 1)
sc3_interactive(sce)

#2clusters best

plotPCA( sce, colour_by = "sc3_3_clusters", size_by = "sc3_3_log2_outlier_score")
plotPCA( sce, colour_by = "sc3_2_clusters", size_by = "sc3_2_log2_outlier_score")
plotPCA( sce, colour_by = "AT1G02050") #LAP6 Cluster2

plotPCA(sce, colour_by = "Batch")

#### Extract cluster information

col_data <- colData(sce)
head(col_data[ , grep("sc3_", colnames(col_data))])
head(col_data)

Clusters<-row_data[ ,colnames(row_data)]

write.csv(Clusters, file="Gene_cluster_markers_all.csv")


###### EXPORTING PCA components
library(magrittr)
assay(sce) %>% log1p %>% rowVars -> vars
names(vars) <- rownames(sce)
vars <- sort(vars, decreasing = TRUE)

sce_sub <- sce[names(vars[1:100]),]
sce_sub

library(Rtsne)
set.seed(5252)

pca_data <- prcomp(t(log1p(assay(sce_sub))))
tsne_data <- Rtsne(pca_data$x[,1:50], pca = FALSE)

reducedDims(sce_sub) <- SimpleList(PCA=pca_data$x, TSNE=tsne_data$Y)
sce_sub

reducedDims(sce_sub)
reducedDimNames(sce_sub)
head(reducedDim(sce_sub, "PCA")[,1:2])

#get PCA results
PCA_results<-reducedDim(sce_sub, "PCA")[,1:2]
PCA_results<-as.matrix(PCA_results)
#PCA_results<-sort.list(PCA_results, decreasing=F)

PCA_R<-t(PCA_results)

#Get normalised log(tpm) data
norm_tpm<-exprs(sce)
norm_tpms<-rbind(norm_tpm, t(PCA_results))
rm(norm_tpm)

#For some reason the PCA loadings are the negative of what they are on the plot

NT<-t(norm_tpms)
NT<-as.data.frame(NT)

NT$PC1<--NT$PC1
NT$PC2<--NT$PC2



#Get non normalised TPMs
tpm_e<-tpm(sce)
tpms_e<-rbind(tpm_e, t(PCA_results))
rm(tpm_e)

NT2<-t(tpms_e)
NT2<-as.data.frame(NT2)

NT2$PC1<--NT2$PC1
NT2$PC2<--NT2$PC2


write.csv(file = "Normalised_tpm_princomps.csv", NT)
write.csv(file = "Raw_tpm_princomps.csv", NT2)

#PC2 correlates with late gene expression
plot(NT$PC2, NT$AT1G01080)

#Melt data for ggplot
NT_melt<-melt(NT, id.vars = c("PC1", "PC2"), measure.vars = )
source("C:/Documents/SLM_wba_plotting/master_theme_ggplot.R")

gene<-"AT2G47400"

NT_gene<-NT_melt[NT_melt$variable==gene,]

ps<-ggplot(NT_gene, aes(PC2, as.numeric(value)))
ps+geom_smooth(color="black", method="loess", size=2)+geom_point(color="red", size =1.5)+
  labs(title = "", y = "Normalised TPM expression", x = "Pseudotime(PC2)") + 
  theme_bw()


####### Plotting
sc3_plot_silhouette.SingleCellExperiment(sce, 2)

```
