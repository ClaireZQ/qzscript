
#########SC3#####################
library(SingleCellExperiment)
library(SC3)
library(scater)

test <- melanoma[,-1]
yan <- test[,which(as.numeric(test[3,]) !=0)]
yan <- yan[4:23687, ]
yan <- yan[,-1]

#log2(yan/10 +1)
ann <- melanoma[1, match(colnames(yan), colnames(melanoma))]
rownames(ann) <- NULL
ann <- t(ann)

ann[,1] <- paste0("Mel",ann[,1],sep="")
ann <- as.data.frame(ann)
colnames(ann) <- "patient" 


sce <- SingleCellExperiment(
    assays = list(
        counts = as.matrix((2^yan - 1)*10) ,
        logcounts = as.matrix(yan)
    ),
    colData = ann
)


plotPCA(sce, colour_by = "patient")

##run SC3
sce <- sc3(sce, ks = 2:10, biology = TRUE)

##
pdf("consensus_6.pdf")
sc3_plot_consensus(sce, k = 6)
dev.off()

# sce <- sc3_estimate_k(sce)
# metadata(sce)$sc3$k_estimation

svm_labels <- colData(sce)$sc3_10_clusters
nonmalig <- melanoma[,match(colnames(yan),colnames(melanoma))]
test <- as.factor(as.numeric(nonmalig[3,]))

ARI_SC3 <- adjustedRandIndex(test, svm_labels)
# > ARI_SC3                                                                                                                                                                    
# [1] 0.3362217

save(sce, svm_labels, test, nonmalig, ARI_SC3, file="sc3.RData")



sample_10 <- metadata(sce)$sc3$consensus$`6`$ `consensus`
rownames(sample_10) <- colnames(yan)
colnames(sample_10) <- colnames(yan)
tsne_out <- Rtsne(sample_10, is.distance=TRUE, theta=0.1, perplexity=100, verbose=TRUE, pca=F, check_duplicates = F)

tsne_out_1 <- data.frame(tsne_out$Y,colData(sce)$sc3_6_clusters)
names(tsne_out_1) <- c("V1","V2", "cluster")


pdf("t-SNE_6_100_2_2.pdf")
ggplot(tsne_out_1, aes(V1,V2, color=cluster)) + geom_point(size=0.25) +
	theme_bw() + xlab("t-SNE1") +ylab("t-SNE2") +
	theme(panel.grid =element_blank(), legend.title=element_blank()) +
	guides(colour = guide_legend(override.aes = list(size=5)))
dev.off()

pdf("t-SNE_6_100_2.pdf")
ggplot(tsne_out_1, aes(V1,V2, color=cluster)) + geom_point(size=0.25) +
	theme_bw() + xlab("t-SNE1") +ylab("t-SNE2") +
	theme(panel.grid =element_blank()) + guides(color=guide_legend(title="SC3", override.aes = list(size=5))) 
dev.off()

#############TSNE+KMEANS##############

for(i in seq(5,30, 5)){
	print(i)
sce <- plotTSNE(sce, rand_seed = 1, return_SCE = TRUE, perplexity = i)
colData(sce)$tSNE_kmeans <- as.character(kmeans(sce@reducedDims$TSNE, centers = 5)$clust)

pdf(paste0("plotTSNE",i, ".pdf", sep=""))
plotTSNE(sce, rand_seed = 1, colour_by = "tSNE_kmeans")
dev.off()

}




sce <- plotTSNE(sce, rand_seed = 1, return_SCE = TRUE, perplexity = 10)
colData(sce)$tSNE_kmeans <- as.character(kmeans(sce@reducedDims$TSNE, centers = 6)$clust)

pdf("plotTSNE_6_new.pdf")
plotTSNE(sce, rand_seed = 1, colour_by = "tSNE_kmeans")
dev.off()

pdf("sc3_plot_silhouette_6.pdf")
sc3_plot_silhouette(sce, k=6)
dev.off()

nonmalig <- melanoma[,match(colnames(yan),colnames(melanoma))]
test <- as.factor(as.numeric(nonmalig[3,]))

test2 <- as.numeric(colData(sce)$tSNE_kmeans)
test2  <- as.factor(test2)

ARI_tsne_kmeans <- adjustedRandIndex(test, test2)
# > ARI_tsne_kmeans                                                                                                                                                            
# [1] 0.3054177  
save(sce, test,test2, ARI_tsne_kmeans, nonmalig, file="tsne_kmeans.RData")



########pcaReduce##############

library(pcaReduce)
library(SC3)
library(scater)
library(SingleCellExperiment)
library(pheatmap)
library(mclust)
set.seed(1234567)
install.packages("~/software/pcaReduce-master",repos = NULL, type="source", dependencies = TRUE,lib='/home/liukuai/R/x86_64-pc-linux-gnu-library/3.4') 

input <- nonmalig[match(hv.genes, rownames(nonmalig)),]
#input <- logcounts(sce[rowData(sce)$sc3_gene_filter, ])
pca.red <- PCAreduce(t(input), nbt = 1, q = 30, method = 'S')[[1]]
colData(sce)$pcaReduce <- as.character(pca.red[,32 - 6])
plotPCA(sce, colour_by = "pcaReduce", legend="none")


adjustedRandIndex(test, colData(sce)$pcaReduce)
# > adjustedRandIndex(test, colData(sce)$pcaReduce)
# [1] 0.233917 




#######seurat##############
load("ARI.RData")

library(Seurat)
library(dplyr)
library(Matrix)

# Examine the memory savings between regular and sparse matrices

embryo <- CreateSeuratObject(raw.data = as.matrix((2^yan - 1)*10), min.cells = 3, min.genes = 200, 
    project = "10X_embryo")

embryo <- NormalizeData(object = embryo, normalization.method = "LogNormalize", 
    scale.factor = 10000)

pdf("FindVariableGenes.pdf")
embryo <- FindVariableGenes(object = embryo, mean.function = ExpMean, dispersion.function = LogVMR, 
    x.low.cutoff = 0.0125, x.high.cutoff = 3, y.cutoff = 0.5, y.high.cutoff = 6)
dev.off()


embryo <- FindVariableGenes(object = embryo, mean.function = ExpMean, dispersion.function = LogVMR, 
    do.plot = FALSE)
hv.genes <- head(rownames(embryo@hvg.info), 4000) ###

#embryo2 <- ScaleData(object = embryo, genes.use = embryo@var.genes, display.progress = FALSE)
embryo2 <- ScaleData(object = embryo, genes.use = hv.genes,display.progress = FALSE)
###############runPCA#########################


#embryo2 <- RunPCA(object = embryo2,  pcs.compute = 30, genes.use = embryo2@var.genes,  do.print = TRUE, pcs.print = 1:5, genes.print = 5)
embryo2 <- RunPCA(object = embryo2,  pcs.compute = 30, genes.use = hv.genes,  do.print = TRUE, pcs.print = 1:5, genes.print = 5)

adjustedRandIndex(test, seurat_cluster)

#############select pca#######################

embryo2 <- JackStraw(object = embryo2, num.replicate = 100, num.pc=30)


PCElbowPlot(object = embryo2)



######find cluster################3

for(i in c(5:8)){
	embryo2 <- FindClusters(object = embryo2, reduction.type = "pca", dims.use = 1:i, 
    	resolution = 0.4, print.output = 0, save.SNN = TRUE)
	embryo2 <- RunTSNE(object = embryo2, dims.use = 1:i, do.fast = TRUE)

	pdf(paste("tSNE_", i, ".pdf", sep = ""))
	TSNEPlot(object = embryo2, no.legend=F, do.lable=T, pt.size=0.5)
	dev.off()

	embryo2.markers <- FindAllMarkers(object = embryo2, only.pos = TRUE, min.pct = 0.25,
       thresh.use = 0.25)
	write.csv(embryo2.markers, file = paste(i,"_PCA_marker.gene.csv"))
	}

# > adjustedRandIndex(test,seurat_cluster)                                                                                                                                     
# [1] 0.8366344

#########find subsets###############################

embryo2 <- FindClusters(object = embryo2, reduction.type = "pca", dims.use = 1:20, resolution = 0.6, print.output = 0, save.SNN = TRUE)
embryo2 <- RunTSNE(object = embryo2, dims.use = 1:20, do.fast = TRUE)

pdf("20_PCA23.pdf", width=8)
TSNEPlot(object = embryo2, do.return = TRUE, no.legend = TRUE, do.label = TRUE, pt.size=0.5)
dev.off()

embryo2.markers <- FindAllMarkers(object = embryo2, only.pos = TRUE, min.pct = 0.25, thresh.use = 0.25)
write.csv(embryo2.markers, file = "10_PCA_marker.gene.csv")







###########SINCERA###############
#input <- logcounts(sce[rowData(sce)$sc3_gene_filter, ])
dat <- apply(input, 1, function(y) scRNA.seq.funcs::z.transform.helper(y))
dd <- as.dist((1 - cor(t(dat), method = "pearson"))/2)
hc <- hclust(dd, method = "average")
num.singleton <- 0
kk <- 1
for (i in 2:dim(dat)[2]) {
    clusters <- cutree(hc, k = i)
    clustersizes <- as.data.frame(table(clusters))
    singleton.clusters <- which(clustersizes$Freq < 2)
    if (length(singleton.clusters) <= num.singleton) {
        kk <- i
    } else {
        break;
    }
}
cat(kk)

pdf("pheatmap.pdf")
pheatmap(
    t(dat),
    cluster_cols = hc,
    cutree_cols = kk,
    kmeans_k = 100,
    show_rownames = FALSE
)
dev.off()

colData(sce)$SINCERA <- as.character(cutree(hc, k = kk))
adjustedRandIndex(test, colData(sce)$SINCERA)
# > adjustedRandIndex(test, colData(sce)$SINCERA)                                                                                                                              
# [1] 0.09310616


########SNN-Cliq###############
distan <- "euclidean"
par.k <- 3
par.r <- 0.7
par.m <- 0.5
# construct a graph
scRNA.seq.funcs::SNN(
    data = t(input),
    outfile = "snn-cliq.txt",
    k = par.k,
    distance = distan
)
# find clusters in the graph
snn.res <- 
    system(
        paste0(
            "python utils/Cliq.py ", 
            "-i snn-cliq.txt ",
            "-o res-snn-cliq.txt ",
            "-r ", par.r,
            " -m ", par.m
        ),
        intern = TRUE
    )
cat(paste(snn.res, collapse = "\n"))
snn.res <- read.table("res-snn-cliq.txt")
# remove files that were created during the analysis
system("rm snn-cliq.txt res-snn-cliq.txt")

colData(deng)$SNNCliq <- as.character(snn.res[,1])
plotPCA(deng, colour_by = "SNNCliq")

adjustedRandIndex(colData(deng)$cell_type2, colData(deng)$SNNCliq)



sc3:0.3362217 
TSNE+KMEANS:0.3054177  
pcaReduce:0.233917 
seurat:0.8366344
SINCERA:0.09310616
#SNN-Cliq:

setwd("~/Downloads/")
ari_score <- read.delim("ARI_score.txt", stringsAsFactors = F, header = F)


pdf("ARI2.pdf", width = 10, height = 10)
ggplot(ari_score, aes(V1,V2, fill=V1)) + geom_bar(stat="identity", width = 0.8) + 
  theme_bw() + xlab("") +ylab("ARI\n") + scale_y_continuous(limits=c(0,1))+
  theme(legend.position = "none",
        axis.text.y = element_text(size = 16),
        panel.border = element_blank(),
        axis.title.y = element_text(size = 22, color = "black", face = "bold"),
        axis.line = element_line(size=0.8, colour = "black"),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        panel.grid =element_blank())

dev.off()

pdf("ARI2.pdf", width = 10, height = 10)
ggplot(ari_score, aes(V1,V2, fill=V1)) + geom_bar(stat="identity", width = 0.6) + 
  theme_bw() + xlab("") +ylab("ARI\n") + scale_y_continuous(limits=c(0,1))+ scale_fill_hue(c=70, l=80) +
  theme(legend.position = "none",
        axis.text.y = element_text(size = 16),
        panel.border = element_blank(),
        axis.title.y = element_text(size = 22, color = "black", face = "bold"),
        axis.line = element_line(size=0.8, colour = "black"),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        panel.grid =element_blank())

dev.off()