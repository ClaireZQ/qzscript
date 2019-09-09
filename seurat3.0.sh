
touch all.csv
echo 'library_id,molecule_h5'>>all.csv
echo 'ems,/home/renjun/project/op/result2/ems/outs/molecule_info.h5'>>all.csv
echo 'normal,/home/renjun/project/op/result2/normal/outs/molecule_info.h5'>>all.csv

cellranger aggr --id=aggop \
                --csv=all.csv \
                --normalize=mapped

library(dplyr)
library(Seurat)

# Load the PBMC dataset
op.data <- Read10X(data.dir = "/home/renjun/project/op/result2/aggop/outs/filtered_feature_bc_matrix/")
# Initialize the Seurat object with the raw (non-normalized data).
op <- CreateSeuratObject(counts = op.data, project = "aggop", min.cells = 3, min.features = 200)

op[["percent.mt"]] <- PercentageFeatureSet(object = op, pattern = "^MT-")

pdf('vlnplot.pdf', width = 10, height = 6)
VlnPlot(object = op, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"))
dev.off()

plot1 <- FeatureScatter(object = op, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(object = op, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")

pdf('nCount_RNA.pdf', width = 10)
CombinePlots(plots = list(plot1, plot2))
dev.off()

op <- subset(x = op, subset = nFeature_RNA > 500 & nFeature_RNA < 3000 & percent.mt < 10)

#Normalizing the data
op <- NormalizeData(object = op, normalization.method = "LogNormalize", scale.factor = 10000)
op <- NormalizeData(object = op)


#Identification of highly variable features (feature selection)
op <- FindVariableFeatures(object = op, selection.method = "vst", nfeatures = 3000)

# Identify the 10 most highly variable genes
top10 <- head(x = VariableFeatures(object = op), 10)

# plot variable features with and without labels
plot1 <- VariableFeaturePlot(object = op)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)

pdf("VariableFeatures.pdf", width = 15)
CombinePlots(plots = list(plot1, plot2))
dev.off()

#scaling the data
all.genes <- rownames(x = op)
op <- ScaleData(object = op, features = all.genes)
op <- ScaleData(object = op, vars.to.regress = "percent.mt")

######pca
op <- RunPCA(object = op, features = VariableFeatures(object = op))
print(x = op[["pca"]], dims = 1:5, nfeatures = 5)

pdf('vizdimpca.pdf')
VizDimLoadings(object = op, dims = 1:2, reduction = "pca")
dev.off()

pdf('pca.pdf')
DimPlot(object = op, reduction = "pca")
dev.off()

pdf('pcaheatmap.pdf', height = 25)
DimHeatmap(object = op, dims = 1:15, cells = 500, balanced = TRUE)
dev.off()

#####tsne
op <- RunTSNE(object = op, features = VariableFeatures(object = op), dims = 1:11)
pdf('tsne.pdf')
DimPlot(object = op, reduction = "tsne")
dev.off()


op <- JackStraw(object = op, num.replicate = 100)
op <- ScoreJackStraw(object = op, dims = 1:20)

pdf('JackStrawPlot.pdf')
JackStrawPlot(object = op, dims = 1:15)
ElbowPlot(object = op)
dev.off()


op <- FindNeighbors(object = op, features = VariableFeatures(object = op), reduction = 'pca', dims = 1:8)
op <- FindClusters(object = op, resolution = 0.5)



op@meta.data[,'orig.ident']=as.character(op@meta.data[,'orig.ident'])
op@meta.data[grep(pattern = "-1",rownames(op@meta.data)),'orig.ident'] <- 1  
op@meta.data[grep(pattern = "-2",rownames(op@meta.data)),'orig.ident'] <- 2
Idents(op) <- op@meta.data[,'orig.ident']

#######find cluster pca16 res0.7
op_cluster <- list()
for (i in c(8,11,16)){
	print(i)
	op_cluster[[i]] <- RunTSNE(object = op, features = VariableFeatures(object = op), dims = 1:i)
	op_cluster[[i]] <- FindNeighbors(object = op_cluster[[i]], features = VariableFeatures(object = op), reduction = 'pca', dims = 1:i)
	op_cluster[[i]] <- FindClusters(object = op_cluster[[i]], resolution = c(0.1,0.3,0.5,0.7,0.9,1.2))
	
	pdf(paste("tSNE_pca", i, ".pdf"))
	print(DimPlot(object = op_cluster[[i]], reduction = "tsne", group.by = 'orig.ident', pt.size = 1))	
	op_x <- cbind(op_cluster[[i]]@reductions$tsne@cell.embeddings, op_cluster[[i]]@meta.data[,5:10])
	for (j in 3:8) {
		print(j)
		p <-ggplot(op_x, aes(x = tSNE_1, y = tSNE_2, colour = op_x[,j])) + geom_point() +
			labs(title = colnames(op_x)[j]) +
			scale_color_discrete(name = colnames(op_x)[j]) 					
		print(p)
	}
	dev.off()	
}

#print(DimPlot(object = op_cluster[[8]], reduction = "tsne", group.by = 'orig.ident',pt.size=2)

######find markers pca16 res0.7
op <- RunTSNE(object = op, features = VariableFeatures(object = op), dims = 1:16)
op <- FindNeighbors(object = op, features = VariableFeatures(object = op), reduction = 'pca', dims = 1:16)
op <- FindClusters(object = op, resolution = 0.7)
	
op.markers <- FindAllMarkers(object = op, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
op.markers %>% group_by(cluster) %>% top_n(n = 2, wt = avg_logFC)

#write.csv(op.markers, file = paste("16","_PCA_marker.gene(endo).csv"))
write.table(op.markers, file = "marker.csv", sep=",")
write.csv(op.markers, file = "marker.csv", row.names = FALSE)



FeaturePlot(object = pbmc, features = c("MS4A1", "GNLY", "CD3E", "CD14", "FCER1A", "FCGR3A", "LYZ", 
    "PPBP", "CD8A"))


op <- RunTSNE(object = op, dims = 1:9)
# note that you can set `label = TRUE` or use the LabelClusters function to help label
# individual clusters
DimPlot(object = op, reduction = "tsne")

#############
op@assays$RNA@counts

top10 <- op.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC)
pdf('cluster_heatmap.pdf',width = 50,height = 25)

DoHeatmap(object = op, features = top10$gene, size = 10) + NoLegend() +
		theme(panel.grid.major = element_blank(),
        axis.text.x = element_text(color="black",face = "bold"),
        axis.text.y = element_text(color="black",face = "bold"))

dev.off() 

pdf('cluster_label.pdf', width = 10)
DimPlot(object = op, reduction = "tsne", label = TRUE)
dev.off() 


 
