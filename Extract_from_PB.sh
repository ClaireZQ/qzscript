library(Seurat)
library(cowplot)

ctrl.data <- read.table(file = "../data/immune_control_expression_matrix.txt.gz", sep = "\t")
stim.data <- read.table(file = "../data/immune_stimulated_expression_matrix.txt.gz", sep = "\t")

# Set up control object
ctrl <- CreateSeuratObject(counts = ctrl.data, project = "IMMUNE_CTRL", min.cells = 5)



library(Seurat)
library(cowplot)

index <- rownames(uterus@meta.data[which(uterus@meta.data$sample =='BA' | uterus@meta.data$sample == 'B'),])
bab_counts <- uterus@assays$RNA@counts[, index]


# Set up control object
BAB <- CreateSeuratObject(counts = bab_counts, project = "BAB", min.cells = 3, min.features=200)


PA <- rownames(uterus@meta.data[which(uterus@meta.data$sample =='PA'),])
P <- rownames(uterus@meta.data[which(uterus@meta.data$sample == 'P'),])
#add sample
BAB@meta.data$sample <- c(rep('B',4517),rep('BA',11274))

#add Batch
BAB[['batch']] <- BAB[['sample']]


#####two matrices comparison
ctrl <- uterus
###liuke data
ctrl$stim <- "CTRL"
ctrl <- subset(ctrl, subset = nFeature_RNA > 500)
ctrl <- NormalizeData(ctrl, verbose = FALSE)
ctrl <- FindVariableFeatures(ctrl, selection.method = "vst", nfeatures = 428)

###zq data
uterus$stim <- "STIM"
uterus <- subset(uterus, subset = nFeature_RNA > 500)
uterus <- NormalizeData(uterus, verbose = FALSE)
uterus <- FindVariableFeatures(uterus, selection.method = "vst", nfeatures = 428)
###perform integration
immune.anchors <- FindIntegrationAnchors(object.list = list(ctrl, uterus), dims = 1:20)  (slowly)
immune.combined <- IntegrateData(anchorset = immune.anchors, dims = 1:20)


###perform an integrated analysis
DefaultAssay(immune.combined) <- 'integrated'
#Run the standard workflow for visualization and clustering
immune.combined <- ScaleData(immune.combined, verbose = FALSE)
immune.combined <- RunPCA(immune.combined, npcs = 30, verbose = FALSE)
#t-SNE and Clustering
immune.combined <- RunTSNE(immune.combined, reduction = 'pca', dims = 1:20)
immune.combined <- FindNeighbors(immune.combined, reduction = 'pca', dims = 1:20)
immune.combined <- FindClusters(immune.combined, resolution = 0.5)

#Visualization
pdf('combined1.pdf', width = 15, height = 8)
p1 <- DimPlot(immune.combined, reduction = 'tsne', group.by = 'stim')
p2 <- DimPlot(immune.combined, reduction = 'tsne', label = TRUE)
plot_grid(p1,p2)
dev.off()


#To visulize the two conditions side-by-side, we can use the split.by argument to show each condition colored by cluster.
pdf('split_liu_zh1.pdf', width = 15, height = 8)
DimPlot(immune.combined, reduction = 'tsne', label = TRUE, split.by = 'stim')
dev.off()

#Identify conserved cell type markers
DefaultAssay(immune.combined) <- 'RNA'
cl10.markers <- FindConservedMarkers(immune.combined, ident.1 = 10, grouping.var = 'stim', verbose = FALSE)
head(cl10.markers)
write.csv(cl10.markers, file = "10ConservedMarkers.csv")
cl00.markers <- FindConservedMarkers(immune.combined, ident.1 = 0, grouping.var = 'stim', verbose = FALSE)
head(cl00.markers)
write.csv(cl00.markers, file = "00ConservedMarkers.csv")
cl01.markers <- FindConservedMarkers(immune.combined, ident.1 = 1, grouping.var = 'stim', verbose = FALSE)
head(cl01.markers)
write.csv(cl01.markers, file = "01ConservedMarkers.csv")
cl02.markers <- FindConservedMarkers(immune.combined, ident.1 = 2, grouping.var = 'stim', verbose = FALSE)
head(cl02.markers)
write.csv(cl02.markers, file = "02ConservedMarkers.csv")
cl03.markers <- FindConservedMarkers(immune.combined, ident.1 = 3, grouping.var = 'stim', verbose = FALSE)
head(cl03.markers)
write.csv(cl03.markers, file = "03ConservedMarkers.csv")
cl04.markers <- FindConservedMarkers(immune.combined, ident.1 = 4, grouping.var = 'stim', verbose = FALSE)
head(cl04.markers)
write.csv(cl04.markers, file = "04ConservedMarkers.csv")
cl05.markers <- FindConservedMarkers(immune.combined, ident.1 = 5, grouping.var = 'stim', verbose = FALSE)
head(cl05.markers)
write.csv(cl05.markers, file = "05ConservedMarkers.csv")
cl06.markers <- FindConservedMarkers(immune.combined, ident.1 = 6, grouping.var = 'stim', verbose = FALSE)
head(cl06.markers)
write.csv(cl06.markers, file = "06ConservedMarkers.csv")
cl13.markers <- FindConservedMarkers(immune.combined, ident.1 = 13, grouping.var = 'stim', verbose = FALSE)
head(cl13.markers)
write.csv(cl13.markers, file = "13ConservedMarkers.csv")

pdf('immune.combined.pdf', width = 15, height = 10)
FeaturePlot(immune.combined, features = c('CD8B', 'CD8A', 'CD4', 'NKG7', 'SELL', 'ADAM19', 'GNLY', 'CD69', 'IL7R'),
	        min.cutoff = 'q9')

immune.combined <- RenameIdents(immune.combined, '0' = 'cluster0', '1' = 'cluster1', '2' = 'cluster2', '3' = 'cluster3',
	                            '4' = 'cluster4', '5' = 'cluster5', '6' = 'cluster6', '7' = 'cluster7', '8' = 'cluster8',
	                            '9' = 'cluster9', '10' = 'cluster10', '11' = 'cluster11', '12' = 'cluster12', '13' = 'cluster13')


pdf('immune.combined1.pdf', width = 8, height = 8)
DimPlot(immune.combined, label = TRUE)
dev.off()

Idents(immune.combined) <- factor(Idents(immune.combined),levels = c('cluster13','cluster12','cluster11','cluster10',
	                        'cluster9','cluster8','cluster7','cluster6','cluster5','cluster4','cluster3','cluster2',
	                        'cluster1','cluster0')) 

markers.to.plot <- c('CD40LG','IL7R','NKG7','CD8A','CD8B','SOX4','SELL','FOXP3',
	                 'CTLA4','HSPA1A','HSPA1B','LEF1','CCR7','CCL5',
	                 'GZMB','FGFBP2','FCGR3A','C1QA','C1QB','RPS12','RPL32',
	                 'STMN1','NCAM1','KRT81','PSG3','S100A16','IDO1','CLEC9A')

pdf('immune_Dotplot1.pdf', width = 15, height = 12)
DimPlot(immune.combined, features = rev(markers.to.plot), cols = BlueAndRed(14), dot.scale = 8, split.by = 'stim')+ RotatedAxis()  ####(错误的语法,但是解决问题是BlueAndRed(14))
DotPlot(immune.combined, features = rev(markers.to.plot), cols = c('blue', 'red'), dot.scale = 8, split.by = 'stim')+ RotatedAxis()
