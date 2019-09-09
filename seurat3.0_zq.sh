

library(dplyr)
library(Seurat)

counts <- Read10X("/home/zhangquan/project/singlecell/uterus/shhjjy/PB/PB_library/outs/filtered_feature_bc_matrix/")
cellinfo <- data.frame(barcode=colnames(counts),sample=c(rep("B1",8661),rep("B2",6569),rep("P1",7738),rep("P2",7978)))

# Initialize the Seurat object with the raw (non-normalized data)
uterus <- CreateSeuratObject(counts = counts, project = "uterus", min.cells = 3, min.features = 200)
uterus
#An object of class Seurat
#21263 features across 30883 samples within 1 assay
#Active assay: RNA (21263 features)

#The [[ operator can add columns to object metadata. This is a great place to stash QC stats
uterus[['percent.mt']] <- PercentageFeatureSet(uterus,pattern= '^MT-')

#Show QC metrics for the first 5 cells
head(uterus@meta.data,5)

#Visualize QC metrics as a violin plot
pdf('vlnplot_Feature_Count_mt1.pdf', width = 10, height = 10)
VlnPlot(uterus, features = c('nFeature_RNA', 'nCount_RNA', 'percent.mt'), ncol = 3)
dev.off()

#FeatureScatter is typically used to visualize feature-feature relationships, but can be used
#for anything caluculated by the object, i.e. columns in object metadata, PC scores etc.
plot1 <- FeatureScatter(object = uterus, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(object = uterus, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")

pdf('nCount_RNA.pdf', width = 10, height = 12)
CombinePlots(plots = list(plot1, plot2))
dev.off()

#add batch effect
uterus@meta.data$sample <- cellinfo$sample[match(rownames(uterus@meta.data),cellinfo$barcode)]
uterus@meta.data$batch <- unlist(sapply(as.character(uterus@meta.data$sample),function(x) unlist(strsplit(x,split="-"))[1] ))

#again show QC metrics for the first 5 cells
head(uterus@meta.data,5)

#filter cells
uterus <- subset(x = uterus, subset = nFeature_RNA > 200 & nFeature_RNA < 7500 & nCount_RNA > 0 & nCount_RNA < 75000 & percent.mt < 5)

#again Visualize QC metrics as a violin plot
pdf('vlnplot_nFeatures_nCount_mt2.pdf', width = 8, height = 8)
VlnPlot(uterus, features = c('nFeature_RNA', 'nCount_RNA', 'percent.mt'), ncol = 3)
dev.off()  

#Normalizing the data
uterus <- NormalizeData(object = uterus, normalization.method = "LogNormalize", scale.factor = 10000)
uterus <- NormalizeData(object = uterus)

#cellcycle genes(liuke)
cc.genes <- read.table("~/project/singlecell/uterus/regev_lab_cell_cycle_genes.txt",header=F,stringsAsFactors=F)
s.genes <- cc.genes[1:43,1]
g2m.genes <- cc.genes[44:97,1]
uterus <- CellCycleScoring(object = uterus, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)
head(uterus@meta.data)
#correction data
uterus <- ScaleData(object = uterus,vars.to.regress=c("batch","nCount_RNA","percent.mt","S.Score", "G2M.Score"), display.progress = FALSE)
#cellcycle genes(Seurat3.0)
exp.mat <- read.table("~/project/singlecell/uterus/nestorawa_forcellcycle_expressionMatrix.txt",header = TRUE,as.is = TRUE, row.names = 1)
s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes
uterus <- CellCycleScoring(uterus, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)
head(uterus[[]])
pdf('RidgePlotcellcycle.pdf',width = 10, height = 10)
RidgePlot(uterus,features = c('PCNA','TOP2A','MCM6','MKI67'),ncol = 2)
dev.off()
uterus <- ScaleData(uterus,vars.to.regress = c('S.Score','G2M.Score'), features = rownames(uterus))


#Identification of highly variable features (feature selection)
uterus <- FindVariableFeatures(object = uterus, selection.method = "vst", nfeatures = 7500)

# Identify the 10 most highly variable genes
top10 <- head(x = VariableFeatures(object = uterus), 10)
top428 <- head(x = VariableFeatures(object = uterus), 428)

# plot variable features with and without labels
plot1 <- VariableFeaturePlot(object = uterus)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE, xnudge = 0, ynudge= 0)

pdf("VariableFeatures.pdf", width = 15,height = 10)
CombinePlots(plots = list(plot1, plot2))
dev.off()

#scaling the data(细胞周期矫正后这步就不需要啦)
all.genes <- rownames(x = uterus)
uterus <- ScaleData(object = uterus, features = all.genes)
uterus <- ScaleData(object = uterus, vars.to.regress = "percent.mt")


###Perform linear dimensional reduction(use PCA)
uterus <- RunPCA(object = uterus, features = VariableFeatures(object = uterus))
print(x = uterus[["pca"]], dims = 1:5, nfeatures = 5)

pdf('vizdimpca.pdf')
VizDimLoadings(object = uterus, dims = 1:2, reduction = "pca")
dev.off()

pdf('pca.pdf')
DimPlot(object = uterus, reduction = "pca")
dev.off()


pdf('pcaheatmap.pdf', height = 20)
DimHeatmap(object = uterus, dims = 1:15, cells = 500, balanced = TRUE)
dev.off()

####uterus <- RunPCA(object = uterus,  pcs.compute = 30, genes.use = top1000,  do.print = TRUE, pcs.print = 1:5, genes.print = 5)

uterus <- JackStraw(object = uterus, num.replicate = 100)
uterus <- ScoreJackStraw(object = uterus, dims = 1:20)

pdf('JackStrawPlot.pdf')
JackStrawPlot(object = uterus, dims = 1:15)
ElbowPlot(object = uterus)
dev.off()


#####find cluster(3.0)#####
uterus <- FindNeighbors(object = uterus, features = VariableFeatures(object = uterus), reduction = 'pca', dims = 1:20)
uterus <- FindClusters(object = uterus, resolution = 0.5)
head(Idents(uterus),5)


######tNSE########
uterus <- RunTSNE(object = uterus, dims = 1:20, do.fast = TRUE, do.label = TRUE)

pdf("tSNE_20PCA(res.0.5).pdf", width=8)
DimPlot(object = uterus, reduction='tsne', pt.size = 1, label = TRUE)
dev.off()

tsne.res20 <- as.data.frame(cbind(uterus@reductions$tsne@cell.embeddings,uterus@meta.data))
colornames <- c('RNA_snn_res.0.5','sample')

p.pca20 <- sapply(colornames,function(x){
	p <- list(ggplot(tsne.res20,aes(x=tSNE_1,y=tSNE_2))+geom_point(aes_string(color = x),size = 1)+
			theme_bw()+labs(x="tSNE-1",y="tSNE-2",title=x,color=""))
	})
pdf("tSNE_pca_20_sample(res0.5).pdf",width=20,height=10)
plot_grid(plotlist=p.pca20,nrow=1)
dev.off()
uterus.markers <- FindAllMarkers(uterus, only.pos = T, min.pct = 0.25, logfc.threshold = 0.25)
write.csv(uterus.markers, file = "20_PCA_marker.gene(res0.5).csv")

new.cluster.ids <- c('Memory CD4+T','CD8+T','NK','B','Naive CD4+T','CD4+FOXP3+Treg','Th1','Macrophageocyte','DC','Platelet','CD14+Mono')
names(new.cluster.ids) <- levels(p.pca20)
p.pca20 <- RenameIdents(p.pca20,new.cluster.ids)
pdf('new.cluster.ids.pdf',width = 15, height = 15)
DimPlot(p.pca20, reduction = 'tsne', label = T, pt.size = 1)
dev.off()

save(uterus,'1000features.rData')



#####find cluster(3.0) use for circulation#####
# uterus_2 <- list()
# pc.nums <- c(7,10,15,20)
# sce_3 <- list()
# pc.nums <- c(5,8,15,20)
B12_4 <- list()
pc.nums <- c(5,10,15,20)
# uterus_5 <- list()
# pc.nums <- c(4,10,13,18)
for (i in pc.nums){
	print(i)
    B12_4[[i]] <- FindNeighbors(object = B12, features = VariableFeatures(object = B12), reduction = 'pca', dims = 1:i)
	B12_4[[i]] <- FindClusters(object = B12, resolution = c(0.4,0.5,0.6,0.8))
	B12_4[[i]] <- RunTSNE(object = B12_4[[i]], dims.use = 1:i, do.fast = TRUE)
}


# tsne.res <- list()
# tsne.res2 <- list()
tsne.res3 <- list()
# tsne.res4 <- list()
for(i in 1:length(pc.nums)){
	tsne.res3[[i]] <- as.data.frame(cbind(B12_4[[pc.nums[i]]]@reductions$tsne@cell.embeddings,B12_4[[pc.nums[i]]]@meta.data))
}
# colornames <- c("nFeature_RNA","nCount_RNA","res.0.4","res.0.5","res.0.6","res.0.8","sample")

colornames <- c('nFeature_RNA','nCount_RNA','RNA_snn_res.0.4','RNA_snn_res.0.5','RNA_snn_res.0.6','RNA_snn_res.0.8','Sample')
p.cluster1 <- sapply(colornames,function(x){
	p <- list(ggplot(tsne.res3[[1]],aes(x=tSNE_1,y=tSNE_2))+geom_point(aes_string(color = x),size=0.1)+
			theme_bw()+labs(x="tSNE-1",y="tSNE-2",title=x,color=""))
	})
#####liuke
p.cluster1 <- sapply(colornames,function(x){
	p <- list(ggplot(tsne.res3[[1]],aes(x=tSNE_1,y=tSNE_2))+geom_point(aes(color=tsne.res3[[1]][,x]),size=0.1)+
			theme_bw()+labs(x="tSNE-1",y="tSNE-2",title=x,color=""))
	})

pdf("tSNE_pca_5.pdf",width=20,height=10)
plot_grid(plotlist=p.cluster1,nrow=2)
dev.off()



p.cluster2 <- sapply(colornames,function(x){
	p <- list(ggplot(tsne.res3[[2]],aes(x=tSNE_1,y=tSNE_2))+geom_point(aes_string(color = x),size=0.1)+
			theme_bw()+labs(x="tSNE-1",y="tSNE-2",title=x,color=""))
	})
pdf("tSNE_pca_10.pdf",width=20,height=8)
plot_grid(plotlist=p.cluster2,nrow=2)
dev.off()




p.cluster3 <- sapply(colornames,function(x){
	p <- list(ggplot(tsne.res3[[3]],aes(x=tSNE_1,y=tSNE_2))+geom_point(aes_string(color = x),size=0.1)+
			theme_bw()+labs(x="tSNE-1",y="tSNE-2",title=x,color=""))
	})
pdf("tSNE_pca_15.pdf",width=20,height=8)
plot_grid(plotlist=p.cluster3,nrow=2)
dev.off()


p.cluster4 <- sapply(colornames,function(x){
	p <- list(ggplot(tsne.res3[[4]],aes(x=tSNE_1,y=tSNE_2))+geom_point(aes_string(color = x),size=0.1)+
			theme_bw()+labs(x="tSNE-1",y="tSNE-2",title=x,color=""))
	})
pdf("tSNE_pca_20.pdf",width=20,height=8)
plot_grid(plotlist=p.cluster4,nrow=2)
dev.off()


#=====================================cellsubset=====================================#(action content)
uterus <- uterus_4[[20]]
uterus@active.ident <- factor(uterus@meta.data$RNA_snn_res.0.5,levels=c(0:10))
names(uterus@active.ident) <- rownames(uterus@meta.data)
pdf("tSNE_pca_20(res0.5).pdf")
TSNEPlot(object = uterus, no.legend=F, do.label=T, pt.size=0.5)
dev.off()
uterus.markers0 <- FindAllMarkers(object = uterus, only.pos = TRUE, min.pct = 0.25,
      logfc.thresh = 0.25)
write.csv(uterus.markers0, file = "20_PCA_marker.gene(res0.5).csv")
markgenes <- c("CD4","IL7R","CD8A","CD8B","NKG7","S100A4","GZMK","GZMB","GZMH","GZMA","LGALS1","CD69","CST7",
	"CCR7","IL2RA","FOXP3")
p.markers <- sapply(markgenes,function(x){
	print(x)
	p <- list(ggplot(tsne.res3[[4]],aes(x=tSNE_1,y=tSNE_2)) + 
		      geom_point(aes_string(color = t(uterus[['RNA']]@scale.data)[,x],size=0.1)) +
			  scale_color_gradient2(low = "grey",high = "darkred",mid="grey")+
			  theme_bw()+labs(x="tSNE-1",y="tSNE-2",title=x,color=""))
	})
pdf("tSNE_pca_20(res0.5_markgenes).pdf",width=30,height=12)
plot_grid(plotlist=p.markers,ncol=6)
dev.off()
######uterus[['RNA']]@scale.data[,x]
######t(sce@scale.data))[,x]


#####Finding differentially expressed features(cluster biomarkers)#####
cluster1.markers <- FindMarkers(uterus, ident.1 = 5, ident.2 = c(0,3), min.pct = 0.25)
head(cluster5.markers, n = 5)

uterus.markers <- FindAllMarkers(uterus, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
uterus.markers %>% group_by(cluster) %>% top_n(n = 2, wt = avg_logFC)
write.csv(uterus.markers,file = 'marker cluster.csv',row.names = FALSE)

new.cluster.ids <- c('Memory CD4+T','CD8+T','NK','B','Naive CD4+T','CD4+FOXP3+Treg','Th1','Macrophageocyte','DC','Platelet','CD14+Mono')
names(new.cluster.ids) <- levels(uterus)
uterus <- RenameIdents(uterus,new.cluster.ids)
pdf('new.cluster.ids.pdf',width = 15, height = 15)
DimPlot(uterus, reduction = 'tsne', label = T, pt.size = 1)
dev.off()




























