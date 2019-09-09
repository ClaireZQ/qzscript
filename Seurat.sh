library(Seurat)

#read counts from cellranger results
counts <- Read10X("/home/liuke/project/single.cell/uterus/aggr/outs/filtered_feature_bc_matrix/hg19/")
cellinfo <- data.frame(barcode=colnames(counts),sample=c(rep("PIL-B",4603),rep("PIL-P",4072),rep("TIL-B",11838),rep("TIL-P",12790)))

# counts_p <- counts[,which(cellinfo$sample %in% c("TIL-P","PIL-P"))]
# counts_b <- counts[,which(cellinfo$sample %in% c("TIL-B","PIL-B"))]


sce <- CreateSeuratObject(raw.data = counts,min.genes=100,min.cells=3,project = "uterus")

#calculate mito genes
mito.genes <- grep(pattern = "^MT-", x = rownames(x = sce@data), value = TRUE)
percent.mito <- Matrix::colSums(sce@raw.data[mito.genes, ])/Matrix::colSums(sce@raw.data)
sce <- AddMetaData(object = sce, metadata = percent.mito, col.name = "percent.mito")


pdf("VlnPlot_nGene.nUMI.percent.mito2.pdf",width=13)
VlnPlot(object = sce, features.plot = c("nGene", "nUMI","percent.mito"), nCol = 3)
dev.off()

#add batch effects
sce@meta.data$sample <- cellinfo$sample[match(rownames(sce@meta.data),cellinfo$barcode)]
sce@meta.data$batch <- unlist(sapply(as.character(sce@meta.data$sample),function(x) unlist(strsplit(x,split="-"))[1] ))

# filter cells
sce <- FilterCells(object = sce, subset.names = c("nUMI", "percent.mito"),
    low.thresholds = c(200, -Inf), high.thresholds = c(8000, 0.05))
# normalized
sce <- NormalizeData(object = sce, normalization.method = "LogNormalize",
    scale.factor = 10000)
# cellcycle genes
cc.genes <- read.table("~/project/single.cell/regev_lab_cell_cycle_genes.txt",header=F,stringsAsFactors=F)
s.genes <- cc.genes[1:43,1]
g2m.genes <- cc.genes[44:97,1]
sce <- CellCycleScoring(object = sce, s.genes = s.genes, g2m.genes = g2m.genes,
    set.ident = TRUE)
head(sce@meta.data)


#correction data
sce <- ScaleData(object = sce,vars.to.regress=c("batch","nUMI","percent.mito","S.Score", "G2M.Score"), display.progress = FALSE)


#==================================== hv genes ==============================
# select genes
sce <- FindVariableGenes(object = sce, mean.function = ExpMean,dispersion.function = LogVMR,
	x.low.cutoff = 0.0125,x.high.cutoff = 3,y.cutoff = 0.5)

# hv.genes <- sce@var.genes
hv.genes <- head(rownames(sce@hvg.info),2000)
# hv.genes <- head(rownames(sce@hvg.info),3000)

###############runPCA#########################

sce <- RunPCA(object = sce,  pcs.compute = 30, genes.use = hv.genes,  do.print = TRUE, pcs.print = 1:5, genes.print = 5)
#############select pca#######################
sce <- JackStraw(object = sce, num.replicate = 100, num.pc=30)

pdf("JackStrawPlot.pdf")
JackStrawPlot(object = sce, PCs = 1:30)
dev.off()
pdf("PCElbowPlot.pdf")
PCElbowPlot(object = sce)
dev.off()

######find cluster################
# sce_2 <- list()
# pc.nums <- c(4,10,15,20)
# sce_3 <- list()
# pc.nums <- c(5,8,15,20)
sce_4 <- list()
pc.nums <- c(5,8,15,20)
# sce_5 <- list()
# pc.nums <- c(4,10,13,18)
for (i in pc.nums){
	print(i)
	sce_4[[i]] <- FindClusters(object = sce, reduction.type = "pca", dims.use = 1:i,
    	resolution = c(0.4,0.5,0.6,0.8), print.output = 0, save.SNN = F)
	sce_4[[i]] <- RunTSNE(object = sce_4[[i]], dims.use = 1:i, do.fast = TRUE)
}


# markgene.data <- as.data.frame(t(sce.lung_0@scale.data[markgene,]),stringsAsFactors=F)

# tsne.res <- list()
# tsne.res2 <- list()
tsne.res3 <- list()
# tsne.res4 <- list()
for(i in 1:length(pc.nums)){
	tsne.res3[[i]] <- as.data.frame(cbind(sce_4[[pc.nums[i]]]@dr$tsne@cell.embeddings,sce_4[[pc.nums[i]]]@meta.data))
}
colornames <- c("nGene","nUMI","res.0.4","res.0.5","res.0.6","res.0.8","sample")
p.cluster1 <- sapply(colornames,function(x){
	p <- list(ggplot(tsne.res3[[1]],aes(x=tSNE_1,y=tSNE_2))+geom_point(aes(color=tsne.res3[[1]][,x]),size=0.1)+
			theme_bw()+labs(x="tSNE-1",y="tSNE-2",title=x,color=""))
	})
pdf("tSNE_pca_5.pdf",width=20,height=8)
plot_grid(plotlist=p.cluster1,nrow=2)
dev.off()
p.cluster2 <- sapply(colornames,function(x){
	p <- list(ggplot(tsne.res3[[2]],aes(x=tSNE_1,y=tSNE_2))+geom_point(aes(color=tsne.res3[[2]][,x]),size=0.1)+
			theme_bw()+labs(x="tSNE-1",y="tSNE-2",title=x,color=""))
	})
pdf("tSNE_pca_8.pdf",width=20,height=8)
plot_grid(plotlist=p.cluster2,nrow=2)
dev.off()
p.cluster3 <- sapply(colornames,function(x){
	p <- list(ggplot(tsne.res3[[3]],aes(x=tSNE_1,y=tSNE_2))+geom_point(aes(color=tsne.res3[[3]][,x]),size=0.1)+
			theme_bw()+labs(x="tSNE-1",y="tSNE-2",title=x,color=""))
	})
pdf("tSNE_pca_15.pdf",width=20,height=8)
plot_grid(plotlist=p.cluster3,nrow=2)
dev.off()
p.cluster4 <- sapply(colornames,function(x){
	p <- list(ggplot(tsne.res3[[4]],aes(x=tSNE_1,y=tSNE_2))+geom_point(aes(color=tsne.res3[[4]][,x]),size=0.1)+
			theme_bw()+labs(x="tSNE-1",y="tSNE-2",title=x,color=""))
	})
pdf("tSNE_pca_20.pdf",width=20,height=8)
plot_grid(plotlist=p.cluster4,nrow=2)
dev.off()

#=========================================== cellsubset ===============================
sce <- sce_2[[20]]
sce@ident <- factor(sce@meta.data$res.0.5,levels=c(0:14))
names(sce@ident) <- rownames(sce@meta.data)
pdf("tSNE_pca_20(res0.5).pdf")
TSNEPlot(object = sce, no.legend=F, do.label=T, pt.size=0.5)
dev.off()
sce.markers0 <- FindAllMarkers(object = sce, only.pos = TRUE, min.pct = 0.25,
      thresh.use = 0.25)
write.csv(sce.markers0, file = "20_PCA_marker.gene(res0.5).csv")
markgenes <- c("CD3D","CD4","CD8A","CD8B","FOS","IL7R","GZMK","CCR7","IL2RA","HSPA6","FOXP3","CTLA4","IFI44L",
	"FGFBP2","NKG7","SPP1","C1QA","GNLY","STMN1","MKI67","FABP4","S100A8","IL8")
p.markers <- sapply(markgenes,function(x){
	print(x)
	p <- list(ggplot(tsne.res2[[4]],aes(x=tSNE_1,y=tSNE_2))+geom_point(aes(color=t(sce@scale.data)[,x]),size=0.1)+
				scale_color_gradient2(low = "grey",high = "darkred",mid="grey")+
				theme_bw()+labs(x="tSNE-1",y="tSNE-2",title=x,color=""))
	})
pdf("tSNE_pca_20(res0.5_markgenes).pdf",width=30,height=12)
plot_grid(plotlist=p.markers,ncol=6)
dev.off()


sce@ident <- factor(sce@meta.data$res.0.8,levels=c(0:16))
names(sce@ident) <- rownames(sce@meta.data)
pdf("tSNE_pca_20(res0.8).pdf")
TSNEPlot(object = sce, no.legend=F, do.label=T, pt.size=0.5)
dev.off()
sce.markers1 <- FindAllMarkers(object = sce, only.pos = TRUE, min.pct = 0.25,
      thresh.use = 0.25)
write.csv(sce.markers1, file = "20_PCA_marker.gene(res0.8).csv")
markgenes <- c("CD3D","CD4","CD8A","CD8B","GZMK","IL7R","CCR7","IL2RA","KLF2","FOS","IFI44L","HSPA6","BATF","ZNF683",
	"FGFBP2","NKG7","SPP1","C1QA","C1QB","GNLY","STMN1","KLRB1","CXCL13","S100A8","S100A9","FABP4")
p.markers <- sapply(markgenes,function(x){
	print(x)
	p <- list(ggplot(tsne.res2[[4]],aes(x=tSNE_1,y=tSNE_2))+geom_point(aes(color=t(sce@scale.data)[,x]),size=0.1)+
				scale_color_gradient2(low = "grey",high = "darkred",mid="grey")+
				theme_bw()+labs(x="tSNE-1",y="tSNE-2",title=x,color=""))
	})
pdf("tSNE_pca_20(res0.8_markgenes).pdf",width=30,height=16)
plot_grid(plotlist=p.markers,ncol=6)
dev.off()


# tnf expression
tsne <- tsne.res[[4]]
tsne$TNF <- sce@scale.data["TNF",]
p.tnf <- ggplot(tsne,aes(x=tSNE_1,y=tSNE_2))+geom_point(aes(color=TNF),size=0.1)+
	theme_bw()+scale_color_gradient2(low = "grey",high = "darkred",mid="grey")+labs(color="TNF")

tsne.t <- tsne.res[[4]][tsne.res[[4]]$batch=="TIL",]
tsne.t$TNF <- sce@scale.data["TNF",match(rownames(tsne.t),colnames(sce@scale.data))]
p.tnf.t <- ggplot(tsne.t,aes(x=tSNE_1,y=tSNE_2))+geom_point(aes(color=TNF),size=0.1)+
	theme_bw()+scale_color_gradient2(low = "grey",high = "darkred",mid="grey")+labs(color="TIL-TNF")

tsne.p <- tsne.res[[4]][tsne.res[[4]]$batch=="PIL",]
tsne.p$TNF <- sce@scale.data["TNF",match(rownames(tsne.p),colnames(sce@scale.data))]
p.tnf.p <- ggplot(tsne.p,aes(x=tSNE_1,y=tSNE_2))+geom_point(aes(color=TNF),size=0.1)+
	theme_bw()+scale_color_gradient2(low = "grey",high = "darkred",mid="grey")+labs(color="PIL-TNF")
pdf("TNF.expression.pdf",width=10,height=8)
plot_grid(p.tnf,p.tnf.t,p.tnf.p,ncol=2)
dev.off()



tsne$TNF_status <- unlist(sapply(tsne$TNF,function(x){
	if(x>0){
		return("Expr")
	} else return("Non-expr")
	}))
tnf.per <- as.data.frame(table(tsne$sample,tsne$TNF_status)/rowSums(table(tsne$sample,tsne$TNF_status)),stringsAsFactors=F)
p.tnf.per <- ggplot(tnf.per[1:4,],aes(x=Var1,y=Freq))+geom_bar(aes(fill=Var1),color="grey30",size=1,stat="identity",width=0.8)+
	 scale_fill_manual(values=c('#d35e14', '#9067a7', '#ab6857', '#7793cb'))+
	 theme_bw()+labs(x="",y="Percent of TNF-expressed cells",fill="")

tsne.sel <- tsne[tsne$TNF>0,]
my_comparisons <- list(c("PIL-P","PIL-B"),c("PIL-P","TIL-P"),c("PIL-P","TIL-B"))
p.tnf.exp <- ggplot(tsne.sel,aes(x=sample,y=TNF))+geom_boxplot(aes(fill=sample),color="grey30",size=1,width=0.8)+
	scale_fill_manual(values=c('#d35e14', '#9067a7', '#ab6857', '#7793cb'))+
	stat_compare_means(comparisons=my_comparisons)+
	theme_bw()+labs(x="",y="TNF expression",fill="")
pdf("TNF.expression2.pdf",width=8.5,height=3.5)
plot_grid(p.tnf.per,p.tnf.exp,ncol=2)
dev.off()



tsne$FOXP1 <- sce@scale.data["FOXP1",]
tsne$FOXP1_status <- unlist(sapply(tsne$FOXP1,function(x){
	if(x>0){
		return("Expr")
	} else return("Non-expr")
	}))
FOXP1.per <- as.data.frame(table(tsne$sample,tsne$FOXP1_status)/rowSums(table(tsne$sample,tsne$FOXP1_status)),stringsAsFactors=F)
p.FOXP1.per <- ggplot(FOXP1.per[1:4,],aes(x=Var1,y=Freq))+geom_bar(aes(fill=Var1),color="grey30",size=1,stat="identity",width=0.8)+
	 scale_fill_manual(values=c('#d35e14', '#9067a7', '#ab6857', '#7793cb'))+
	 theme_bw()+labs(x="",y="Percent of FOXP1-expressed cells",fill="")

tsne.sel <- tsne[tsne$FOXP1>0,]
p.FOXP1.exp <- ggplot(tsne.sel,aes(x=sample,y=FOXP1))+geom_boxplot(aes(fill=sample),color="grey30",size=1,width=0.8)+
	scale_fill_manual(values=c('#d35e14', '#9067a7', '#ab6857', '#7793cb'))+
	stat_compare_means(comparisons=my_comparisons)+
	theme_bw()+labs(x="",y="FOXP1 expression",fill="")
pdf("FOXP1.expression2.pdf",width=8.5,height=3.5)
plot_grid(p.FOXP1.per,p.FOXP1.exp,ncol=2)
dev.off()


# CD4 and CD8
tsne$CD <- rep("others",dim(tsne)[1])
tsne$CD[which(tsne$res.0.5 %in% c(0,2,3,4,5,9,10))] <- "CD4"
tsne$CD[which(tsne$res.0.5 %in% c(1,6,7,8,11))] <- "CD8"

# CD4
tsne.cd4 <- tsne[tsne$CD=="CD4",]

#tnf
tnf.per.cd4 <- as.data.frame(table(tsne.cd4$sample,tsne.cd4$TNF_status)/rowSums(table(tsne.cd4$sample,tsne.cd4$TNF_status)),stringsAsFactors=F)
p.tnf.per.cd4 <- ggplot(tnf.per.cd4[1:4,],aes(x=Var1,y=Freq))+geom_bar(aes(fill=Var1),color="grey30",size=1,stat="identity",width=0.8)+
	 scale_fill_manual(values=c('#d35e14', '#9067a7', '#ab6857', '#7793cb'))+
	 theme_bw()+labs(x="",y="Percent of TNF-expressed cells",fill="")

tsne.sel.cd4 <- tsne.cd4[tsne.cd4$TNF>0,]
my_comparisons <- list(c("PIL-P","PIL-B"),c("PIL-P","TIL-P"),c("PIL-P","TIL-B"),c("TIL-B","TIL-P"),c("TIL-B","PIL-B"))
p.tnf.exp.cd4 <- ggplot(tsne.sel.cd4,aes(x=sample,y=TNF))+geom_boxplot(aes(fill=sample),color="grey30",size=1,width=0.8)+
	scale_fill_manual(values=c('#d35e14', '#9067a7', '#ab6857', '#7793cb'))+
	stat_compare_means(comparisons=my_comparisons)+
	theme_bw()+labs(x="",y="TNF expression",fill="")
pdf("TNF.expression2(cd4).pdf",width=8.5,height=3.5)
plot_grid(p.tnf.per.cd4,p.tnf.exp.cd4,ncol=2)
dev.off()

#foxp1
FOXP1.per.cd4 <- as.data.frame(table(tsne.cd4$sample,tsne.cd4$FOXP1_status)/rowSums(table(tsne.cd4$sample,tsne.cd4$FOXP1_status)),stringsAsFactors=F)
p.FOXP1.per.cd4 <- ggplot(FOXP1.per.cd4[1:4,],aes(x=Var1,y=Freq))+geom_bar(aes(fill=Var1),color="grey30",size=1,stat="identity",width=0.8)+
	 scale_fill_manual(values=c('#d35e14', '#9067a7', '#ab6857', '#7793cb'))+
	 theme_bw()+labs(x="",y="Percent of FOXP1-expressed cells",fill="")

tsne.sel.cd4 <- tsne.cd4[tsne.cd4$FOXP1>0,]
my_comparisons <- list(c("PIL-P","PIL-B"),c("PIL-P","TIL-P"),c("PIL-P","TIL-B"),c("TIL-B","TIL-P"),c("TIL-B","PIL-B"))
p.FOXP1.exp.cd4 <- ggplot(tsne.sel.cd4,aes(x=sample,y=FOXP1))+geom_boxplot(aes(fill=sample),color="grey30",size=1,width=0.8)+
	scale_fill_manual(values=c('#d35e14', '#9067a7', '#ab6857', '#7793cb'))+
	stat_compare_means(comparisons=my_comparisons)+
	theme_bw()+labs(x="",y="FOXP1 expression",fill="")
pdf("FOXP1.expression2(cd4).pdf",width=8.5,height=3.5)
plot_grid(p.FOXP1.per.cd4,p.FOXP1.exp.cd4,ncol=2)
dev.off()



# CD8
tsne.cd8 <- tsne[tsne$CD=="CD8",]

#tnf
tnf.per.cd8 <- as.data.frame(table(tsne.cd8$sample,tsne.cd8$TNF_status)/rowSums(table(tsne.cd8$sample,tsne.cd8$TNF_status)),stringsAsFactors=F)
p.tnf.per.cd8 <- ggplot(tnf.per.cd8[1:4,],aes(x=Var1,y=Freq))+geom_bar(aes(fill=Var1),color="grey30",size=1,stat="identity",width=0.8)+
	 scale_fill_manual(values=c('#d35e14', '#9067a7', '#ab6857', '#7793cb'))+
	 theme_bw()+labs(x="",y="Percent of TNF-expressed cells",fill="")

tsne.sel.cd8 <- tsne.cd8[tsne.cd8$TNF>0,]
my_comparisons <- list(c("PIL-P","PIL-B"),c("PIL-P","TIL-P"),c("PIL-P","TIL-B"),c("TIL-B","TIL-P"),c("TIL-B","PIL-B"))
p.tnf.exp.cd8 <- ggplot(tsne.sel.cd8,aes(x=sample,y=TNF))+geom_boxplot(aes(fill=sample),color="grey30",size=1,width=0.8)+
	scale_fill_manual(values=c('#d35e14', '#9067a7', '#ab6857', '#7793cb'))+
	stat_compare_means(comparisons=my_comparisons)+
	theme_bw()+labs(x="",y="TNF expression",fill="")
pdf("TNF.expression2(cd8).pdf",width=8.5,height=3.5)
plot_grid(p.tnf.per.cd8,p.tnf.exp.cd8,ncol=2)
dev.off()

#foxp1
FOXP1.per.cd8 <- as.data.frame(table(tsne.cd8$sample,tsne.cd8$FOXP1_status)/rowSums(table(tsne.cd8$sample,tsne.cd8$FOXP1_status)),stringsAsFactors=F)
p.FOXP1.per.cd8 <- ggplot(FOXP1.per.cd8[1:4,],aes(x=Var1,y=Freq))+geom_bar(aes(fill=Var1),color="grey30",size=1,stat="identity",width=0.8)+
	 scale_fill_manual(values=c('#d35e14', '#9067a7', '#ab6857', '#7793cb'))+
	 theme_bw()+labs(x="",y="Percent of FOXP1-expressed cells",fill="")

tsne.sel.cd8 <- tsne.cd8[tsne.cd8$FOXP1>0,]
my_comparisons <- list(c("PIL-P","PIL-B"),c("PIL-P","TIL-P"),c("PIL-P","TIL-B"),c("TIL-B","TIL-P"),c("TIL-B","PIL-B"))
p.FOXP1.exp.cd8 <- ggplot(tsne.sel.cd8,aes(x=sample,y=FOXP1))+geom_boxplot(aes(fill=sample),color="grey30",size=1,width=0.8)+
	scale_fill_manual(values=c('#d35e14', '#9067a7', '#ab6857', '#7793cb'))+
	stat_compare_means(comparisons=my_comparisons)+
	theme_bw()+labs(x="",y="FOXP1 expression",fill="")
pdf("FOXP1.expression2(cd8).pdf",width=8.5,height=3.5)
plot_grid(p.FOXP1.per.cd8,p.FOXP1.exp.cd8,ncol=2)
dev.off()


pdf("TNF.and.FOXP1.expression.pdf",width=20,height=12)
plot_grid(p.tnf.per,p.tnf.exp,p.FOXP1.per,p.FOXP1.exp,
	p.tnf.per.cd4,p.tnf.exp.cd4,p.FOXP1.per.cd4,p.FOXP1.exp.cd4,
	p.tnf.per.cd8,p.tnf.exp.cd8,p.FOXP1.per.cd8,p.FOXP1.exp.cd8,ncol=4)
dev.off()



#=================== comparison of foxp1 + and - ========================




tsne$FOXP1_sample <- paste(tsne$FOXP1_status,tsne$sample,sep="_")
comparisons <- list(c("Expr_TIL-P","Non-expr_TIL-P"),
					c("Expr_PIL-P","Non-expr_PIL-P"),
					c("Expr_PIL-P","Expr_TIL-P"),
					c("Non-expr_PIL-P","Non-expr_TIL-P"))

seurat_object <- sce
seurat_object@ident <- factor(tsne$FOXP1_sample)
names(seurat_object@ident) <- rownames(seurat_object)
foxp1.compare <- lapply(comparisons,function(x){
	FindMarkers(object = seurat_object, ident.1 = x[1], ident.2 = x[2])
	})

for(i in 1:4){
	write.xlsx(foxp1.compare[[i]],file="foxp1_compare.xlsx",sheetName=names(foxp1.compare)[i],col.names=TRUE, row.names=TRUE, append=T, showNA=TRUE, password=NULL)
}

# #=========================================== cellsubset ===============================
# sce <- sce_5[[18]]
# sce@ident <- factor(sce@meta.data$res.0.5,levels=c(0:12))
# names(sce@ident) <- rownames(sce@meta.data)
# pdf("tSNE_pca_18(res0.5).pdf")
# TSNEPlot(object = sce, no.legend=F, do.label=T, pt.size=0.5)
# dev.off()
# sce.markers <- FindAllMarkers(object = sce, only.pos = TRUE, min.pct = 0.25,
#       thresh.use = 0.25)
# write.csv(sce.markers, file = "18_PCA_marker.gene(res0.5).csv")
# markgenes <- c("CD3D","CD4","CD8A","CD8B","FOS","IL7R","GZMK","CCR7","IL2RA","HSPA6","FOXP3","CTLA4","IFI44L",
# 	"FGFBP2","NKG7","SPP1","C1QA","GNLY","STMN1","MKI67","FABP4","S100A8","IL8")
# p.markers <- sapply(markgenes,function(x){
# 	print(x)
# 	p <- list(ggplot(tsne.res4[[4]],aes(x=tSNE_1,y=tSNE_2))+geom_point(aes(color=t(sce@scale.data)[,x]),size=0.1)+
# 				scale_color_gradient2(low = "grey",high = "darkred",mid="grey")+
# 				theme_bw()+labs(x="tSNE-1",y="tSNE-2",title=x,color=""))
# 	})
# pdf("tSNE_pca_18(res0.5_markgenes).pdf",width=30,height=12)
# plot_grid(plotlist=p.markers,ncol=6)
# dev.off()


# #===================================== monocle ================================
# HSMM_p <- newCellDataSet(as.matrix(counts_p), lowerDetectionLimit = 0.5, expressionFamily=negbinomial.size())
# HSMM <- estimateSizeFactors(HSMM)

# HSMM <- detectGenes(HSMM, min_expr = 0.1)
# expressed_genes <- row.names(subset(fData(HSMM),
#     num_cells_expressed >= 50))