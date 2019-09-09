pb@meta.data$nSample_RNA <- c(rep('B1',8661),rep('B2',6569),rep('P1',7738),rep('P2',7978))

#如果给细胞分群后重新画图,用以下的语句可以去除分群信息,群信息主要传给了Idents
> DimPlot(sce,reduction='pca',label = F,group.by=NULL)
> unique(Idents(sce))
[1] all
Levels: all
> Idents(sce) <- 'all'





install.packages('dplyr')
library(dplyr)
library(Seurat)

