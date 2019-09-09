
matrix similar comparison
float <- mtx_similar(scRNA1:np.ndarray, scRNA2:np.ndarray)
farr1 = scRNA1.ravel()
farr2 = scRNA2.ravel()
len1 = length(farr1)
len2 = length(farr2)

if(len1 > len2){
      farr1 = farr1[:len2]
else:
      farr2 = farr2[:len1]


numer = np.sum(farr1 * farr2)
denom = np.sqrt(np.sum(farr1**2) * np.sum(farr2*2))
similar = numer/denom
(similar + 1)/2



load("liu")
sce1 <- uterus
load('zhang')
sce2 <- uterus
sce1 <- as.vector(sce1)
sce2 <- as.vector(sce2)
len1 <- length(sce1)
len2 <- length(sce2)
if (len1 > len2)
{
	sce1 <- sce1[0:len2]
} else {
	sce2 <- sce2[0:len1]
}

numer = sum(sce1*sce2)
denom = sqrt(sum(sce1**2)*sum(sce2**2))
similar = numer / denom
re = (similar+1)/2

clusters <- uterus@meta.data$seurat_clusters == 0
clusters0 <- rownames(uterus@meta.data)[clusters]
clusters00 <- uterus@assays$RNA@scale.data[,clusters0]
 
clusters1 <- rownames(uterus@meta.data)[uterus@meta.data$seurat_clusters==1]
clusters01 <- uterus@assays$RNA@scale.data[,clusters1]
old <- uterus
mat <- matrix(NA,15,12)
for(j in 0:14)
{	clustersj <- rownames(old@meta.data)[old@meta.data$seurat_clusters == j]
  	clustersj <- old@assays$RNA@scale.data[,clustersj]
  	sce1 <- as.vector(clustersj)
  	len1 <- length(sce1)
for(i in 0:11){     
	clustersi <- rownames(uterus@meta.data)[uterus@meta.data$seurat_clusters == i]
	clustersi <- uterus@assays$RNA@scale.data[,clustersi]
	sce2 <- as.vector(clustersi)
	len2 <- length(sce2)
	if (len1 > len2)
{
	sce1 <- sce1[0:len2]
} else {
	sce2 <- sce2[0:len1]
}

numer = sum(sce1*sce2)
denom = sqrt(sum(sce1**2)*sum(sce2**2))
similar = numer / denom
re = (similar+1)/2
mat[(j+1),(i+1)] <- re
}
}



if(sce1.shape!=sce2.shape){
	minx = min(sce1.shape[0],sce2.shape[0])
	miny = min(sce1.shape[1],sce2.shape[1])
	differ = min1[0:minx,0:miny]-min2[0:minx,0:miny]
} else {differ = sce1 - sce2
	numera = sum(differ**2)
	denom=sum(sce1**2)
	similar = 1-(numera/denom)
	re = similar

}






























