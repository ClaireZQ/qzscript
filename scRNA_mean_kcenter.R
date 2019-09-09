####use kmeans to compare similarity of two matrixes 
mat1 <- matrix(NA,nr = 9,nc = 428)
colnames(mat1) <- old@assays$RNA@var.features
rownames(mat1) <- paste0('cluster',0:8)
for(i in 0:8){
	names <- rownames(old@meta.data)[old@meta.data$seurat_clusters == i]
	a <- old@assays$RNA@scale.data[old@assays$RNA@var.features,names]	
	mat1[i+1,] <- apply(a,1,mean)
}
mat1 <- matrix(NA,9,2)
colnames(mat1) <- c("B1",'B2')
rownames(mat1) <- paste0('cl',0:8)
for (i in 0:9) {
	names1 <- rownames(P12@meta.data)[P12@meta.data$seurat_clusters == i & P12@meta.data$sample == 'P1']
	counts1 <- length(names1)
    names2 <- rownames(P12@meta.data)[P12@meta.data$seurat_clusters == i & P12@meta.data$sample == 'P2']
    counts2 <- length(names2)
    mat[i+1,]<-c(counts1,counts2)
} 





mat2 <- matrix(NA,nr = 8,nc = 428)
colnames(mat2) <- uterus@assays$RNA@var.features
rownames(mat2) <- paste0('cluster',0:7)
for(j in 0:7){
	names <- rownames(uterus@meta.data)[uterus@meta.data$seurat_clusters == j]
	b <- uterus@assays$RNA@scale.data[uterus@assays$RNA@var.features,names]	
	mat2[j+1,] <- apply(b,1,mean)
}
int <- intersect(colnames(mat1),colnames(mat2))
tmp1 <- mat1[,int]
tmp2 <- mat2[,int]
# a <- tmp1[rownames(tmp1)=='clusteri',]字符串中的i是不变的,只是可以提取以行来填充的矩阵
# b <- tmp2[rownames(tmp2)=='clusterj',]
final_mat1 <- matrix(NA,15,12)
rownames(final_mat1) <- paste0('cluster',0:14)
colnames(final_mat1) <- paste0('cluster',0:11)
 for (i in 0:14) {
  for (j in 0:11) {
   a <- tmp1[i+1,]
   b <- tmp2[j+1,]
   sce1 <- as.vector(a)
   sce2 <- as.vector(b)
   len1 <- length(sce1)
   len2 <- length(sce2)
   numer = sum(sce1*sce2)
   denom = sqrt(sum(sce1**2)*sum(sce2**2))
   similar = numer / denom
   re = (similar+1)/2
   final_mat[i+1,j+1] <- re 
  }      
 }

for (i in 0:14) {
  for (j in 0:11) {
   a <- tmp1[i+1,]
   b <- tmp2[j+1,]
   sce1 <- as.vector(a)
   sce2 <- as.vector(b)
   len1 <- length(sce1)
   len2 <- length(sce2)
   numer = sum(sce1*sce2)
   denom = sqrt(sum(sce1**2)*sum(sce2**2))
   similar = numer / denom
   final_mat1[i+1,j+1] <- similar
  }      
 }


argene <- union(colnames(mat1),colnames(mat2))
createVector <- function(x,y,argene){
	vectorx <- matrix(0,nr=1,nc=length(argene)); colnames(vectorx) <- argene 
	vectory <- matrix(0,nr=1,nc=length(argene)); colnames(vectory) <- argene
	for( i in colnames(x)){
		print(i)
		vectorx[1,i] <- x[1,i]
	}
	for( j in colnames(y)){
		vectory[1,j] <- y[1,j]
	}
	
	tmp <- list(vectorx,vectory)
	return (tmp)	

}


x <- matrix(NA, nr = 9, nc = length(argene))
y <- matrix(NA, nr = 8, nc = length(argene))
result <- matrix(NA, 9, 8)

for(i in 0:8){
	for(j in 0:7){
		a <- createVector(mat1[i+1,],mat2[j+1,],argene)
		x1 <- a[[1]]
		y1 <- a[[2]]
		x[i+1,] <- x1
		y[j+1,] <- y1
		sce1 <- as.vector(x)
		sce2 <- as.vector(y)
		len1 <- length(sce1)
        len2 <- length(sce2)
        numer = sum(sce1*sce2)
        denom = sqrt(sum(sce1**2)*sum(sce2**2))
        similar = numer / denom
		result[i+1,j+1] <- similar 

	}
}
 






















