
# CreateTmpMatrix -------------------------------------
#'                                                    -
#' Arguments                                          -
#' @ x: expression matrix with Cells * Genes          -
#' @ Y: another expression matrix with Cells * Genes  -
#'                                                    -
# -----------------------------------------------------
CreateTmpMatrix <- function(x, y){
  tmp <- matrix(0, nr = length(union(rownames(x),rownames(y))), nc = length(union(colnames(x),colnames(y))))
  rown <- union(rownames(x), rownames(y))
  coln <- union(colnames(x), colnames(y))
  rownames(tmp) <- rown; colnames(tmp) <- coln
  return (tmp)
}

# CompletionMatrix ------------------------------------
#'                                                    -
#' Arguments                                          -
#' @ tmp: CreateTmpMatrix(x, y)                       -
#' @ x: expression matrix with Cells * Genes          -
#'                                                    -
# -----------------------------------------------------
CompletionMatrix <- function(tmp, x){
  for(i in rownames(x)){
    for(j in colnames(x)){
      tmp[i,j] <- x[i,j]
    }
  }
  return (tmp)
}




tmp <- CreateTmpMatrix(cll0,clu0)
cll0_final <- CompletionMatrix(tmp,cll0)
clu0_final <- CompletionMatrix(tmp,clu0)

 #extract 428genes 组成矩阵
a <- cll0[old@assays$RNA@var.features,]
b <- clu0[uterus@assays$RNA@var.features,]

a[[i]],b[[j]]



a <- list()
for(i in 0:14){
  names <- rownames(old@meta.data)[old@meta.data$seurat_clusters == i]
  a[[i+1]] <- old@assays$RNA@scale.data[old@assays$RNA@var.features,names]
}

b <- list()
for (j in 0:11) {
   names <- rownames(uterus@meta.data)[uterus@meta.data$seurat_clusters == j]
   b[[j+1]] <- uterus@assays$RNA@scale.data[uterus@assays$RNA@var.features,names]
}
a <- list()
for(i in 0:14){
  names <- rownames(old@meta.data)[old@meta.data$seurat_clusters == i]
  a[[i+1]] <- old@assays$RNA@counts[old@assays$RNA@var.features,names]
}

b <- list()
for (j in 0:11) {
   names <- rownames(uterus@meta.data)[uterus@meta.data$seurat_clusters == j]
   b[[j+1]] <- uterus@assays$RNA@counts[uterus@assays$RNA@var.features,names]
}

mat <- matrix(NA,15,12)

for(i in 0:14){
  for(j in 0:11){
     tmp <- CreateTmpMatrix(a[[i+1]],b[[j+1]])
     final_a <- CompletionMatrix(tmp,a[[i+1]])
     final_b <- CompletionMatrix(tmp,b[[j+1]])
     sce1 <- as.vector(final_a)
     sce2 <- as.vector(final_b)
     len1 <- length(sce1)
     len2 <- length(sce2)
     

     numer = sum(sce1*sce2)
     denom = sqrt(sum(sce1**2)*sum(sce2**2))
     similar = numer / denom
     re = (similar+1)/2
     mat[i+1,j+1] <- re 



  }
}


SimilarNorm <- function(final_a,final_b){
  differ = final_a-final_b
  dist = norm(differ, type = 'F')
  len1 = norm(final_a)
  len2 = norm(final_b)
  denom = (len1 + len2) / 2
  similar = 1 - (dist/denom)
  return(similar)
}

fab <- function(final_a,final_b){
  x1 <- final_a - final_b
  x2 <- t(final_a - final_b)
  x3 <- x1%*%x2
  x4 <- sum(as.vector(x3))
  x5 <- sqrt(x4)
  return(x5)
}
for(i in 0:14){
  for(j in 0:11){
     tmp <- CreateTmpMatrix(a[[i+1]],b[[j+1]])
     final_a <- CompletionMatrix(tmp,a[[i+1]])
     final_b <- CompletionMatrix(tmp,b[[j+1]])
     
mat[i+1,j+1] <- SimilarNorm(final_a,final_b)
}
}






fab <- function(final_a,final_b){
  x1 <- final_a - final_b
  x2 <- t(final_a - final_b)
  x3 <- x1%*%x2
  x4 <- sum(as.vector(x3))
  x5 <- sqrt(x4)
  return(x5)
}










a <- list()
for(i in 0:1){
  names <- rownames(old@meta.data)[old@meta.data$seurat_clusters == i]
  a[[i+1]] <- old@assays$RNA@scale.data[old@assays$RNA@var.features,names]
}

b <- list()
j=0
   names <- rownames(uterus@meta.data)[uterus@meta.data$seurat_clusters == j]
   b[[i+1]] <- uterus@assays$RNA@scale.data[uterus@assays$RNA@var.features,names]






for(i in 1:6){
  write.xlsx(a[[i]],file = 'a1.xlsx',
    sheetName = names(a)[i],
    col.names = T, row.names = T,
    append = T)
}



names(a) <- pate0('matrix', 1:10)


    tmp <- CreateTmpMatrix(a[[1]],b[[1]])
     final_a <- CompletionMatrix(tmp,a[[1]])
     final_b <- CompletionMatrix(tmp,b[[1]])

