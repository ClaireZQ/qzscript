
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




