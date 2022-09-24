interaction_overlap <- function(A, method = "horn"){
  
  dis.matrix <- as.matrix(vegan::vegdist(A,method = method))
  overlap.matrix <- 1 - dis.matrix
  
  sp.names <- as.character(1:nrow(A))
  if(!is.null(rownames(A))){
    sp.names <- rownames(A)
  }
  
  overlap.pairs <- list()
  
  for(i in 1:nrow(overlap.matrix)){
    for(j in 1:ncol(overlap.matrix)){
      overlap.pairs[[length(overlap.pairs)+1]] <-
        data.frame(sp1 = sp.names[i],sp2 = sp.names[j], overlap = overlap.matrix[i,j])
    }
  }
  overlap.df <- dplyr::bind_rows(overlap.pairs)
  return(overlap.df)
}