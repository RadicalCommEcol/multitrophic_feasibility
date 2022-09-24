# verteces_unit_ball computes the vertices of the simplicial cone when it intersects the unit ball
verteces_unit_ball <- function(A_int){
  
  dimensions <- ncol(A_int)
  
  num_species <- 1:dimensions
  
  r_values <- data.frame(matrix(ncol = (2*dimensions), nrow = dimensions))
  r_names <- c(paste0("N_", num_species),paste0("r_", num_species)) 
  colnames(r_values) <- r_names
  r_values$color_sol <- NA
  
  for(i in 1:dimensions){
    
    sup_A_col_i <- A_int[,i]
    
    Nv_i_aux <- sqrt(1/sum( sup_A_col_i* sup_A_col_i))
    
    # We add a zero in x_aux at position i
    if(i==1){
      Nv_i <- c(Nv_i_aux,rep(0,dimensions-1))
    }else if(i==dimensions){
      Nv_i <- c(rep(0,dimensions-1),Nv_i_aux)
    }else{
      Nv_i <- c(rep(0,i-1),Nv_i_aux,rep(0,dimensions-i))
    }
    
    Nv_i_mat <- matrix(Nv_i,nrow = dimensions, ncol = 1)
    rv_i_mat <- (-1)*A_int %*%  Nv_i_mat
    rv_i <- as.numeric(rv_i_mat)
    
    r_values[i,1:dimensions] <- Nv_i
    r_values[i,(dimensions+1):(2*dimensions)] <- rv_i
    r_values$color_sol[i] <- paste0("vertex ",i)
    
  }
  
  return(r_values)
  
}