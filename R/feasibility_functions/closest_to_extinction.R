
closest_to_extinction <- function(center, A_int){
  
  dimensions <- ncol(A_int)
  curves_data <- cone_verteces_director_vertices(A_int)
  
  dir_ver_mat <- as.matrix(curves_data[,(dimensions+3):ncol(curves_data)])
  
  # we check the direction of the director vectors
  correction_direction <- NULL
  
  for(i in 1: dimensions){
    
    if(all((dir_ver_mat %*% as.numeric(curves_data[i,(1:dimensions)]) %>% round(12))>=0)){
      correction_direction <- c(correction_direction,1)
    }else{
      correction_direction <- c(correction_direction,-1)
    }
    
  }
  
  acr_distance_center_side <- NULL
  
  for(i in 1:nrow(curves_data)){
    
    normal_vector_i <- correction_direction[i] * curves_data[i,c((dimensions+3):ncol(curves_data))]
    sin_theta <- sum(normal_vector_i*center)
    acr_distance_center_side <- c(acr_distance_center_side,asin(sin_theta))
    
  }
  
  min_arc_dist <- min(acr_distance_center_side)
  vertex_min_arc_dist_index <- which(round(acr_distance_center_side,12) == round(min_arc_dist,12))
  
  tangent_points_matrix <- matrix(rep(0,dimensions*length(vertex_min_arc_dist_index)),
                                  nrow = length(vertex_min_arc_dist_index),
                                  ncol = dimensions)
  
  colnames(tangent_points_matrix) <- paste0("r_T_",1:dimensions)
  
  for(i in vertex_min_arc_dist_index){
    
    normal_vector_i <- correction_direction[i] * curves_data[i,c((dimensions+3):ncol(curves_data))]
    sin_theta <- sum(normal_vector_i*center)
    tangent_aux <- center - sin_theta*normal_vector_i %>% as.numeric()
    module_tangent_aux <- sqrt(sum(tangent_aux*tangent_aux))
    tangent_points_matrix[i,] <- tangent_aux / module_tangent_aux
    
  }
  
  N_tangent_points_matrix <- matrix(rep(0,dimensions*length(vertex_min_arc_dist_index)),
                                    nrow = length(vertex_min_arc_dist_index),
                                    ncol = dimensions)
  
  colnames(N_tangent_points_matrix) <- paste0("N_T_",1:dimensions)
  
  for(i in 1:nrow(N_tangent_points_matrix)){
    
    r_T_i <- tangent_points_matrix[i,] %>% as.numeric()
    r_T_i_column <- matrix(r_T_i, nrow = dimensions, ncol = 1)
    N_T_i_colum <- (-1)*pracma::inv(A_int) %*% r_T_i_column %>% as.numeric()
    N_tangent_points_matrix[i,] <- N_T_i_colum
    
  }
  
  return(cbind(tangent_points_matrix,N_tangent_points_matrix))
  
}
