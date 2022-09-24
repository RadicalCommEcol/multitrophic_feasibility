
further_away_species <- function(center, A_int){
  
  dimensions <- ncol(A_int)
  curves_data <- cone_verteces_director_vertices(A_int)
  acr_distance_center_vertex <- NULL
  for(i in 1:nrow(curves_data)){
    
    vertex_i <- curves_data[i,c(1:dimensions)]
    acr_distance_center_vertex <- c(acr_distance_center_vertex,acos(sum(vertex_i*center)))
    
  }
  
  max_arc_dist <- max(acr_distance_center_vertex)
  vertex_max_arc_dist_index <- which(acr_distance_center_vertex == max_arc_dist)
  # vertex_max_arc_dist <- curves_data[vertex_max_arc_dist_index,c(1:dimensions)] %>% as.numeric()
  # population_max_arc_dist <- (-1) * pracma::inv(A_int) %*% matrix(vertex_max_arc_dist, nrow = dimensions, ncol = 1)
  
  return(vertex_max_arc_dist_index)
  
}
