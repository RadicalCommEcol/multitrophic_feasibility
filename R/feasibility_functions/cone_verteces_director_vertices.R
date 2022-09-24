# cone_verteces_director_vertices computes the vertices of the simplicial cone when it 
# intersects the unit ball, as well as their  normal vectors.
cone_verteces_director_vertices <- function(A_int){
  
  dimensions <- ncol(A_int)
  all_verteces_data <- verteces_unit_ball(A_int) %>% mutate(V_name=NA)
  all_verteces_data$V_name <- 1:dimensions
  all_verteces_r <- all_verteces_data[,(dimensions+1):ncol(all_verteces_data)]
  
  new_columns_names <- paste0("dir_vec_",1:dimensions)
  all_verteces_r[,new_columns_names] <- NA
  
  for(i in 1:dimensions){
    
    #other_vertex_data <<- all_verteces_r[-i,1:dimensions]
    
    other_vertex_data <- all_verteces_r[-i,1:dimensions]
    #print(other_vertex_data)
    # dir_vect_initialization <- rep(1,dimensions)#runif(dimensions)
    # director_vector_i_solver <- nleqslv(dir_vect_initialization,director_vector_hyperplane_solver)
    # director_vector_i <- director_vector_i_solver[["x"]]
    #print(director_vector_i)
    
    director_vector_i_aux <- pracma::crossn(as.matrix(other_vertex_data))
    module_director_vector_i_aux <- sqrt(sum(director_vector_i_aux*director_vector_i_aux))
    director_vector_i <- director_vector_i_aux/module_director_vector_i_aux
    all_verteces_r[i,new_columns_names] <- director_vector_i
    
  }
  
  return(all_verteces_r)
  
}
