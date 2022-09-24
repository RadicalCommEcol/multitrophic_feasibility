
# Given 3 vertices on a (d-1)-unit ball (Vi, Vj,... Vn), this function estimates the following info_
# 1) their incenter;
# 2) the great circle distance from I to the incircle in the spherical polygon defined by Vi,
# Vj... Vk;
# the area of the spherical/area of the sphere

incenter_inradius_isoprob_calculation <- function(A_int){
  
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

  inv_dir_ver_mat <- pracma::inv(dir_ver_mat)
  I_aux_mat <- inv_dir_ver_mat %*% matrix(correction_direction,nrow = dimensions, ncol = 1)
  I_aux <- as.numeric(I_aux_mat)
  sin_theta <- sqrt(1/sum(I_aux*I_aux)) # After imposing normalization
  # cos_theta <- sqrt(1-sin_theta*sin_theta)
  I_aux2 <- I_aux*sin_theta
  I <- I_aux2/sqrt(sum(I_aux2*I_aux2))
  
  dir_ver_mat_zero <- as.matrix(curves_data[,(dimensions+3):ncol(curves_data)])
  
  dir_ver_mat_zero %*% I
  
  # probability of the isotropic area: Area of the spherical cap with angle theta/ Area of full circle
  a_par <- 0.5*(dimensions-1)
  b_par <- 0.5
  probability_sp_cap <- 0.5*zipfR::Rbeta(sin_theta*sin_theta,a_par,b_par)
  
  return(list(I,asin(sin_theta),probability_sp_cap))
  
}
