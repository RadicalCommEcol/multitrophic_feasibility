
small_omega_iso <- function(A_int){
  
  # coordinates of FD's Vertices
  verteces_unit_ball(A_int)
  
  # Get main curves data: vertices and normal vectors for their hyperplanes
  dimensions <- ncol(A_int)
  
  # Estimation of the incenter position and of the inradius for all the three vertices (of the 
  # simplicial cone) combinations
  incenter_inradius_isoprob <- incenter_inradius_isoprob_calculation(A_int)

  isoprobability <- incenter_inradius_isoprob[[3]]
  
  return(isoprobability^(1/nrow(A_int)))
  
}

