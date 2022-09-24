
tangent_points_spherical_triangle <- function(A_int,Vi,Vj,Vk){
  
  dimensions <- ncol(A_int)
  
  incenter_inradius <- incenter_inradius_isoprob_calculation(A_int)
  I <- incenter_inradius[[1]]
  inradius <- incenter_inradius[[2]]
  cos_theta <- cos(inradius)
  
  set.seed(123)
  T0 <- runif(3+3*dimensions)
  incenter_evaluation <-  nleqslv(T0,incenter_tangent_points_solver,
                                  control=list(maxit = 1e8, xtol=1e-10, allowSingular =T))
  incenter_data <- incenter_evaluation[["x"]]
  incenter_data # Here we obtain the results before the iteration that makes them zero. These data is
  #useful to get coordinates of the tangent points.
  incenter_evaluation[["termcd"]]
  
  T1 <- incenter_data[(4+0*dimensions):(3+1*dimensions)]
  T2 <- incenter_data[(4+1*dimensions):(3+2*dimensions)]
  T3 <- incenter_data[(4+2*dimensions):(3+3*dimensions)]
  
  return(list(T1,T2,T3))
  
}

