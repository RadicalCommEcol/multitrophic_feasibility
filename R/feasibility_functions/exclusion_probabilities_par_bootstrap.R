
exclusion_probabilities_par_bootstrap <- function(A_int, replicates = 1e3){
  
  dimensions <- ncol(A_int)
  
  incenter_inradius_isoprob <- incenter_inradius_isoprob_calculation(A_int)
  
  I <- incenter_inradius_isoprob[[1]]
  
  all_verteces_data_aux <- verteces_unit_ball(A_int)
  all_verteces_data <- all_verteces_data_aux[,c((dimensions+1):(2*dimensions))]
  
  feasibility_parts <- matrix(rep(0,nrow(A_int)*replicates),
                              nrow = nrow(A_int), ncol = replicates)
  
  names_sp <- colnames(A_int)
  
  for (i in 1:nrow(A_int)) {
    
    cat(names_sp[i],"\n")
    
    A_int_mod <- A_int
    
    A_int_mod[,i] <- -I
    
    feasibility_parts[i, ] <- Omega_song_et_al_par_bootstrap(A_int_mod, 
                                                             replicates,
                                                             parallelize = FALSE)
    
  }
  
  return(feasibility_parts)
  
}
