
feasibility_metrics <- function(A,
                                omega.replicates,
                                bootstrap.replicates,
                                noise.replicates,
                                noise.threshold = c(0,1),
                                noise.type = c("random","targeted")){
  
  # is A invertible?
  f <- function(m) class(try(solve(t(m) %*% m),
                             silent = T))[[1]] == "matrix"
  # the operation fails, A is not invertible
  if (f(A) == FALSE) {
    
    # add noise to the matrix, replicated bootstrap.replicates times
    perturbed_matrices <- replicate(n = noise.replicates,
                                    expr = MultitrophicFun::add_noise_matrix(A,
                                                            noise.threshold = noise.threshold,
                                                            noise.type = noise.type),
                                    simplify = FALSE)
    res.noise.list <- list()
    
    for(i.noise.rep in 1:noise.replicates){
      my.res <- isotropy_metrics_4_int_matrix(A_int = perturbed_matrices[[i.noise.rep]],
                                              number_Omega_replicates = omega.replicates,
                                              number_boot_replicates = bootstrap.replicates)
      my.res$omega_isotropic <- small_omega_iso(A_int = perturbed_matrices[[i.noise.rep]])
      my.res$noise.replicate <- i.noise.rep
      
      res.noise.list[[i.noise.rep]] <- my.res
      
    }
    
    res.noise.df <- bind_rows(res.noise.list)
    
    # TODO check that error propagation works this or other way
    # this is likely too simple
    res <- res.noise.df %>%
      summarise(small_omega_mean = mean(small_omega_mean,na.rm = TRUE),
                small_omega_lowerCI = mean(small_omega_lowerCI,na.rm = TRUE),
                small_omega_upperCI = mean(small_omega_upperCI,na.rm = TRUE),
                omega_isotropic = mean(omega_isotropic,na.rm = TRUE),
                isotropy_index_mean = mean(isotropy_index_mean,na.rm = TRUE),
                isotropy_index_lowerCI = mean(isotropy_index_lowerCI,na.rm = TRUE),
                isotropy_index_upperCI = mean(isotropy_index_upperCI,na.rm = TRUE))
    
  }else{
    
    res <- isotropy_metrics_4_int_matrix(A_int = A,
                                         number_Omega_replicates = omega.replicates,
                                         number_boot_replicates = bootstrap.replicates)
    res$omega_isotropic <- small_omega_iso(A_int = A)
    
  }
  
  res.clean <- res[,c("small_omega_mean","small_omega_lowerCI","small_omega_upperCI","omega_isotropic",
                           "isotropy_index_mean","isotropy_index_lowerCI","isotropy_index_upperCI")]
  names(res.clean)[1:3] <- c("omega_mean","omega_lowerCI","omega_upperCI")

  # -------------------------------------------------------------------------
  return(res.clean)
}

