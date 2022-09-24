
exclusion_probabilities <- function(A,
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
      my.res <- prob_extinction_4_int_matrix(A_int = perturbed_matrices[[i.noise.rep]],
                                              number_Omega_replicates = omega.replicates,
                                              number_boot_replicates = bootstrap.replicates)

      res.noise.list[[i.noise.rep]] <- my.res
      
    }
    
    res.noise.df <- bind_rows(res.noise.list)
    
    # TODO check that error propagation works this or other way
    # this is likely too simple
    res <- res.noise.df %>% group_by(species) %>%
      summarise(prob_excl_mean = mean(prob_excl_mean),
                prob_excl_lowerCI = mean(prob_excl_lowerCI),
                prob_excl_upperCI = mean(prob_excl_upperCI))
    
  }else{
    
    res <- prob_extinction_4_int_matrix(A_int = A,
                                         number_Omega_replicates = omega.replicates,
                                         number_boot_replicates = bootstrap.replicates)
    
  }
  
  # -------------------------------------------------------------------------
  return(res)
}

