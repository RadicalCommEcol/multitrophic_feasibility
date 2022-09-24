
boot_prob_excl_Omega_raw <- function(A_int, number_Omega_replicates,
                                 number_boot_replicates){
  
  
  simple_mean <- function(x, indices){
    return(sum(x[indices])/length(indices))
  }
  
  
  dimensions <- ncol(A_int)
  
  incenter_inradius_isoprob <- incenter_inradius_isoprob_calculation(A_int)
  
  I <- incenter_inradius_isoprob[[1]]
  
  center <- I
  
  name_colums <- c("species", paste0("prob_excl_rep_",1:number_boot_replicates),
                   paste0("Omega_rep_",1:number_boot_replicates))
  
  caracoles_prob_exclusion_plants_mat <- matrix(rep(0, dimensions*(1 + 2*number_boot_replicates)),
                                                nrow = dimensions, 
                                                ncol = (1 + 2*number_boot_replicates))
  
  bootstrap_raw_data <-
    exclusion_probabilities_par_bootstrap(A_int, replicates = number_Omega_replicates)
  
  for (sp_i in 1:dimensions) {
    
    bootmean_area_excl_sp_i <- boot(data = bootstrap_raw_data[sp_i,], 
                                    statistic = simple_mean, 
                                    R = number_boot_replicates)
    
    caracoles_prob_exclusion_plants_mat[sp_i, 2:(1+number_boot_replicates)] <- 
      bootmean_area_excl_sp_i[["t"]] 
    
  }
  
  Omega_values <- colSums(caracoles_prob_exclusion_plants_mat[, c(2:(1+number_boot_replicates))])
  
  caracoles_prob_exclusion_plants_mat[, c((2+number_boot_replicates):(1 + 2*number_boot_replicates))] <-
    matrix(rep(Omega_values,dimensions), nrow = number_boot_replicates) %>% t()
  
  for (sp_i in 1:nrow(caracoles_prob_exclusion_plants_mat)) {
    caracoles_prob_exclusion_plants_mat[sp_i, c(2:(1+number_boot_replicates))] <-
      caracoles_prob_exclusion_plants_mat[sp_i, c(2:(1+number_boot_replicates))]/
      caracoles_prob_exclusion_plants_mat[sp_i, c((2+number_boot_replicates):(1 + 2*number_boot_replicates))]
  }
  
  caracoles_prob_exclusion_plants_aux <- as_tibble(data.frame(caracoles_prob_exclusion_plants_mat))
  colnames(caracoles_prob_exclusion_plants_aux) <- name_colums

  caracoles_prob_exclusion_plants_aux$species = colnames(A_int)
  
  return(caracoles_prob_exclusion_plants_aux)
  
}
