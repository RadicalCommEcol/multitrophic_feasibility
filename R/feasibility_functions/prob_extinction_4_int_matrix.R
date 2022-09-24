prob_extinction_4_int_matrix <- function(A_int, number_Omega_replicates,
                                          number_boot_replicates){
  
  prob_excl_Omega_df <- boot_prob_excl_Omega_raw(A_int, number_Omega_replicates,
                                                 number_boot_replicates)
  
  prob_excl_agg_info <- tibble(species = colnames(A_int))
  prob_excl_agg_info$prob_excl_mean <- NA
  prob_excl_agg_info$prob_excl_lowerCI <- NA
  prob_excl_agg_info$prob_excl_upperCI <- NA
  
  prob_excl_df <- prob_excl_Omega_df[,2:(number_boot_replicates-1)]
  
  for(i in 1:nrow(prob_excl_df)){
    
    probability_excl_Sp <- prob_excl_df[i, ] %>% as.numeric()
    
    probability_excl_CI <- quantile(probability_excl_Sp, prob=c(.025,.975)) %>% as.numeric()
    
    prob_excl_agg_info$prob_excl_mean[i] <- mean(probability_excl_Sp, na.rm = T)
    prob_excl_agg_info$prob_excl_lowerCI[i] <- probability_excl_CI[1]
    prob_excl_agg_info$prob_excl_upperCI[i] <- probability_excl_CI[2]
    
  }
  
  
  return(prob_excl_agg_info)
  
}
