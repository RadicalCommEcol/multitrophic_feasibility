
isotropy_metrics_4_int_matrix <- function(A_int, number_Omega_replicates,
                                          number_boot_replicates){
  
  prob_excl_Omega_df <- boot_prob_excl_Omega_raw(A_int, number_Omega_replicates,
                                             number_boot_replicates)
  
  isotropy_agg_info <- tibble(interaction_matrix = 1)
  small_omega_mean <- NA
  small_omega_lowerCI <- NA
  small_omega_upperCI <- NA
  isotropy_index_mean <- NA
  isotropy_index_lowerCI <- NA
  isotropy_index_upperCI <- NA
  PV_index_mean <- NA
  PV_index_lowerCI <- NA
  PV_index_upperCI <- NA
  
  
  Omega_rep <- prob_excl_Omega_df[1,(1+number_boot_replicates):(2*number_boot_replicates)] %>%
    as.numeric()
  small_omega_rep <- Omega_rep^(1/nrow(A_int))
  isotropy_rep <- NULL
  PV_rep <- NULL
  
  for(i in 1:number_boot_replicates){
    
    probability_vector <- prob_excl_Omega_df[,i+1] %>% pull()
    
    isotropy_rep <-  c(isotropy_rep, anisotropy_index(probability_vector))
    PV_rep <- c(PV_rep, PV(probability_vector))
    
  }
  
  small_omega_rep_CI <- quantile(small_omega_rep, prob=c(.025,.975)) %>% as.numeric()
  isotropy_rep_CI <- quantile(isotropy_rep, prob=c(.025,.975)) %>% as.numeric()
  PV_rep_CI <- quantile(PV_rep, prob=c(.025,.975)) %>% as.numeric()
  
  isotropy_agg_info$small_omega_mean <- mean(small_omega_rep, na.rm = T)
  isotropy_agg_info$small_omega_lowerCI <- small_omega_rep_CI[1]
  isotropy_agg_info$small_omega_upperCI <- small_omega_rep_CI[2]
  isotropy_agg_info$isotropy_index_mean <- mean(isotropy_rep, na.rm = T)
  isotropy_agg_info$isotropy_index_lowerCI <- isotropy_rep_CI[1]
  isotropy_agg_info$isotropy_index_upperCI <- isotropy_rep_CI[2]
  isotropy_agg_info$PV_index_mean <- mean(PV_rep, na.rm = T)
  isotropy_agg_info$PV_index_lowerCI <- PV_rep_CI[1]
  isotropy_agg_info$PV_index_upperCI <- PV_rep_CI[2]
  
  return(isotropy_agg_info %>% select(-interaction_matrix))
  
}
