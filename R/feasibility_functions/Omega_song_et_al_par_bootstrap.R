Omega_song_et_al_par_bootstrap <- function(A_int, replicates = 1e3, parallelize = FALSE) {
  
  Sigma <- solve(t(A_int) %*% A_int,tol = 1e-17)
  S <- ncol(A_int)
  
  # m <- matrix(0, S, 1)
  # a <- matrix(0, S, 1)
  # b <- matrix(Inf, S, 1)
  
  #set.seed(123)
  
  if(parallelize){
  
  res_list <- foreach (i = 1:replicates, .combine=c) %dopar% {
    
    d <- try(mvtnorm::pmvnorm(lower = rep(0, S),
                              upper = rep(Inf, S),
                              mean = rep(0, S), sigma = Sigma),silent = TRUE)
    out <- ifelse(class(d) == "try-error",0,d[1])
    
    # print(out)

  }
  }else{
    res_list <- foreach (i = 1:replicates, .combine=c) %do% {
      
      d <- try(mvtnorm::pmvnorm(lower = rep(0, S),
                                upper = rep(Inf, S),
                                mean = rep(0, S), sigma = Sigma),silent = TRUE)
      out <- ifelse(class(d) == "try-error",0,d[1])
      
      # print(out)
      
    }
  }

  return(res_list)
  
}


