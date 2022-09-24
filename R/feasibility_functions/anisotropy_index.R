anisotropy_index <- function(probability_vector){
  
  res <- NA
  
  # check the input parameter: probability_vector
  if(any(is.na(probability_vector))){
    cat("Probabilities should not contain NAs.\n")
  }else if(!all(class(probability_vector) %in% c("numeric"))){
    cat("Probabilities should be numbers.\n")
  }else if(all(class(probability_vector) %in% c("matrix","array"))){
    cat("The input parameter should be vector.\n")
  }else if(!all(probability_vector >=0)){
    cat("Probabilities should be non-negaive numbers.\n")
  }else if(round(sum(probability_vector),2) != 1){
    cat("The sum of all probabilities should be equal to one.\n")
  }else{
    
    # If the vector with probabilities is OK, then the index is calculated.

    if(length(probability_vector)>1){
    number_species <- length(probability_vector)
    
    probability_vector_positive <- probability_vector[probability_vector > 0]
    
    log_prob <- log(probability_vector_positive)
    
    res <- -sum(probability_vector_positive*log_prob)/log(number_species)
    }
  }# if-else
  return(res)
}