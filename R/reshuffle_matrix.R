# this function reshuffles matrix elements keeping its names
# it does so differently depending on whether the matrix is square or not.
# if square, it reshuffles topology, assigning links to any matrix position
# In this case, diagonal elements are kept but their strenghts are also reshuffled.
# if bipartite, it simply reshuffles everything but makes sure all columns
# have at least one link if this is met in the original matrix.
reshuffle_matrix <- function(A){
  
  randA <- A
  
  if(identical(colnames(A),rownames(A))){
    
    randA[randA != 0] <- 0
    
    # how many links in the non diagonal
    non.diag.links <- sum(A != 0) - nrow(A)
    
    # this is the vector of interactions strenghts reshuffled
    randomized.strengths <- sample(A[A != 0])
    
    vec.upper.tri <- as.vector(upper.tri(A))
    vec.lower.tri <- as.vector(lower.tri(A))
    # vector of non diagonal positions
    vec.non.diag <- vec.upper.tri|vec.lower.tri
    
    # randomly sample the number of links in the non diagonal positions
    non.diag.pos <- sample(which(vec.non.diag == T),non.diag.links)
    
    randA[non.diag.pos] <- randomized.strengths[1:non.diag.links]
    diag(randA) <- randomized.strengths[(non.diag.links+1):length(randomized.strengths)]
    # how many links in the upper triangle
    # randupper <- sum(A[upper.tri(A)] != 0)
    # 
    # # diagonal positions have links
    # diag.link.num <- nrow(A)
    # 
    # # randomly sample positions from the upper triangle
    # vec.upper.tri <- as.vector(upper.tri(A))
    # # conversion works
    # # matrix(vec.upper.tri,nrow = diag.link.num)
    # upper.tri.pos <- sample(which(vec.upper.tri == T),randupper)
    # 
    # # this is the vector of interactions strenghts reshuffled
    # randomized.strengths <- sample(A[A != 0])
    # 
    # # assign the first reshuffled strengths to the upper triangle
    # # this works because a matrix is a vector
    # randA[upper.tri.pos] <- randomized.strengths[1:randupper]
    # 
    # # get the positions of the lower triangle
    # randlower <- t(randA) #- diag(randA)
    # # assign diagonal strengths
    # diag(randA) <- randomized.strengths[(randupper+1):(randupper+diag.link.num)]
    # # assign the remaining interaction strengths
    # randA[which(randlower != 0)] <- 
    #   randomized.strengths[(randupper+diag.link.num+1):length(randomized.strengths)]
  }else{
    # randA <- A
    
    randA[] <- sample(A)  
    
    # make sure that every column (insect sp) is assigned an interaction
    # but only if the original data allows it (i.e all columns have at least one)
    if(sum(colSums(A) == 0)){
      
      empty.cols <- sum(colSums(randA) == 0)
      count <- 0
      while(empty.cols > 0 | count < 1e3){
        randA[] <- sample(A)  
        empty.cols <- sum(colSums(randA) == 0)
        count <- count+1
      }# while
    }# if
    # while()
    
  }
  
  return(randA)
}