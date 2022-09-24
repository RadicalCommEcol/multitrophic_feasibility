
#' combine interaction matrices in a single block matrix
#' 
#' In the process, generate intraguild matrices for herbivore and pollinator guilds.
#' This needs additional info: phenology, nesting, and larval food requirements.
#' Note that plant phenology and animal phenology is in different formats. 
#' 
#' Plant phenology is divided in three categories: 
#' early, middle, late. 
#' 
#' Animal phenology is given by four columns per sp: 
#' min.month, min.day, max.month, max.day. 
#' These identify the activity period of each sp each year.
#' 
#' Nesting information is given in categories, one category per taxon;
#' larval feeding information likewise. Taxa in the same category are assumed
#' to compete.
#' 
#' If randomize, inter-guild links are reshuffled keeping row and column totals 
#' fixed.
#' 
#' If mean field, intraguild matrices are mean field matrices, with values given
#' by mean.field.offdiag and mean.field.diag
#'
#' @param pp.all.years nested list, by year and plot, with number of spatial associations among plant sp
#' @param ph.all.years nested list, by year and plot, with number of visits of herbivores (columns) to plants (rows)
#' @param pfv.all.years nested list, by year and plot, with number of visits of floral visitors (columns) to plants (rows)
#' @param plant.phenology dataframe with two columns. ID and pheno.cat, giving the phenology category of each plant: early, middle, late
#' @param animal.phenology dataframe with ate least six columns: ID, year, min.month, min.day, max.month, max.day
#' @param animal.nesting.info dataframe with two columns: ID and nesting (categorical)
#' @param animal.larval.info dataframe with two columns: ID and larval.food.requirements (categorical)
#' @param randomize TRUE or FALSE
#' @param time.limit time limit in seconds for running bipartite::swap.web, which is the underlying randomization function.
#' this is set because that function can get stuck. If it runs above the limit, it starts again.
#' @param intraguild.type factor: mean_field, phenology, nesting, larvae, phenology_nesting_larvae
#' @param mean.field.offdiag numeric
#' @param mean.field.diag numeric
#'
#' @return nested list with three levels. In the first level, element 1 is the list of matrices, 
#' element 2 is the names of the community. In the second level, the list of matrices is nested by year and plot, 
#' with each element being a block matrix for the full community
#' 
#' @export
#'
#' @examples
aux_combine_matrices <- function(pp.all.years,
                                 ph.all.years,
                                 pfv.all.years,
                                 plant.phenology = NULL,
                                 animal.phenology,
                                 animal.nesting.info,
                                 animal.larval.info,
                                 # taxo.in,
                                 # include.overlap = TRUE,
                                 randomize = FALSE,
                                 time.limit = 10,
                                 intraguild.type = c("mean_field","phenology",
                                                     "nesting",
                                                     "larvae",
                                                     "phenology_nesting_larvae"),
                                 mean.field.offdiag = NULL,
                                 mean.field.diag = NULL,
                                 verbose = FALSE){
  
  
  # -------------------------------------------------------------------------
  
  # aux function
  range01 <- function(x){(x-min(x))/(max(x)-min(x))}
  
  # -------------------------------------------------------------------------
  # simple function to randomize a matrix keeping the diagonal elements fixed
  # it simply reshuffles the non-diagonal elements randomly, therefore
  # does not keep marginals or node degrees, but it does keep connectance.
  # Include a flag for non-square matrices, in which diagonal elements are not
  # intra-specific coefs.
  
  reshuffle_matrix <- function(A, keep.diag = TRUE){
    A_rand <- A
    
    if(keep.diag){
      non.diag <- A[row(A) != col(A)]
      A_rand[row(A_rand) != col(A_rand)] <- sample(non.diag,length(non.diag))
    }else{
      A_rand[] <- sample(A,length(A))
    }
    return(A_rand)
  }
  
  # -------------------------------------------------------------------------
  # auxiliary function for randomizing matrices while preventing it from
  # getting stuck
  # swap.fun <- function(A,time.limit, verbose = FALSE){
  #   
  #   res <- tryCatch(
  #     {
  #       R.utils::withTimeout(bipartite::swap.web(1,(A*1e3)), timeout = time.limit)                    
  #     }
  #     ,TimeoutException = function(ex) NA
  #   )
  #   # trycatch returns a list, check if error and return the matrix
  #   if(sum(is.na(res[[1]])) == 0) res <- res[[1]]/1e3
  #   
  #   # if it does not work, start over again
  #   # up to #counter times
  #   counter <- 0
  #   while(sum(is.na(res))>0 & counter < 5){
  #     if(verbose) print(paste("swap.web restarted for the ",counter, " time",sep=""))
  #     res <- tryCatch(
  #       {
  #         R.utils::withTimeout(bipartite::swap.web(1,(A*1e3)), timeout = time.limit)                    
  #       }
  #       ,TimeoutException = function(ex) NA
  #     )
  #     counter <- counter + 1
  #     # trycatch returns a list, check if error and return the matrix
  #     if(sum(is.na(res[[1]])) == 0) res <- res[[1]]/1e3
  #   }
  #   return <- res
  # }
  # -------------------------------------------------------------------------
  
  # constants
  years <- names(pp.all.years)
  plots <- 1:length(pp.all.years[[1]])
  
  # -------------------------------------------------------------------------
  # convert every intraguild matrix to a square matrix
  # this is inherited from the previous version, and maintained simply to
  # obtain the rownames and colnames when building the intraguild matrices
  # note that I don't use the actual values of these "observed" visit matrices anymore
  
  # 1 - floral visitors
  fv.all.years <- pfv.all.years
  
  for(i.year in 1:length(years)){
    for(i.plot in 1:length(plots)){
      
      A <- pfv.all.years[[i.year]][[i.plot]]
      
      visits_matrix <- matrix(0,nrow = ncol(A),ncol = ncol(A),
                              dimnames = list(colnames(A),colnames(A)))
      # diag: colsum
      # non-diag: first, select non-zero rows (rows with shared visits)
      # second, sum of the minimums (min number of shared visits)
      for(i in 1:nrow(visits_matrix)){
        for(j in 1:ncol(visits_matrix)){
          if(i == j){
            visits_matrix[i,j] <- colSums(A)[i]
          }else{
            # only selected columns
            subm <- A[,c(i,j)]
            # # only non-zero rows
            # subm <- subm[apply(subm,1,FUN = function(x){all(x != 0)}),]
            # note: the above is not necessary, actually
            # sum of the minimum values
            visits_matrix[i,j] <- sum(apply(subm,1,FUN = function(x){min(x)}))
          }
        }# for j
      }# for i
      
      fv.all.years[[i.year]][[i.plot]] <- visits_matrix
      
    }# for i.plot
  }# for i.year
  
  # 2 - herbivores
  h.all.years <- ph.all.years
  
  for(i.year in 1:length(years)){
    for(i.plot in 1:length(plots)){
      
      A <- ph.all.years[[i.year]][[i.plot]]
      
      visits_matrix <- matrix(0,nrow = ncol(A),ncol = ncol(A),
                              dimnames = list(colnames(A),colnames(A)))
      # diag: colsum
      # non-diag: first, select non-zero rows (rows with shared visits)
      # second, sum of the minimums (min number of shared visits)
      for(i in 1:nrow(visits_matrix)){
        for(j in 1:ncol(visits_matrix)){
          if(i == j){
            visits_matrix[i,j] <- colSums(A)[i]
          }else{
            # only selected columns
            subm <- A[,c(i,j)]
            # # only non-zero rows
            # subm <- subm[apply(subm,1,FUN = function(x){all(x != 0)}),]
            # note: the above is not necessary, actually
            # sum of the minimum values
            visits_matrix[i,j] <- sum(apply(subm,1,FUN = function(x){min(x)}))
          }
        }# for j
      }# for i
      
      h.all.years[[i.year]][[i.plot]] <- visits_matrix
      
    }# for i.plot
  }# for i.year
  
  # -------------------------------------------------------------------------
  # generate intraguild matrices depending on the mechanism
  
  # store them in a pre-allocated list
  p_intraguild <- list()
  fv_intraguild <- list()
  h_intraguild <- list()
  
  for(i.year in 1:length(years)){
    
    p_intraguild[[i.year]] <- list()
    fv_intraguild[[i.year]] <- list()
    h_intraguild[[i.year]] <- list()
    
    for(i.plot in 1:length(plots)){
      p_intraguild[[i.year]][[i.plot]] <- NA
      fv_intraguild[[i.year]][[i.plot]] <- NA
      h_intraguild[[i.year]][[i.plot]] <- NA
    }
  }
  names(p_intraguild) <- years
  names(fv_intraguild) <- years
  names(h_intraguild) <- years
  
  if(intraguild.type == "mean_field"){
    
    for(i.year in 1:length(years)){
      for(i.plot in 1:length(plots)){
        
        # plants
        pp.num <- nrow(pp.all.years[[i.year]][[i.plot]]) * 
          ncol(pp.all.years[[i.year]][[i.plot]])
        
        p_intraguild[[i.year]][[i.plot]] <- 
          # matrix(data = rep(mean(pp.all.years[[i.year]][[i.plot]]),pp.num),
          matrix(data = rep(mean.field.offdiag,pp.num),
                 nrow = nrow(pp.all.years[[i.year]][[i.plot]]),
                 dimnames = list(rownames(pp.all.years[[i.year]][[i.plot]]),
                                 colnames(pp.all.years[[i.year]][[i.plot]])))
        
        diag(p_intraguild[[i.year]][[i.plot]]) <- mean.field.diag
        
        # floral visitors and herbivores
        plot.fv.names <- colnames(fv.all.years[[i.year]][[i.plot]])
        plot.h.names <- colnames(h.all.years[[i.year]][[i.plot]])
        
        fv_intraguild[[i.year]][[i.plot]] <- 
          matrix(data = rep(mean.field.offdiag,length(plot.fv.names)^2),
                 nrow = length(plot.fv.names),
                 dimnames = list(plot.fv.names,plot.fv.names))
        
        diag(fv_intraguild[[i.year]][[i.plot]]) <- mean.field.diag
        
        h_intraguild[[i.year]][[i.plot]] <-
          matrix(data = rep(mean.field.offdiag,length(plot.h.names)^2),
                 nrow = length(plot.h.names),
                 dimnames = list(plot.h.names,plot.h.names))
        
        diag(h_intraguild[[i.year]][[i.plot]]) <- mean.field.diag
        
      }# for i.plot
    }# for i.year
    
  }else if(intraguild.type == "phenology"){
    
    # the overlap is slightly different for plants and for animal communities
    # plants have three phenological modes, see "plant.phenology"
    # so that sp in same category have overlap 1, adjacent 0.5, and other 0
    
    # the phenological overlap is constant across years for plants
    # i.e. early, middle, and late species are conserved in different years
    
    # competition is given, for plants, by spatial associations weighted
    # by the phenological overlap
    for(i.year in 1:length(years)){
      
      for(i.plot in 1:length(plots)){
        my.plants <- sort(unique(row.names(pp.all.years[[i.year]][[i.plot]])))
        p.overlap.matrix <- matrix(0,nrow = nrow(pp.all.years[[i.year]][[i.plot]]),
                                   ncol = ncol(pp.all.years[[i.year]][[i.plot]]),
                                   dimnames = list(my.plants,my.plants))
        
        # fill up plant overlap matrix
        for(i.plant in 1:nrow(p.overlap.matrix)){
          for(j.plant in 1:ncol(p.overlap.matrix)){
            icat <- plant.phenology$pheno.cat[plant.phenology$ID == rownames(p.overlap.matrix)[i.plant]]
            jcat <- plant.phenology$pheno.cat[plant.phenology$ID == rownames(p.overlap.matrix)[j.plant]]
            
            if(icat == jcat){
              my.overlap <- 1
            }else if(icat == "middle" | jcat == "middle"){
              my.overlap <- 0.5
            }else{
              my.overlap <- 0
            }
            
            p.overlap.matrix[i.plant,j.plant] <- my.overlap
            
          }# for j
        }# for i
        
        # weight the plant observations by phenological overlap
        p_intraguild[[i.year]][[i.plot]] <- pp.all.years[[i.year]][[i.plot]] * p.overlap.matrix
        
      }# for i.plot
    }# for i.year
    
    # for animals, obtain relative phenological overlap
    pheno.matrix <- list()
    
    for(i.year in 1:length(years)){
      my.sp.data <- subset(animal.phenology,!is.na(min.month) & year == years[i.year])
      pheno.matrix[[i.year]] <- get_phenologic_overlap(sp.data = my.sp.data)
    }# for i.year
    
    # here, apply the phenology mask to each pollinator and each herbivore community
    # thus, the competition coefficient is *only* the phenological overlap
    for(i.year in 1:length(years)){
      
      for(i.plot in 1:length(plots)){
        
        plot.fv.names <- colnames(fv.all.years[[i.year]][[i.plot]])
        plot.h.names <- colnames(h.all.years[[i.year]][[i.plot]])
        
        fv.plot.mask <- pheno.matrix[[i.year]][plot.fv.names,plot.fv.names]
        h.plot.mask <- pheno.matrix[[i.year]][plot.h.names,plot.h.names]
        
        fv_intraguild[[i.year]][[i.plot]] <- fv.plot.mask
        h_intraguild[[i.year]][[i.plot]] <- h.plot.mask
        
      }# for i.plot
    }# for i.year
    
  }else if(intraguild.type == "nesting"){
    
    # update fv_overlap and h_overlap:
    # only taxa with same nesting requirements compete
    # for plants, mean field
    
    for(i.year in 1:length(years)){
      for(i.plot in 1:length(plots)){
        
        # plants
        pp.num <- nrow(pp.all.years[[i.year]][[i.plot]]) * 
          ncol(pp.all.years[[i.year]][[i.plot]])
        
        p_intraguild[[i.year]][[i.plot]] <- 
          # matrix(data = rep(mean(pp.all.years[[i.year]][[i.plot]]),pp.num),
          matrix(data = rep(mean.field.offdiag,pp.num),
                 nrow = nrow(pp.all.years[[i.year]][[i.plot]]),
                 dimnames = list(rownames(pp.all.years[[i.year]][[i.plot]]),
                                 colnames(pp.all.years[[i.year]][[i.plot]])))
        
        diag(p_intraguild[[i.year]][[i.plot]]) <- mean.field.diag
        
        # animals
        my.fv <- sort(unique(row.names(fv.all.years[[i.year]][[i.plot]])))
        fv.nest.overlap.matrix <- matrix(0,nrow = nrow(fv.all.years[[i.year]][[i.plot]]),
                                         ncol = ncol(fv.all.years[[i.year]][[i.plot]]),
                                         dimnames = list(my.fv,my.fv))
        
        my.h <- sort(unique(row.names(h.all.years[[i.year]][[i.plot]])))
        h.nest.overlap.matrix <- matrix(0,nrow = nrow(h.all.years[[i.year]][[i.plot]]),
                                        ncol = ncol(h.all.years[[i.year]][[i.plot]]),
                                        dimnames = list(my.h,my.h))
        
        # fill up nest overlap matrices and apply the mask
        # floral visitors
        for(i.fv in 1:nrow(fv.nest.overlap.matrix)){
          for(j.fv in 1:ncol(fv.nest.overlap.matrix)){
            icat <- animal.nesting.info$nesting[animal.nesting.info$ID == rownames(fv.nest.overlap.matrix)[i.fv]]
            jcat <- animal.nesting.info$nesting[animal.nesting.info$ID == rownames(fv.nest.overlap.matrix)[j.fv]]
            
            if(!is.na(jcat) & !is.na(icat)){
              if(icat == jcat){
                my.overlap <- mean.field.offdiag
              }else{
                my.overlap <- 0
              }
            }else{
              my.overlap <- 0
            }
            
            fv.nest.overlap.matrix[i.fv,j.fv] <- my.overlap
            
          }# for j
        }# for i
        
        diag(fv.nest.overlap.matrix) <- mean.field.diag
        
        # herbivores
        for(i.h in 1:nrow(h.nest.overlap.matrix)){
          for(j.h in 1:ncol(h.nest.overlap.matrix)){
            icat <- animal.nesting.info$nesting[animal.nesting.info$ID == rownames(h.nest.overlap.matrix)[i.h]]
            jcat <- animal.nesting.info$nesting[animal.nesting.info$ID == rownames(h.nest.overlap.matrix)[j.h]]
            
            if(!is.na(jcat) & !is.na(icat)){
              if(icat == jcat){
                my.overlap <- mean.field.offdiag
              }else{
                my.overlap <- 0
              }
            }else{
              my.overlap <- 0
            }
            
            h.nest.overlap.matrix[i.h,j.h] <- my.overlap
            
          }# for j
        }# for i
        
        diag(h.nest.overlap.matrix) <- mean.field.diag
        
        # only nesting overlap
        
        fv_intraguild[[i.year]][[i.plot]] <- fv.nest.overlap.matrix
        h_intraguild[[i.year]][[i.plot]] <- h.nest.overlap.matrix
        
      }# for i.plot
    }# for i.year
    
  }else if(intraguild.type == "larvae"){
    
    # equivalent to nesting
    
    for(i.year in 1:length(years)){
      for(i.plot in 1:length(plots)){
        
        # plants
        pp.num <- nrow(pp.all.years[[i.year]][[i.plot]]) * 
          ncol(pp.all.years[[i.year]][[i.plot]])
        
        p_intraguild[[i.year]][[i.plot]] <- 
          # matrix(data = rep(mean(pp.all.years[[i.year]][[i.plot]]),pp.num),
          matrix(data = rep(mean.field.offdiag,pp.num),
                 nrow = nrow(pp.all.years[[i.year]][[i.plot]]),
                 dimnames = list(rownames(pp.all.years[[i.year]][[i.plot]]),
                                 colnames(pp.all.years[[i.year]][[i.plot]])))
        
        diag(p_intraguild[[i.year]][[i.plot]]) <- mean.field.diag
        
        # animals
        my.fv <- sort(unique(row.names(fv.all.years[[i.year]][[i.plot]])))
        fv.larval.overlap.matrix <- matrix(0,nrow = nrow(fv.all.years[[i.year]][[i.plot]]),
                                           ncol = ncol(fv.all.years[[i.year]][[i.plot]]),
                                           dimnames = list(my.fv,my.fv))
        
        my.h <- sort(unique(row.names(h.all.years[[i.year]][[i.plot]])))
        h.larval.overlap.matrix <- matrix(0,nrow = nrow(h.all.years[[i.year]][[i.plot]]),
                                          ncol = ncol(h.all.years[[i.year]][[i.plot]]),
                                          dimnames = list(my.h,my.h))
        
        # fill up nest overlap matrices and apply the mask
        # floral visitors
        for(i.fv in 1:nrow(fv.larval.overlap.matrix)){
          for(j.fv in 1:ncol(fv.larval.overlap.matrix)){
            icat <- animal.larval.info$larval.food.requirements[animal.larval.info$ID == rownames(fv.larval.overlap.matrix)[i.fv]]
            jcat <- animal.larval.info$larval.food.requirements[animal.larval.info$ID == rownames(fv.larval.overlap.matrix)[j.fv]]
            
            if(!is.na(jcat) & !is.na(icat)){
              if(icat == jcat){
                my.overlap <- mean.field.offdiag
              }else{
                my.overlap <- 0
              }
            }else{
              my.overlap <- 0
            }
            
            fv.larval.overlap.matrix[i.fv,j.fv] <- my.overlap
            
          }# for j
        }# for i
        
        diag(fv.larval.overlap.matrix) <- mean.field.diag
        
        # herbivores
        for(i.h in 1:nrow(h.larval.overlap.matrix)){
          for(j.h in 1:ncol(h.larval.overlap.matrix)){
            icat <- animal.larval.info$larval.food.requirements[animal.larval.info$ID == rownames(h.larval.overlap.matrix)[i.h]]
            jcat <- animal.larval.info$larval.food.requirements[animal.larval.info$ID == rownames(h.larval.overlap.matrix)[j.h]]
            
            if(!is.na(jcat) & !is.na(icat)){
              if(icat == jcat){
                my.overlap <- mean.field.offdiag
              }else{
                my.overlap <- 0
              }
            }else{
              my.overlap <- 0
            }
            
            h.larval.overlap.matrix[i.h,j.h] <- my.overlap
            
          }# for j
        }# for i
        
        diag(h.larval.overlap.matrix) <- mean.field.diag
        
        # only larval overlap
        
        fv_intraguild[[i.year]][[i.plot]] <- fv.larval.overlap.matrix
        h_intraguild[[i.year]][[i.plot]] <- h.larval.overlap.matrix
        
      }# for i.plot
    }# for i.year
    
  }else if(intraguild.type == "nesting_larvae"){
    
    for(i.year in 1:length(years)){
      for(i.plot in 1:length(plots)){
        
        # plants 
        # normalise plant observations to make them comparable with animals
        # apply works weird, joins by column, hence the transpose
        p_intraguild[[i.year]][[i.plot]] <- t(apply(pp.all.years[[i.year]][[i.plot]],MARGIN = 1,FUN = range01,simplify = T))
        diag(p_intraguild[[i.year]][[i.plot]]) <- mean.field.diag
        
        # nesting and larval competition
        
        my.fv <- sort(unique(row.names(fv.all.years[[i.year]][[i.plot]])))
        fv.larval.overlap.matrix <- matrix(0,nrow = nrow(fv.all.years[[i.year]][[i.plot]]),
                                           ncol = ncol(fv.all.years[[i.year]][[i.plot]]),
                                           dimnames = list(my.fv,my.fv))
        
        my.h <- sort(unique(row.names(h.all.years[[i.year]][[i.plot]])))
        h.larval.overlap.matrix <- matrix(0,nrow = nrow(h.all.years[[i.year]][[i.plot]]),
                                          ncol = ncol(h.all.years[[i.year]][[i.plot]]),
                                          dimnames = list(my.h,my.h))
        
        # fill up larval overlap matrices and apply the mask
        # floral visitors
        for(i.fv in 1:nrow(fv.larval.overlap.matrix)){
          for(j.fv in 1:ncol(fv.larval.overlap.matrix)){
            icat <- animal.larval.info$larval.food.requirements[animal.larval.info$ID == rownames(fv.larval.overlap.matrix)[i.fv]]
            jcat <- animal.larval.info$larval.food.requirements[animal.larval.info$ID == rownames(fv.larval.overlap.matrix)[j.fv]]
            
            if(!is.na(jcat) & !is.na(icat)){
              if(icat == jcat){
                my.overlap <- mean.field.offdiag
              }else{
                my.overlap <- 0
              }
            }else{
              my.overlap <- 0
            }
            
            fv.larval.overlap.matrix[i.fv,j.fv] <- my.overlap
            
          }# for j
        }# for i
        
        diag(fv.larval.overlap.matrix) <- mean.field.diag
        
        # herbivores
        for(i.h in 1:nrow(h.larval.overlap.matrix)){
          for(j.h in 1:ncol(h.larval.overlap.matrix)){
            icat <- animal.larval.info$larval.food.requirements[animal.larval.info$ID == rownames(h.larval.overlap.matrix)[i.h]]
            jcat <- animal.larval.info$larval.food.requirements[animal.larval.info$ID == rownames(h.larval.overlap.matrix)[j.h]]
            
            if(!is.na(jcat) & !is.na(icat)){
              if(icat == jcat){
                my.overlap <- mean.field.offdiag
              }else{
                my.overlap <- 0
              }
            }else{
              my.overlap <- 0
            }
            
            h.larval.overlap.matrix[i.h,j.h] <- my.overlap
            
          }# for j
        }# for i
        
        diag(h.larval.overlap.matrix) <- mean.field.diag
        
        # nest
        fv.nest.overlap.matrix <- matrix(0,nrow = nrow(fv.all.years[[i.year]][[i.plot]]),
                                         ncol = ncol(fv.all.years[[i.year]][[i.plot]]),
                                         dimnames = list(my.fv,my.fv))
        
        h.nest.overlap.matrix <- matrix(0,nrow = nrow(h.all.years[[i.year]][[i.plot]]),
                                        ncol = ncol(h.all.years[[i.year]][[i.plot]]),
                                        dimnames = list(my.h,my.h))
        
        # fill up nest overlap matrices and apply the mask
        # floral visitors
        for(i.fv in 1:nrow(fv.nest.overlap.matrix)){
          for(j.fv in 1:ncol(fv.nest.overlap.matrix)){
            icat <- animal.nesting.info$nesting[animal.nesting.info$ID == rownames(fv.nest.overlap.matrix)[i.fv]]
            jcat <- animal.nesting.info$nesting[animal.nesting.info$ID == rownames(fv.nest.overlap.matrix)[j.fv]]
            
            if(!is.na(jcat) & !is.na(icat)){
              if(icat == jcat){
                my.overlap <- mean.field.offdiag
              }else{
                my.overlap <- 0
              }
            }else{
              my.overlap <- 0
            }
            
            fv.nest.overlap.matrix[i.fv,j.fv] <- my.overlap
            
          }# for j
        }# for i
        
        diag(fv.nest.overlap.matrix) <- mean.field.diag
        
        # herbivores
        for(i.h in 1:nrow(h.nest.overlap.matrix)){
          for(j.h in 1:ncol(h.nest.overlap.matrix)){
            icat <- animal.nesting.info$nesting[animal.nesting.info$ID == rownames(h.nest.overlap.matrix)[i.h]]
            jcat <- animal.nesting.info$nesting[animal.nesting.info$ID == rownames(h.nest.overlap.matrix)[j.h]]
            
            if(!is.na(jcat) & !is.na(icat)){
              if(icat == jcat){
                my.overlap <- mean.field.offdiag
              }else{
                my.overlap <- 0
              }
            }else{
              my.overlap <- 0
            }
            
            h.nest.overlap.matrix[i.h,j.h] <- my.overlap
            
          }# for j
        }# for i
        
        diag(h.nest.overlap.matrix) <- mean.field.diag
        
        # final step: phenology + nesting + larval
        
        fv_intraguild[[i.year]][[i.plot]] <- (fv.larval.overlap.matrix + fv.nest.overlap.matrix)/2
        h_intraguild[[i.year]][[i.plot]] <- (h.larval.overlap.matrix + h.nest.overlap.matrix)/2
        
      }# for i.plot
    }# for i.year
    
  }else if(intraguild.type == "nesting_larvae_phenology"){
    
    # first obtain phenological overlap, 
    # then obtain nesting and larval proxies,
    # lastly combine them: ((nest+larval)/2) * phenology
    # in this way, phenology modulates the coefficients from nest+larval
    # for plants, phenology * number of associations
    for(i.year in 1:length(years)){
      for(i.plot in 1:length(plots)){
        my.plants <- sort(unique(row.names(pp.all.years[[i.year]][[i.plot]])))
        p.overlap.matrix <- matrix(0,nrow = nrow(pp.all.years[[i.year]][[i.plot]]),
                                   ncol = ncol(pp.all.years[[i.year]][[i.plot]]),
                                   dimnames = list(my.plants,my.plants))
        
        # fill up plant overlap matrix
        for(i.plant in 1:nrow(p.overlap.matrix)){
          for(j.plant in 1:ncol(p.overlap.matrix)){
            icat <- plant.phenology$pheno.cat[plant.phenology$ID == rownames(p.overlap.matrix)[i.plant]]
            jcat <- plant.phenology$pheno.cat[plant.phenology$ID == rownames(p.overlap.matrix)[j.plant]]
            
            if(icat == jcat){
              my.overlap <- 1
            }else if(icat == "middle" | jcat == "middle"){
              my.overlap <- 0.5
            }else{
              my.overlap <- 0
            }
            
            p.overlap.matrix[i.plant,j.plant] <- my.overlap
            
          }# for j
        }# for i
        
        # normalise plant observations to make them comparable with animals
        # apply works weird, joins by column, hence the transpose
        norm.plants.matrix <- t(apply(pp.all.years[[i.year]][[i.plot]],MARGIN = 1,FUN = range01,simplify = T))
        diag(norm.plants.matrix) <- mean.field.diag
        
        # weight the plant observations by phenological overlap
        p_intraguild[[i.year]][[i.plot]] <- norm.plants.matrix * p.overlap.matrix
        
      }# for i.plot
    }# for i.year
    
    # animals
    pheno.matrix <- list()
    
    for(i.year in 1:length(years)){
      my.sp.data <- subset(animal.phenology,!is.na(min.month) & year == years[i.year])
      pheno.matrix[[i.year]] <- get_phenologic_overlap(sp.data = my.sp.data)
    }# for i.year
    
    for(i.year in 1:length(years)){
      
      for(i.plot in 1:length(plots)){
        
        # first, obtain the phenology overlap for each pollinator and each herbivore community
        # values 0-1, asymmetrical
        
        plot.fv.names <- colnames(fv.all.years[[i.year]][[i.plot]])
        plot.h.names <- colnames(h.all.years[[i.year]][[i.plot]])
        
        fv.plot.mask <- pheno.matrix[[i.year]][plot.fv.names,plot.fv.names]
        h.plot.mask <- pheno.matrix[[i.year]][plot.h.names,plot.h.names]
        
        # obtain nesting and larval intraguild proxies
        # for floral visitors and herbivores
        
        my.fv <- sort(unique(row.names(fv.all.years[[i.year]][[i.plot]])))
        fv.larval.overlap.matrix <- matrix(0,nrow = nrow(fv.all.years[[i.year]][[i.plot]]),
                                           ncol = ncol(fv.all.years[[i.year]][[i.plot]]),
                                           dimnames = list(my.fv,my.fv))
        
        my.h <- sort(unique(row.names(h.all.years[[i.year]][[i.plot]])))
        h.larval.overlap.matrix <- matrix(0,nrow = nrow(h.all.years[[i.year]][[i.plot]]),
                                          ncol = ncol(h.all.years[[i.year]][[i.plot]]),
                                          dimnames = list(my.h,my.h))
        
        # fill up larval overlap matrices and apply the mask
        # floral visitors
        for(i.fv in 1:nrow(fv.larval.overlap.matrix)){
          for(j.fv in 1:ncol(fv.larval.overlap.matrix)){
            icat <- animal.larval.info$larval.food.requirements[animal.larval.info$ID == rownames(fv.larval.overlap.matrix)[i.fv]]
            jcat <- animal.larval.info$larval.food.requirements[animal.larval.info$ID == rownames(fv.larval.overlap.matrix)[j.fv]]
            
            if(!is.na(jcat) & !is.na(icat)){
              if(icat == jcat){
                my.overlap <- mean.field.offdiag
              }else{
                my.overlap <- 0
              }
            }else{
              my.overlap <- 0
            }
            
            fv.larval.overlap.matrix[i.fv,j.fv] <- my.overlap
            
          }# for j
        }# for i
        
        diag(fv.larval.overlap.matrix) <- mean.field.diag
        
        # herbivores
        for(i.h in 1:nrow(h.larval.overlap.matrix)){
          for(j.h in 1:ncol(h.larval.overlap.matrix)){
            icat <- animal.larval.info$larval.food.requirements[animal.larval.info$ID == rownames(h.larval.overlap.matrix)[i.h]]
            jcat <- animal.larval.info$larval.food.requirements[animal.larval.info$ID == rownames(h.larval.overlap.matrix)[j.h]]
            
            if(!is.na(jcat) & !is.na(icat)){
              if(icat == jcat){
                my.overlap <- mean.field.offdiag
              }else{
                my.overlap <- 0
              }
            }else{
              my.overlap <- 0
            }
            
            h.larval.overlap.matrix[i.h,j.h] <- my.overlap
            
          }# for j
        }# for i
        
        diag(h.larval.overlap.matrix) <- mean.field.diag
        
        # nest
        fv.nest.overlap.matrix <- matrix(0,nrow = nrow(fv.all.years[[i.year]][[i.plot]]),
                                         ncol = ncol(fv.all.years[[i.year]][[i.plot]]),
                                         dimnames = list(my.fv,my.fv))
        
        h.nest.overlap.matrix <- matrix(0,nrow = nrow(h.all.years[[i.year]][[i.plot]]),
                                        ncol = ncol(h.all.years[[i.year]][[i.plot]]),
                                        dimnames = list(my.h,my.h))
        
        # nesting
        # floral visitors
        for(i.fv in 1:nrow(fv.nest.overlap.matrix)){
          for(j.fv in 1:ncol(fv.nest.overlap.matrix)){
            icat <- animal.nesting.info$nesting[animal.nesting.info$ID == rownames(fv.nest.overlap.matrix)[i.fv]]
            jcat <- animal.nesting.info$nesting[animal.nesting.info$ID == rownames(fv.nest.overlap.matrix)[j.fv]]
            
            if(!is.na(jcat) & !is.na(icat)){
              if(icat == jcat){
                my.overlap <- mean.field.offdiag
              }else{
                my.overlap <- 0
              }
            }else{
              my.overlap <- 0
            }
            
            fv.nest.overlap.matrix[i.fv,j.fv] <- my.overlap
            
          }# for j
        }# for i
        
        diag(fv.nest.overlap.matrix) <- mean.field.diag
        
        # herbivores
        for(i.h in 1:nrow(h.nest.overlap.matrix)){
          for(j.h in 1:ncol(h.nest.overlap.matrix)){
            icat <- animal.nesting.info$nesting[animal.nesting.info$ID == rownames(h.nest.overlap.matrix)[i.h]]
            jcat <- animal.nesting.info$nesting[animal.nesting.info$ID == rownames(h.nest.overlap.matrix)[j.h]]
            
            if(!is.na(jcat) & !is.na(icat)){
              if(icat == jcat){
                my.overlap <- mean.field.offdiag
              }else{
                my.overlap <- 0
              }
            }else{
              my.overlap <- 0
            }
            
            h.nest.overlap.matrix[i.h,j.h] <- my.overlap
            
          }# for j
        }# for i
        
        diag(h.nest.overlap.matrix) <- mean.field.diag
        
        # final step: phenology modulates the (nesting + larval) matrix
        
        fv_intraguild[[i.year]][[i.plot]] <- ((fv.larval.overlap.matrix + fv.nest.overlap.matrix)/2)* fv.plot.mask
        h_intraguild[[i.year]][[i.plot]] <- ((h.larval.overlap.matrix + h.nest.overlap.matrix)/2)* h.plot.mask
        
      }# for i.plot
    }# for i.year
    
  }else{
    message("function aux_combine_matrices ERROR: please set a valid intraguild.type")
    return(NULL)
  }
  
  
  # -------------------------------------------------------------------------
  # round links to three decimal places
  # this makes no difference in the feasibility domains obtained 
  # and is necessary for circumventing the swap.web algorithm limitations
  # of working with integers when randomizing.
  for(i.year in 1:length(years)){
    for(i.plot in plots){
      p_intraguild[[i.year]][[i.plot]] <- round(p_intraguild[[i.year]][[i.plot]],3)
      fv_intraguild[[i.year]][[i.plot]] <- round(fv_intraguild[[i.year]][[i.plot]],3)
      h_intraguild[[i.year]][[i.plot]] <- round(h_intraguild[[i.year]][[i.plot]],3)
      
      pfv.all.years[[i.year]][[i.plot]] <- round(pfv.all.years[[i.year]][[i.plot]],3)
      ph.all.years[[i.year]][[i.plot]] <- round(ph.all.years[[i.year]][[i.plot]],3)
    }
  }
  
  # -------------------------------------------------------------------------
  # randomize here, when all square matrices are built.
  # this can be done here because randomizations of observed interactions 
  # i.e. plant-herb and plant-fv, do not affect the intraguild matrices
  # already calculated. For plant-plant interactions, the matrix used 
  # is p_overlap, and it does not influence other matrices, so I can
  # randomize it as well.
  
  # In this version I randomize using bipartite::swap.web, to constrain
  # both row and col marginals and connectance.
  # It works only with integers, so I multiply the original matrices by 1e3, 
  # since they have three decimal places, and then divide again.
  
  if(randomize){
    # pp.all.years.null <- list()
    p_intraguild.null <- vector(mode = "list", length = length(years))
    fv_intraguild.null <- vector(mode = "list", length = length(years))
    h_intraguild.null <- vector(mode = "list", length = length(years))
    ph.all.years.null <- vector(mode = "list", length = length(years))
    pfv.all.years.null <- vector(mode = "list", length = length(years))
    
    # Flag to check whether randomization is successful
    na.matrices <- list()
    
    for(i.year in 1:length(years)){
      na.matrices[[i.year]] <- list() 
      for(i.plot in 1:length(plots)){
        na.matrices[[i.year]][[i.plot]] <- rep(FALSE,5)
      }
    }
    
    
    for(i.year in 1:length(years)){
      
      # pp.all.years.null[[i.year]] <- list()
      p_intraguild.null[[i.year]] <- vector(mode = "list", length = length(plots))
      fv_intraguild.null[[i.year]] <- vector(mode = "list", length = length(plots))
      h_intraguild.null[[i.year]] <- vector(mode = "list", length = length(plots))
      ph.all.years.null[[i.year]] <- vector(mode = "list", length = length(plots))
      pfv.all.years.null[[i.year]] <- vector(mode = "list", length = length(plots))
      
      for(i.plot in plots){
        
        if(verbose){
          print(paste(i.year,"-",i.plot,"-started",sep=""))
        }
        
        # intraguild matrices
        # plants
        p_intraguild.null[[i.year]][[i.plot]] <- reshuffle_matrix(p_intraguild[[i.year]][[i.plot]])
          # (bipartite::swap.web(1,(p_intraguild[[i.year]][[i.plot]])*1e3)[[1]])/1e3
        #   swap.fun(p_intraguild[[i.year]][[i.plot]],time.limit)
        # 
        # if(all(!is.na(p_intraguild.null[[i.year]][[i.plot]]))){
        #   colnames(p_intraguild.null[[i.year]][[i.plot]]) <- 
        #     colnames(p_intraguild[[i.year]][[i.plot]])
        #   rownames(p_intraguild.null[[i.year]][[i.plot]]) <- 
        #     rownames(p_intraguild[[i.year]][[i.plot]])
        # }else{
        #   na.matrices[[i.year]][[i.plot]][1] <- TRUE
        # }
        # floral visitors
        fv_intraguild.null[[i.year]][[i.plot]] <- reshuffle_matrix(fv_intraguild[[i.year]][[i.plot]])
          # (bipartite::swap.web(1,(fv_intraguild[[i.year]][[i.plot]])*1e3)[[1]])/1e3
        #   swap.fun(fv_intraguild[[i.year]][[i.plot]],time.limit)
        # 
        # if(all(!is.na(fv_intraguild.null[[i.year]][[i.plot]]))){
        #   colnames(fv_intraguild.null[[i.year]][[i.plot]]) <- 
        #     colnames(fv_intraguild[[i.year]][[i.plot]])
        #   rownames(fv_intraguild.null[[i.year]][[i.plot]]) <- 
        #     rownames(fv_intraguild[[i.year]][[i.plot]])
        # }else{
        #   na.matrices[[i.year]][[i.plot]][2] <- TRUE
        #   if(verbose) print(paste("INTRAGUILD FV year ",i.year,", plot ",i.plot," FAILED",sep=""))
        # }
        
        # herbivores
        h_intraguild.null[[i.year]][[i.plot]] <- reshuffle_matrix(h_intraguild[[i.year]][[i.plot]])
          # (bipartite::swap.web(1,(h_intraguild[[i.year]][[i.plot]])*1e3)[[1]])/1e3
        #   swap.fun(h_intraguild[[i.year]][[i.plot]],time.limit)
        # 
        # if(all(!is.na(h_intraguild.null[[i.year]][[i.plot]]))){
        #   colnames(h_intraguild.null[[i.year]][[i.plot]]) <- 
        #     colnames(h_intraguild[[i.year]][[i.plot]])
        #   rownames(h_intraguild.null[[i.year]][[i.plot]]) <- 
        #     rownames(h_intraguild[[i.year]][[i.plot]])
        # }else{
        #   na.matrices[[i.year]][[i.plot]][3] <- TRUE
        # }
        
        # interguild matrices
        # plant-herbivores
        ph.all.years.null[[i.year]][[i.plot]] <- reshuffle_matrix(ph.all.years[[i.year]][[i.plot]],
                                                                  keep.diag = F)
          # (bipartite::swap.web(1,(ph.all.years[[i.year]][[i.plot]])*1e3)[[1]])/1e3
        #   swap.fun(ph.all.years[[i.year]][[i.plot]],time.limit)
        # 
        # if(all(!is.na(ph.all.years.null[[i.year]][[i.plot]]))){
        #   colnames(ph.all.years.null[[i.year]][[i.plot]]) <- 
        #     colnames(ph.all.years[[i.year]][[i.plot]])
        #   rownames(ph.all.years.null[[i.year]][[i.plot]]) <- 
        #     rownames(ph.all.years[[i.year]][[i.plot]])
        # }else{
        #   na.matrices[[i.year]][[i.plot]][4] <- TRUE
        # }
        
        # plant-floral visitors
        pfv.all.years.null[[i.year]][[i.plot]] <- reshuffle_matrix(pfv.all.years[[i.year]][[i.plot]],
                                                                   keep.diag = F)
          # (bipartite::swap.web(1,(pfv.all.years[[i.year]][[i.plot]])*1e3)[[1]])/1e3
        #   swap.fun(pfv.all.years[[i.year]][[i.plot]],time.limit)
        # 
        # if(all(!is.na(pfv.all.years.null[[i.year]][[i.plot]]))){
        #   colnames(pfv.all.years.null[[i.year]][[i.plot]]) <- 
        #     colnames(pfv.all.years[[i.year]][[i.plot]])
        #   rownames(pfv.all.years.null[[i.year]][[i.plot]]) <- 
        #     rownames(pfv.all.years[[i.year]][[i.plot]])
        # }else{
        #   na.matrices[[i.year]][[i.plot]][5] <- TRUE
        # }
        
      }# for i.plot
    }# for i.year
    
    # names(pp.all.years.null) <- years
    names(p_intraguild.null) <- years
    names(fv_intraguild.null) <- years
    names(h_intraguild.null) <- years
    names(ph.all.years.null) <- years
    names(pfv.all.years.null) <- years
    
    # bad practice, but i don't want to change all the rest
    # pp.all.years <- pp.all.years.null
    p_intraguild <- p_intraguild.null
    fv_intraguild <- fv_intraguild.null
    h_intraguild <- h_intraguild.null
    ph.all.years <- ph.all.years.null
    pfv.all.years <- pfv.all.years.null
    
    # if the swap.web function failed, return NULL
    # if(any(na.matrices)){
    #   message("function aux_combine_matrices: randomization process failed, returning NULL")
    #   # return(NULL)
    # }
    
  }
  
  # -------------------------------------------------------------------------
  # now, standardize the intraguild AND interguild matrices
  # there are different options
  
  # here I implement a standardization independent for each community/interaction type
  # so that each matrix is assumed to have the same amount of "interaction effort" 
  # basically meaning each matrix sums to 1
  
  # natural or log values? there is strong variability among them
  # in principle, natural.
  
  pp.all.years.norm <- p_intraguild
  fv.all.years.norm <- fv_intraguild
  h.all.years.norm <- h_intraguild
  pfv.all.years.norm <- pfv.all.years
  ph.all.years.norm <- ph.all.years
  
  for(i.year in 1:length(years)){
    for(i.plot in 1:length(plots)){
      
      # pp.all.years
      pp.all.years.norm[[i.year]][[i.plot]] <- p_intraguild[[i.year]][[i.plot]]/sum(p_intraguild[[i.year]][[i.plot]])
      
      # fv.all.years
      fv.all.years.norm[[i.year]][[i.plot]] <- fv_intraguild[[i.year]][[i.plot]]/sum(fv_intraguild[[i.year]][[i.plot]])
      
      # h.all.years
      h.all.years.norm[[i.year]][[i.plot]] <- h_intraguild[[i.year]][[i.plot]]/sum(h_intraguild[[i.year]][[i.plot]])
      
      # pfv.all.years
      pfv.all.years.norm[[i.year]][[i.plot]] <- pfv.all.years[[i.year]][[i.plot]]/sum(pfv.all.years[[i.year]][[i.plot]])
      
      # ph.all.years
      ph.all.years.norm[[i.year]][[i.plot]] <- ph.all.years[[i.year]][[i.plot]]/sum(ph.all.years[[i.year]][[i.plot]])
      
    }# for i.plot
  }# for i.year
  
  # impose signs ------------------------------------------------------------
  
  # only positive should be plant-floral visitor
  # and herb-plant!!
  pp_sign <- -1
  pfv_sign <- 1
  ph_sign <- -1
  hp_sign <- 1
  hh_sign <- -1
  fvfv_sign <- -1
  
  # combine all matrices ----------------------------------------------------
  
  community_matrices <- list()
  sp.names <- list()
  
  for(i.year in 1:length(years)){
    
    # year-plot list, each element a block matrix
    community_matrices[[i.year]] <- list()
    
    # also store the set of plants, herbivores, floral visitors
    # per year
    plants.year <- NULL
    herbivores.year <- NULL
    floral.visitors.year <- NULL
    
    for(i.plot in 1:length(plots)){
      
      # components of the block matrix
      pp.matrix <- pp_sign * pp.all.years.norm[[i.year]][[i.plot]]
      ph.matrix <- ph_sign * ph.all.years.norm[[i.year]][[i.plot]]
      pfv.matrix <- pfv_sign * pfv.all.years.norm[[i.year]][[i.plot]]
      fv.overlap.matrix <- fvfv_sign * fv.all.years.norm[[i.year]][[i.plot]]
      h.overlap.matrix <- hh_sign * h.all.years.norm[[i.year]][[i.plot]]
      
      if(all(is.na(pp.matrix)) | all(is.na(ph.matrix)) | 
         all(is.na(pfv.matrix)) | 
         all(is.na(fv.overlap.matrix)) | 
         all(is.na(h.overlap.matrix))){
        
        community_matrices[[i.year]][[i.plot]] <- NA
      }else{
        
        # get the set of species for this community
        valid.names <- get_valid_sp(pp.matrix = pp.matrix,
                                    ph.matrix = ph.matrix,
                                    pfv.matrix = pfv.matrix)
        
        # construct the empty matrix
        block.matrix <- build_block_matrix(plants = valid.names[[1]],
                                           herbivores = valid.names[[2]],
                                           floral.visitors = valid.names[[3]])
        
        # fill the matrix
        # Note that in the function "fill_block_matrix" there is an argument
        # switch_herb_sign that is set to TRUE by default (so not specified here), 
        # that sets plant effects on herbivores as positive, but not the 
        # other way around
        block.matrix.filled <- fill_block_matrix(block.matrix = block.matrix,
                                                 pp.matrix = pp.matrix,
                                                 ph.matrix = ph.matrix,
                                                 pfv.matrix = pfv.matrix,
                                                 fv.overlap.matrix = fv.overlap.matrix,
                                                 h.overlap.matrix = h.overlap.matrix)
        
        # store the matrix
        community_matrices[[i.year]][[i.plot]] <- block.matrix.filled
        
        # store the valid names of this community
        plants.year <- unique(c(plants.year,valid.names[[1]]))
        herbivores.year <- unique(c(herbivores.year,valid.names[[2]]))
        floral.visitors.year <- unique(c(floral.visitors.year,valid.names[[3]]))
      }# if-else is.na
      
    }# for i.plot
    
    if(all(is.na(community_matrices[[i.year]][[i.plot]]))){
      sp.names[[i.year]] <- NA
    }else{
      sp.names[[i.year]] <- list(plants = sort(plants.year),
                                 herbivores = sort(herbivores.year),
                                 floral.visitors = sort(floral.visitors.year))
    }
  }# for i.year
  
  names(community_matrices) <- years
  names(sp.names) <- years
  
  return(list(community_matrices = community_matrices,
              sp.names = sp.names))
  
}
