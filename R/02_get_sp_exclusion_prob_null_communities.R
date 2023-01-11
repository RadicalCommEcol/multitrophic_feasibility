
# obtain species probabilities of exclusion
# from community matrices. This code is computationally intensive and 
# it is parallelized. Note that this script obtains the metrics for the 
# randomized communities.

# NOTE: to ensure reproducibility, you will have to modify the path inside 
# the foreach loop (lines 121-123)

# INPUT
# community matrices: "results/community_matrices.RData"
# community names: "results/community_names.RData"

# OUTPUT
# species exclusion prob dataframe: "results/exclusion_probabilities_null.csv"

# note: it may be more efficient to save each randomization separately, because otherwise
# if the run crashes just for one community, it will wipe out the entire set of results
# hence the write.csv2 bit inside the foreach loop.

# read community matrices -------------------------------------------------
library(tidyverse)
# devtools::install_github("RadicalCommEcol/MultitrophicFun")
library(MultitrophicFun)
library(matlib) # to multiply matrices
library(nleqslv) # to solve non-linear equations
library(zipfR) # beta incomplete function to estimate the area of d-dimensional spherical caps 
library(pracma) # to solve n-dimensional cross products
library(boot) # to bootstrap
library(CValternatives) # to estimate PV index devtools::install_github("T-Engel/CValternatives")

library(foreach)
library(doParallel)

# source auxiliary functions
list.files("R/feasibility_functions/", full.names = TRUE) %>% map(source)

# set number of cores -----------------------------------------------------

workers <- 4
cl <- makeCluster(workers)
# register the cluster for using foreach
registerDoParallel(cl)

# -------------------------------------------------------------------------
# which version do we read?
# vers <- "_v2"
vers <- ""

# -------------------------------------------------------------------------

file.name <- "community_matrices"
file.name <- paste(file.name,vers,".RData",sep="")

# community_matrices[[intraguild.type]][[year]][[plot]]
load(paste("results/",file.name,sep=""))

# positions for the blocks are dynamic, species pool is different for 
# each plot/year
load(paste("results/community_names",vers,".RData",sep=""))

intraguild.types <- names(community_matrices)

# set important constants -------------------------------------------------
# number of null replicates
null.replicates <- 10

# IMPORTANT - set this if running more than once, e.g. to add more replicates
run.number <- 2

# replicates for the feasibility calculations
omega.replicates <- 100
bootstrap.replicates <- 100

# number of noise replicates for the feasibility calculation
# this is for the case in which the original matrix is not invertible
# so that I add a small amount of noise to make it so.
noise.replicates <- 10
noise.type <- "random"
# how much noise to add relative to the minimum observed values
# e.g. an order of magnitude of 1% 
# relative to the minimum of the matrix elements
noise.relative.magnitude <- .1

years <- as.numeric(names(community_matrices[[1]]))
plots <- 1:9
guild.combinations <- c("all")

# this is the combined id to loop over in parallel
id <- expand.grid(intraguild.types,years,plots,guild.combinations)
id.char <- paste(id[,1],"_",id[,2],"_",id[,3],"_",id[,4],sep="")

# in case some are already calculated
null.files <- list.files("results/null_fd/",pattern = "sp_exc*")

if(length(null.files)>0){
  null.files.2 <- substr(null.files,24,nchar(null.files))
  null.files.3 <- substr(null.files.2,1,nchar(null.files.2)-4)
  
  null.run.files <- grepl(paste("run",run.number,sep=""),null.files.3)
  null.files.4 <- null.files.3[!null.run.files]
  null.files.run.2 <- null.files.3[null.run.files]
  null.files.run.3 <- substr(null.files.run.2,1,nchar(null.files.run.2)-11)
  
  null.files.all <- unique(c(null.files.4,null.files.run.3))
  
  id.char <- id.char[which(!id.char %in% null.files.all)]
}
# calculate feasibility metrics  -------------------------------------------

comb.fun <- function(...) {
  mapply('rbind', ..., SIMPLIFY=FALSE)
}

feasibility.metrics.null <- foreach(i.id = 1:length(id.char),
                                    .combine=comb.fun, 
                                    .packages = c("tidyverse","foreach","matlib",
                                                  "nleqslv","zipfR","pracma",
                                                  "boot","CValternatives",
                                                  "MultitrophicFun")) %dopar% {
                                                    
                                                    # cat(id.char[i.id],"- started\n")
                                                    
                                                    list.files("/home/david/Work/Projects/EBD/multitrophic_feasibility/R/feasibility_functions/", 
                                                               full.names = TRUE) %>% map(source)
                                                    source("/home/david/Work/Projects/EBD/multitrophic_feasibility/R/reshuffle_matrix.R")
                                                    
                                                    f <- function(m) class(try(solve(t(m) %*% m),
                                                                               silent = T))[[1]] == "matrix"
                                                    
                                                    # first, recover the parameters of each matrix
                                                    
                                                    my.type <- NA
                                                    my.year <- NA
                                                    my.plot <- NA
                                                    my.guild <- NA
                                                    
                                                    # intraguild matrix type
                                                    if(grepl("mean_field",id.char[i.id])){
                                                      my.type <- "mean_field"
                                                    }else if(grepl("nesting_larvae_phenology",id.char[i.id])){
                                                      my.type <- "nesting_larvae_phenology"
                                                    }else if(grepl("nesting_larvae",id.char[i.id])){
                                                      my.type <- "nesting_larvae"
                                                    }else if(grepl("phenology",id.char[i.id])){
                                                      my.type <- "phenology"
                                                    }else if(grepl("nesting",id.char[i.id])){
                                                      my.type <- "nesting"
                                                    }else if(grepl("larvae",id.char[i.id])){
                                                      my.type <- "larvae"
                                                    }
                                                    
                                                    # guild
                                                    if(grepl("plants-floral visitors",id.char[i.id])){
                                                      my.guild <- "plants-floral visitors"
                                                    }else if(grepl("plants-herbivores",id.char[i.id])){
                                                      my.guild <- "plants-herbivores"
                                                    }else if(grepl("floral visitors",id.char[i.id])){
                                                      my.guild <- "floral visitors"
                                                    }else if(grepl("herbivores",id.char[i.id])){
                                                      my.guild <- "herbivores"
                                                    }else if(grepl("plants",id.char[i.id])){
                                                      my.guild <- "plants"
                                                    }else{
                                                      my.guild <- "all"
                                                    }
                                                    
                                                    # year
                                                    if(grepl("2019",id.char[i.id])){
                                                      my.year <- "2019"
                                                    }else if(grepl("2020",id.char[i.id])){
                                                      my.year <- "2020"
                                                    }
                                                    
                                                    # plot
                                                    my.plot1 <- stringr::str_remove(id.char[i.id],my.type)
                                                    my.plot2 <- stringr::str_remove(my.plot1,my.guild)
                                                    my.plot3 <- stringr::str_remove(my.plot2,my.year)
                                                    my.plot <- as.numeric(stringr::str_remove_all(my.plot3,"_"))
                                                    
                                                    # recover the matrix
                                                    year.plot.matrix <- 
                                                      community_matrices[[my.type]][[my.year]][[my.plot]]
                                                    
                                                    plant.positions <- which(rownames(year.plot.matrix) %in% 
                                                                               sp.names[[my.year]][["plants"]])
                                                    herb.positions <- which(rownames(year.plot.matrix) %in% 
                                                                              sp.names[[my.year]][["herbivores"]])
                                                    fv.positions <- which(rownames(year.plot.matrix) %in% 
                                                                            sp.names[[my.year]][["floral.visitors"]])
                                                    
                                                    if(my.guild == "plants"){
                                                      
                                                      my.matrix <- year.plot.matrix[plant.positions,plant.positions]
                                                      
                                                    }else if(my.guild == "plants-floral visitors"){
                                                      
                                                      my.matrix <- year.plot.matrix[c(plant.positions,fv.positions),
                                                                                    c(plant.positions,fv.positions)]
                                                      
                                                    }else if(my.guild == "floral visitors"){
                                                      
                                                      my.matrix <- year.plot.matrix[fv.positions,fv.positions]
                                                      
                                                    }else if(my.guild == "herbivores"){
                                                      
                                                      my.matrix <- year.plot.matrix[herb.positions,herb.positions]
                                                      
                                                    }else if(my.guild == "plants-herbivores"){
                                                      
                                                      my.matrix <- year.plot.matrix[c(plant.positions,herb.positions),
                                                                                    c(plant.positions,herb.positions)]
                                                      
                                                    }else if(my.guild == "all"){
                                                      
                                                      my.matrix <- year.plot.matrix
                                                    }
                                                    
                                                    # -------------------------------------------------------------------------
                                                    # positions in the subsetted matrix
                                                    my.plant.positions <- which(rownames(my.matrix) %in% 
                                                                                  sp.names[[my.year]][["plants"]])
                                                    my.herb.positions <- which(rownames(my.matrix) %in% 
                                                                                 sp.names[[my.year]][["herbivores"]])
                                                    my.fv.positions <- which(rownames(my.matrix) %in% 
                                                                               sp.names[[my.year]][["floral.visitors"]])
                                                    
                                                    # -------------------------------------------------------------------------
                                                    # randomize
                                                    null.reps.list <- list()
                                                    for(i.null in 1:null.replicates){
                                                      
                                                      my.matrix.null <- my.matrix
                                                      
                                                      # this is to make sure that no interactions are assigned
                                                      # to pollinators-herbivores blocks. In fact,
                                                      # interactions are reshuffled only within their block.
                                                      
                                                      if(my.guild %in% c("plants",
                                                                         "floral visitors",
                                                                         "herbivores")){
                                                        
                                                        my.matrix.null <- reshuffle_matrix(my.matrix)
                                                        
                                                      }else if(my.guild == "plants-floral visitors"){
                                                        pp.block <- reshuffle_matrix(my.matrix.null[my.plant.positions,my.plant.positions])
                                                        fv.block <- reshuffle_matrix(my.matrix.null[my.fv.positions,my.fv.positions])
                                                        pfv.block <- reshuffle_matrix(my.matrix.null[my.plant.positions,my.fv.positions])
                                                        
                                                        my.matrix.null[my.plant.positions,my.plant.positions] <- pp.block
                                                        my.matrix.null[my.fv.positions,my.fv.positions] <- fv.block
                                                        my.matrix.null[my.plant.positions,my.fv.positions] <- pfv.block
                                                        
                                                        # the block matrix is symmetric regarding plant-insect coefs
                                                        my.matrix.null[my.fv.positions,my.plant.positions] <- t(pfv.block)
                                                        
                                                      }else if(my.guild == "plants-herbivores"){
                                                        pp.block <- reshuffle_matrix(my.matrix.null[my.plant.positions,my.plant.positions])
                                                        h.block <- reshuffle_matrix(my.matrix.null[my.herb.positions,my.herb.positions])
                                                        ph.block <- reshuffle_matrix(my.matrix.null[my.plant.positions,my.herb.positions])
                                                        
                                                        my.matrix.null[my.plant.positions,my.plant.positions] <- pp.block
                                                        my.matrix.null[my.herb.positions,my.herb.positions] <- h.block
                                                        my.matrix.null[my.plant.positions,my.herb.positions] <- ph.block
                                                        
                                                        my.matrix.null[my.herb.positions,my.plant.positions] <- t(ph.block)
                                                        
                                                      }else{# all
                                                        pp.block <- reshuffle_matrix(my.matrix.null[my.plant.positions,my.plant.positions])
                                                        fv.block <- reshuffle_matrix(my.matrix.null[my.fv.positions,my.fv.positions])
                                                        pfv.block <- reshuffle_matrix(my.matrix.null[my.plant.positions,my.fv.positions])
                                                        h.block <- reshuffle_matrix(my.matrix.null[my.herb.positions,my.herb.positions])
                                                        ph.block <- reshuffle_matrix(my.matrix.null[my.plant.positions,my.herb.positions])
                                                        
                                                        my.matrix.null[my.plant.positions,my.plant.positions] <- pp.block
                                                        my.matrix.null[my.fv.positions,my.fv.positions] <- fv.block
                                                        my.matrix.null[my.plant.positions,my.fv.positions] <- pfv.block
                                                        my.matrix.null[my.herb.positions,my.herb.positions] <- h.block
                                                        my.matrix.null[my.plant.positions,my.herb.positions] <- ph.block
                                                        
                                                        my.matrix.null[my.herb.positions,my.plant.positions] <- t(ph.block)
                                                        my.matrix.null[my.fv.positions,my.plant.positions] <- t(pfv.block)
                                                        
                                                      }# if-else guild combinations
                                                      
                                                      # -------------------------------------------------------------------------
                                                      # obtain feasibility domains 
                                                      
                                                      A <- -my.matrix.null
                                                      
                                                      my.noise.threshold <- c(0,min(abs(A[which(A != 0)])) * 
                                                                                noise.relative.magnitude)
                                                      
                                                      # check that the randomized matrix is invertible, 
                                                      # otherwise just return NA
                                                      # In the previous code, I add random noise to the matrices
                                                      # to make them invertible, but that is time consuming.
                                                      # this random noise addition is actually integrated in the
                                                      # feasibility_metrics function.
                                                      
                                                      if(nrow(A)>2 & f(A) == "TRUE"){
                                                        sp.exclusions <- try(exclusion_probabilities(A = A,
                                                                                                     omega.replicates = omega.replicates,
                                                                                                     bootstrap.replicates = bootstrap.replicates,
                                                                                                     noise.replicates = noise.replicates,
                                                                                                     noise.threshold = noise.threshold,
                                                                                                     noise.type = noise.type))
                                                      }else{
                                                        sp.exclusions <- data.frame(species = NA,
                                                                                    prob_excl_mean = NA,
                                                                                    prob_excl_lowerCI = NA,
                                                                                    prob_excl_upperCI = NA)
                                                      }
                                                      
                                                      if(inherits(sp.exclusions,"try-error")){
                                                        sp.exclusions <- data.frame(species = NA,
                                                                                    prob_excl_mean = NA,
                                                                                    prob_excl_lowerCI = NA,
                                                                                    prob_excl_upperCI = NA)
                                                      }
                                                      
                                                      sp.exclusions$year <- my.year
                                                      sp.exclusions$plot <- my.plot
                                                      sp.exclusions$guild <- my.guild
                                                      sp.exclusions$intraguild.type <- my.type
                                                      
                                                      null.reps.list[[i.null]] <- sp.exclusions
                                                      
                                                    }
                                                    null.reps.df <- bind_rows(null.reps.list)
                                                    
                                                    # return
                                                    write.csv2(null.reps.df, file = paste("results/null_fd/sp_exclusion_prob_NULL_",id.char[i.id],"_",
                                                                                          null.replicates,"rep_","run",run.number,".csv",sep=""))
                                                    null.reps.df
                                                    
                                                  }

# exc.null.all <- list.files("results/null_fd/",pattern = "sp_exc*",full.names = T) %>% map_dfr(read_csv2)
# exc.null.clean <- exc.null.all %>% 
#   filter(!is.na(species)) %>% 
#   select(-1) %>%
#   select(year,plot,guild,intraguild.type,species,prob_excl_mean,
#          prob_excl_lowerCI,prob_excl_upperCI)
# names(exc.null.clean)[3] <- "guilds"
# 
# write.csv2(exc.null.clean,"results/exclusion_probabilities_null.csv",
#            row.names = F)

# feasibility.df <- feasibility.metrics[[1]]
# exclusions.df <- feasibility.metrics[[2]]

# test for displaying some matrix features, 
# mainly diagonal dominances, diagonal values, 
# which seem to strongly drive fd

# tt <- feasibility.metrics[[1]] %>%
#   group_by(guild,intraguild.type) %>%
#   summarise(mean.zeros = mean(n.zeros),mean.diag = mean(diag.mean),
#             mean.diag.dom = mean(diag.dom))
# ggplot(tt, aes(x = intraguild.type, y = mean.diag)) + 
#   geom_bar(aes(fill = guild),stat = "identity",position=position_dodge()) + 
#   theme_bw() + 
#   NULL

# store results -----------------------------------------------------------
# 
# fd.name <- "feasibility_domain_null"
# fd.name <- paste("results/",fd.name,vers,".csv",sep="")
# 
# # exc.name <- "exclusion_probabilities_observed"
# # exc.name <- paste("results/",exc.name,vers,".csv",sep="")
# 
# write.csv2(x = feasibility.metrics.null,file = fd.name,
#            row.names = FALSE)

stopCluster(cl)



