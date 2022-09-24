
# obtain feasibility domains (and species probabilities of exclusion)
# from community matrices. This code is computationally intensive and 
# it is parallelized. Note that this script obtains the metrics for the 
# observed communities, the following scripts do the same for randomized ones:

# "get_feasibility_domain_null_communities"
# "get_sp_exclusion_prob_null_communities"

# INPUT
# community matrices: "results/community_matrices.RData"
# community names: "results/community_names.RData"

# OUTPUT
# feasibility domain dataframe: "results/feasibility_domain_observed.csv"
# species exclusion prob dataframe: "results/exclusion_probabilities_observed.csv"

# read community matrices -------------------------------------------------
library(tidyverse)
# devtools::install_github("RadicalCommEcol/MultitrophicFun")
library(MultitrophicFun)
library(matlib) # to multiply matrices
library(nleqslv) # to solve non-linear equations
library(zipfR) # beta incomplete function to estimate the area of d-dimensional spherical caps 
library(pracma) # to solve n-dimensional cross products
library(boot) # to bootstrap
library(CValternatives) # to estimate PV index

library(foreach)
library(doParallel)

# source auxiliary functions
list.files("R/feasibility_functions/", full.names = TRUE) %>% map(source)

# set number of cores -----------------------------------------------------

workers <- 10
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
# null.replicates <- length(community_matrices_null[[1]])

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
guild.combinations <- c("plants","floral visitors","herbivores",
                        "plants-floral visitors","plants-herbivores",
                        "all")

# this is the combined id to loop over in parallel
id <- expand.grid(intraguild.types,years,plots,guild.combinations)
id.char <- paste(id[,1],"_",id[,2],"_",id[,3],"_",id[,4],sep="")

# calculate feasibility metrics  -------------------------------------------

comb.fun <- function(...) {
  mapply('rbind', ..., SIMPLIFY=FALSE)
}

feasibility.metrics <- foreach(i.id = 1:length(id.char),
                               .combine=comb.fun, 
                               .packages = c("tidyverse","foreach","matlib",
                                             "nleqslv","zipfR","pracma",
                                             "boot","CValternatives",
                                             "MultitrophicFun")) %dopar% {
                                               
                                               # cat(id.char[i.id],"- started\n")
                                               
                                               list.files("/home/david/Work/Projects/EBD/multitrophic_feasibility/R/feasibility_functions/", 
                                                          full.names = TRUE) %>% map(source)
                                               
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
                                               # obtain feasibility domains and species exclusion probabilities
                                               
                                               A <- -my.matrix
                                               
                                               my.noise.threshold <- c(0,min(abs(A[which(A != 0)])) * 
                                                                         noise.relative.magnitude)
                                               if(nrow(A)>2){
                                                 omega.df <- feasibility_metrics(A = A,
                                                                                 omega.replicates = omega.replicates,
                                                                                 bootstrap.replicates = bootstrap.replicates,
                                                                                 noise.replicates = noise.replicates,
                                                                                 noise.threshold = my.noise.threshold,
                                                                                 noise.type = noise.type)
                                               }else{
                                                 omega.df <- data.frame(omega_mean = NA,
                                                                        omega_lowerCI = NA,
                                                                        omega_upperCI = NA,
                                                                        omega_isotropic = NA,
                                                                        isotropy_index_mean = NA,
                                                                        isotropy_index_lowerCI = NA,
                                                                        isotropy_index_upperCI = NA)
                                               }
                                               # this is a test for checking some features of the matrices
                                               # omega.df <- data.frame(diag.mean = mean(diag(A)),
                                               #                        diag.pos = sum(diag(A) > 0),
                                               #                        diag.dom = mean(diag(A)-rowSums(A)),
                                               #                        n.zeros = sum(A == 0),
                                               #                        total.strength = sum(A))
                                               
                                               omega.df$year <- my.year
                                               omega.df$plot <- my.plot
                                               omega.df$guild <- my.guild
                                               omega.df$intraguild.type <- my.type
                                               
                                               if(nrow(A)>2){
                                                 sp.exclusions <- exclusion_probabilities(A = A,
                                                                                          omega.replicates = omega.replicates,
                                                                                          bootstrap.replicates = bootstrap.replicates,
                                                                                          noise.replicates = noise.replicates,
                                                                                          noise.threshold = noise.threshold,
                                                                                          noise.type = noise.type)
                                               }else{
                                                 sp.exclusions <- data.frame(species = NA,
                                                                             prob_excl_mean = NA,
                                                                             prob_excl_lowerCI = NA,
                                                                             prob_excl_upperCI = NA)
                                               }

                                               sp.exclusions$year <- my.year
                                               sp.exclusions$plot <- my.plot
                                               sp.exclusions$guild <- my.guild
                                               sp.exclusions$intraguild.type <- my.type
                                               
                                               # cat(id.char[i.id],"- completed\n")
                                               # write.csv2(omega.df, file = paste("results/fd_",id.char[i.id],".csv",sep=""))
                                               # write.csv2(sp.exclusions, file = paste("results/sp_exclusions_",id.char[i.id],".csv",sep=""))
                                               
                                               # return
                                               list(omega.df,sp.exclusions)
                                               
                                             }

# store results -----------------------------------------------------------

fd.name <- "feasibility_domain_observed"
fd.name <- paste("results/",fd.name,vers,".csv",sep="")

exc.name <- "exclusion_probabilities_observed"
exc.name <- paste("results/",exc.name,vers,".csv",sep="")

write.csv2(x = feasibility.metrics[[1]],file = fd.name,
           row.names = FALSE)
write.csv2(x = feasibility.metrics[[2]],file = exc.name,
           row.names = FALSE)

stopCluster(cl)



