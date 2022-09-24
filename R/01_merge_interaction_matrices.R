
# combine plant-plant, plant-floral visitor, plant-herbivore matrices
# in a block matrix. Also, calculate intraguild matrices for
# herbivores and floral visitors
# 
# In parallel, store all names of the different species of each guild

# INPUT: 
# plant phenology data: "data/plant_phenology_categories.csv"
# animal phenology data: "data/species_phenology_taxonomy.csv"
# animal resource requirements: "data/species_nest_larval_info.csv"

# OUTPUT:
# list with block community matrices: "results/community_matrices.RData"
# list with species composition in each local community: "results/community_names.RData"

# -------------------------------------------------------------------------

library(tidyverse)
# devtools::install_github("RadicalCommEcol/MultitrophicFun")
library(MultitrophicFun)

source("R/aux_combine_matrices.R")

# -------------------------------------------------------------------------
# add a version suffix?
vers.out <- ""
vers <- ""

# -------------------------------------------------------------------------
# which intraguild matrix types?
intraguild.types <- c("mean_field",
                      "nesting_larvae",
                      "nesting_larvae_phenology"
                      )

# -------------------------------------------------------------------------
# mean field coefficients
mean.field.offdiag <- .2 # Saavedra et al. 2013
mean.field.diag <- 1

# compute null matrices? --------------------------------------------------
# null matrices are obtained for each intraguild matrix type

include.null <- TRUE
replicates <- 100

# read data ---------------------------------------------------------------

years <- c(2019,2020)
plots <- 1:9

plant.phenology <- read.csv2("data/plant_phenology_categories.csv")
animal.phenology <- read.csv2("data/species_phenology_taxonomy.csv")
animal.info <- read.csv2("data/species_nest_larval_info.csv")
animal.nesting.info <- animal.info[,c("ID","nesting")]
animal.larval.info <- animal.info[,c("ID","larval.food.requirements")]

pp.all.years <- list()
ph.all.years <- list()
pfv.all.years <- list()

for(i.year in 1:length(years)){
  
  load(paste("./data/plant_plant_matrices_",
             years[i.year],vers,".RData",sep=""))
  load(paste("./data/plant_floral_visitor_matrices_",
             years[i.year],vers,".RData",sep=""))
  load(paste("./data/plant_herbivore_matrices_",
             years[i.year],vers,".RData",sep=""))
  
  pp.all.years[[i.year]] <- p_p
  ph.all.years[[i.year]] <- p_h
  pfv.all.years[[i.year]] <- p_fv
}

names(pp.all.years) <- years
names(ph.all.years) <- years
names(pfv.all.years) <- years

# remove empty rows and columns -------------------------------------------
# keep only species that appear in a given year and plot

for(i.year in 1:length(years)){
  for(i.plot in 1:length(plots)){
    for(i.guild in c("pp","ph","pfv")){
      
      if(i.guild == "pp"){
        my.matrix <- pp.all.years[[i.year]][[i.plot]]
        
        my.valid.rows <- apply(my.matrix,1,sum)
        my.valid.rows <- names(my.valid.rows)[which(my.valid.rows != 0)]
        
        my.valid.cols <- apply(my.matrix,2,sum)
        my.valid.cols <- names(my.valid.cols)[which(my.valid.cols != 0)]
        
        # slightly different from ph,pfv, because this is a plant-plant
        # matrix, so one sp cannot be only on rows/cols.
        my.valid.sp <- intersect(my.valid.rows,my.valid.cols)
        my.matrix <- my.matrix[my.valid.sp,my.valid.sp]
        
        pp.all.years[[i.year]][[i.plot]] <- my.matrix
        
      }else if(i.guild == "ph"){
        my.matrix <- ph.all.years[[i.year]][[i.plot]]
        
        my.valid.rows <- apply(my.matrix,1,sum)
        my.valid.rows <- names(my.valid.rows)[which(my.valid.rows != 0)]
        
        my.valid.cols <- apply(my.matrix,2,sum)
        my.valid.cols <- names(my.valid.cols)[which(my.valid.cols != 0)]
        
        ph.all.years[[i.year]][[i.plot]] <- my.matrix[my.valid.rows,
                                                      my.valid.cols]
        
      }else if(i.guild == "pfv"){
        my.matrix <- pfv.all.years[[i.year]][[i.plot]]
        
        my.valid.rows <- apply(my.matrix,1,sum)
        my.valid.rows <- names(my.valid.rows)[which(my.valid.rows != 0)]
        
        my.valid.cols <- apply(my.matrix,2,sum)
        my.valid.cols <- names(my.valid.cols)[which(my.valid.cols != 0)]
        
        pfv.all.years[[i.year]][[i.plot]] <- my.matrix[my.valid.rows,
                                                      my.valid.cols]
      }# if i.guild
    }# for i.guild
  }# for i.plot
}# for i.year

# -------------------------------------------------------------------------
# obtain the block matrix associated with each intraguild matrix type

community_matrices <- list()
community_matrices_null <- list()

for(i.type in 1:length(intraguild.types)){
  
  community_matrices[[i.type]] <- list()
  
  # -------------------------------------------------------------------------
  # obtain "observed" matrix of i.type
  cat(i.type," started\n",sep="")
  my.observed.matrix <- aux_combine_matrices(pp.all.years = pp.all.years,
                                             ph.all.years = ph.all.years,
                                             pfv.all.years = pfv.all.years,
                                             plant.phenology = plant.phenology,
                                             animal.phenology = animal.phenology,
                                             animal.nesting.info = animal.nesting.info,
                                             animal.larval.info = animal.larval.info,
                                             randomize = FALSE,
                                             intraguild.type = intraguild.types[i.type],
                                             mean.field.offdiag = mean.field.offdiag,
                                             mean.field.diag = mean.field.diag)

  # -------------------------------------------------------------------------
  community_matrices[[i.type]] <- my.observed.matrix[[1]]
  
  # retrieve the names as well
  # this only needs to be done once, as species composition does not change
  if(i.type == 1){
    sp.names <- my.observed.matrix[[2]]
  }

}# for i.type
names(community_matrices) <- intraguild.types
# -------------------------------------------------------------------------

save(community_matrices,
     file = paste("results/community_matrices",vers.out,".RData",sep=""))

save(sp.names,
     file = paste("results/community_names",vers.out,".RData",sep=""))
