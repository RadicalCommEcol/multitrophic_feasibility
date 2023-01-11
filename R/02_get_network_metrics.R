
# obtain community-level network metrics for the observed and
# randomized communities

# NOTE: to ensure reproducibility, you will have to modify the path inside 
# the foreach loop (lines 482-486)

# INPUT
# community matrices: "results/community_matrices.RData"
# community names: "results/community_names.RData"
# structural zeros: "data/potential_interactions.csv"

# OUTPUT
# community-level metrics: "results/community_metrics.csv" and
# "results/community_metrics_null.csv"

# NOTE:
# in the main analyses we only consider overlap and diagonal dominance,
# but here a few more metrics are obtained. This code is parallelized, and 
# can be computationally demanding.

# If the modularity function raises errors, be sure to have infomap installed, 
# see comment right below.

# -------------------------------------------------------------------------

library(tidyverse)
library(robustbase)
library(igraph)
# devtools::install_github("RadicalCommEcol/MultitrophicFun")
library(MultitrophicFun)

source("R/aux_linkRankModularity.R")
source("R/aux_allowed_interactions.R")

# temp
source("R/interaction_overlap.R")

# see https://ecological-complexity-lab.github.io/infomap_ecology_package/
# Infomap should be installed on the parent folder of the project
library(infomapecology)
source("R/aux_run_infomap_monolayer2.R")
# for parallel computing
source("R/aux_run_infomap_monolayer3.R")
source("R/infomap_modularity.R")
check_infomap()

# -------------------------------------------------------------------------
include.observed <- TRUE

# include null matrices ---------------------------------------------------

include.null <- TRUE
init.replicates <- 1
end.replicates <- 100

# This is a legacy from a code in which I tried different null models
null.models <- c("topology")

# read community matrices -------------------------------------------------

# vers <- "_v2"
vers <- ""

# -------------------------------------------------------------------------

load(file = paste("results/community_matrices",vers,".RData",sep=""))

# load sp names as well
load(file = paste("results/community_names",vers,".RData",sep=""))
plants <- sort(unique(c(sp.names[["2019"]][["plants"]],sp.names[["2020"]][["plants"]])))
fv <- sort(unique(c(sp.names[["2019"]][["floral.visitors"]],sp.names[["2020"]][["floral.visitors"]])))
herb <- sort(unique(c(sp.names[["2019"]][["herbivores"]],sp.names[["2020"]][["herbivores"]])))

# read structural zeros ---------------------------------------------------

sz <- read.csv2(paste("data/potential_interactions",vers,".csv",sep=""))

# sz are divided by year, so join them together in a single set
# such that pairs that are structural zero all years are kept, and
# those that have at least one observation in one year are not structural zeros

sz.combined <- sz %>% pivot_wider(names_from = year,values_from = c(structural.zero,observed.freq))

sz.combined$structural.zero_2019[which(is.na(sz.combined$structural.zero_2019))] <- TRUE
sz.combined$structural.zero_2020[which(is.na(sz.combined$structural.zero_2020))] <- TRUE
sz.combined$structural.zero <- sz.combined$structural.zero_2019 * sz.combined$structural.zero_2020

sz.combined$observed.freq_2019[which(is.na(sz.combined$observed.freq_2019))] <- 0
sz.combined$observed.freq_2020[which(is.na(sz.combined$observed.freq_2020))] <- 0
sz.combined$observed.freq <- sz.combined$observed.freq_2019 + sz.combined$observed.freq_2020
sz.combined$structural.zero[sz.combined$observed.freq > 0] <- 0

sz.combined <- sz.combined[,c("plant","animal","type","overlap","observed.freq","structural.zero")]
sz.combined$structural.zero <- as.logical(sz.combined$structural.zero)

# define range ------------------------------------------------------------

years <- c(2019,2020)
plots <- 1:9

guild.combinations <- c("plants",
                        "floral visitors",
                        "herbivores",
                        "plants-floral visitors",
                        "plants-herbivores",
                        "all")

intraguild.types <- names(community_matrices)

# define metrics ----------------------------------------------------------

# constants for the modularity algorithm
lr.damping <- .85
lr.alg <- "prpack"

# minimum number of nodes to calculate modularity
min.node.modularity <- 2

metric.names <- c("richness",
                  "connectance",
                  "diagonally_dominant_sp",
                  "avg_diagonal_dominance",
                  "avg_intraguild_niche_overlap",
                  "avg_interguild_niche_overlap",
                  "degree_distribution",
                  "qual_modularity",
                  "complexity")

community.metrics <- expand.grid(year = years, 
                                 plot = plots, 
                                 guilds = guild.combinations,
                                 intraguild.type = intraguild.types,
                                 metric = metric.names,value = NA)

module_members <- NULL

# -------------------------------------------------------------------------
# auxiliary function for calculating community-level overlap
# from a list of pairwise values
# this calculates the average total overlap, i.e. the total overlap of 
# a species with the rest of the community, averaged.

avg_total_overlap <- function(df){
  res <- df %>% 
    group_by(sp1) %>%
    summarise(total = sum(overlap, na.rm = T))
  return(mean(res$total))
}


# -------------------------------------------------------------------------
# TESTs
# i.type <- 1
# i.year <- 1
# i.plot <- 1
# i.guild <- 3
# 
# my.type <- "mean_field"
# my.year <- "2019"
# my.plot <- 1
# my.guild <- "herbivores"
# my.null <- "topology"
# my.rep <- 1

# calculate metrics -------------------------------------------------------
if(include.observed){
  
  for(i.type in 1:length(intraguild.types)){
    for(i.year in 1:length(years)){
      for(i.plot in plots){
        
        year.plot.matrix <- 
          community_matrices[[i.type]][[as.character(years[i.year])]][[i.plot]]
        
        plant.positions <- which(rownames(year.plot.matrix) %in% 
                                   sp.names[[i.year]][["plants"]])
        herb.positions <- which(rownames(year.plot.matrix) %in% 
                                  sp.names[[i.year]][["herbivores"]])
        fv.positions <- which(rownames(year.plot.matrix) %in% 
                                sp.names[[i.year]][["floral.visitors"]])
        
        for(i.guild in 1:length(guild.combinations)){
          
          cat(intraguild.types[i.type],"-",years[i.year],"-",i.plot,"-",guild.combinations[i.guild],"\n")
          
          if(guild.combinations[i.guild] == "plants"){
            
            my.matrix <- year.plot.matrix[plant.positions,plant.positions]
            
          }else if(guild.combinations[i.guild] == "plants-floral visitors"){
            
            my.matrix <- year.plot.matrix[c(plant.positions,fv.positions),
                                          c(plant.positions,fv.positions)]
            
          }else if(guild.combinations[i.guild] == "floral visitors"){
            
            my.matrix <- year.plot.matrix[fv.positions,fv.positions]
            
          }else if(guild.combinations[i.guild] == "herbivores"){
            
            my.matrix <- year.plot.matrix[herb.positions,herb.positions]
            
          }else if(guild.combinations[i.guild] == "plants-herbivores"){
            
            my.matrix <- year.plot.matrix[c(plant.positions,herb.positions),
                                          c(plant.positions,herb.positions)]
            
          }else if(guild.combinations[i.guild] == "all"){
            
            my.matrix <- year.plot.matrix
          }
          
          # convert to igraph object ------------------------------------------------
          
          my.network <- graph.adjacency(my.matrix,
                                        mode="directed",
                                        weighted=TRUE)
          
          # richness ----------------------------------------------------------------
          
          my.richness <- length(unique(c(rownames(my.matrix),colnames(my.matrix))))
          
          # diagonal dominance ------------------------------------------------------
          
          # diagonal elements
          d <- diag(my.matrix)
          # rowsums
          nd <- rowSums(my.matrix)
          # but without the diagonal
          nd <- nd - d
          
          # number of diagonally dominant sp
          my.diag.dom.sp <- sum(abs(d) - abs(nd) > 0)
          
          # average degree of diagonal dominance
          my.avg.diag.dom <- mean(abs(d) - abs(nd))
          
          # niche overlap -----------------------------------------------------------
          # differentiate intra and inter-guild overlap
          
          if(guild.combinations[i.guild] %in% c("plants","floral visitors","herbivores")){
            
            # in these cases, only intra-guild overlap
            pair.intra.overlap <- interaction_overlap(abs(my.matrix))
            ov.clean <- subset(pair.intra.overlap, sp1 != sp2 & !is.nan(overlap) & !is.na(overlap))
            
            my.avg.intra.overlap <- avg_total_overlap(ov.clean)
            my.avg.inter.overlap <- NA
            # # from https://stackoverflow.com/questions/32669609/identifying-unique-pairs-of-values-from-two-columns-in-a-dataframe
            # ov.unique <-  ov.clean[!duplicated(t(apply(ov.clean, 1, sort))),]
            # 
            # my.avg.intra.overlap <- mean(ov.unique$overlap,na.rm = T)
            # my.avg.inter.overlap <- NA
            
          }else if(guild.combinations[i.guild] == "plants-floral visitors"){
            
            my.intra.plants <- year.plot.matrix[plant.positions,plant.positions]
            
            intra.plants.overlap <- interaction_overlap(abs(my.intra.plants))
            intra.plants.clean <- subset(intra.plants.overlap, sp1 != sp2 & !is.nan(overlap) & !is.na(overlap))
            intra.plants.unique <-  intra.plants.clean[!duplicated(t(apply(intra.plants.clean, 1, sort))),]
            
            my.intra.fv <- year.plot.matrix[fv.positions,fv.positions]
            
            intra.fv.overlap <- interaction_overlap(abs(my.intra.fv))
            intra.fv.clean <- subset(intra.fv.overlap, sp1 != sp2 & !is.nan(overlap) & !is.na(overlap))
            intra.fv.unique <-  intra.fv.clean[!duplicated(t(apply(intra.fv.clean, 1, sort))),]
            
            my.inter.plants.fv <- year.plot.matrix[plant.positions,fv.positions]
            
            inter.plants.fv.overlap <- interaction_overlap(abs(my.inter.plants.fv))
            inter.plants.fv.clean <- subset(inter.plants.fv.overlap, sp1 != sp2 & !is.nan(overlap) & !is.na(overlap))
            inter.plants.fv.unique <-  inter.plants.fv.clean[!duplicated(t(apply(inter.plants.fv.clean, 1, sort))),]
            
            my.inter.fv.plants <- year.plot.matrix[fv.positions,plant.positions]
            
            inter.fv.plants.overlap <- interaction_overlap(abs(my.inter.fv.plants))
            inter.fv.plants.clean <- subset(inter.fv.plants.overlap, sp1 != sp2 & !is.nan(overlap) & !is.na(overlap))
            inter.fv.plants.unique <-  inter.fv.plants.clean[!duplicated(t(apply(inter.fv.plants.clean, 1, sort))),]
            
            my.avg.intra.overlap <- avg_total_overlap(bind_rows(intra.plants.clean,
                                                                intra.fv.clean))
            my.avg.inter.overlap <- avg_total_overlap(bind_rows(inter.plants.fv.clean,
                                                                inter.fv.plants.clean))
            
          }else if(guild.combinations[i.guild] == "plants-herbivores"){
            
            my.intra.plants <- year.plot.matrix[plant.positions,plant.positions]
            
            intra.plants.overlap <- interaction_overlap(abs(my.intra.plants))
            intra.plants.clean <- subset(intra.plants.overlap, sp1 != sp2 & !is.nan(overlap) & !is.na(overlap))
            intra.plants.unique <-  intra.plants.clean[!duplicated(t(apply(intra.plants.clean, 1, sort))),]
            
            my.intra.herb <- year.plot.matrix[herb.positions,herb.positions]
            
            intra.herb.overlap <- interaction_overlap(abs(my.intra.herb))
            intra.herb.clean <- subset(intra.herb.overlap, sp1 != sp2 & !is.nan(overlap) & !is.na(overlap))
            intra.herb.unique <-  intra.herb.clean[!duplicated(t(apply(intra.herb.clean, 1, sort))),]
            
            my.inter.plants.herb <- year.plot.matrix[plant.positions,herb.positions]
            
            inter.plants.herb.overlap <- interaction_overlap(abs(my.inter.plants.herb))
            inter.plants.herb.clean <- subset(inter.plants.herb.overlap, sp1 != sp2 & !is.nan(overlap) & !is.na(overlap))
            inter.plants.herb.unique <-  inter.plants.herb.clean[!duplicated(t(apply(inter.plants.herb.clean, 1, sort))),]
            
            my.inter.herb.plants <- year.plot.matrix[herb.positions,plant.positions]
            
            inter.herb.plants.overlap <- interaction_overlap(abs(my.inter.herb.plants))
            inter.herb.plants.clean <- subset(inter.herb.plants.overlap, sp1 != sp2 & !is.nan(overlap) & !is.na(overlap))
            inter.herb.plants.unique <-  inter.herb.plants.clean[!duplicated(t(apply(inter.herb.plants.clean, 1, sort))),]
            
            my.avg.intra.overlap <- avg_total_overlap(bind_rows(intra.plants.clean,
                                                                intra.herb.clean))
            my.avg.inter.overlap <- avg_total_overlap(bind_rows(inter.plants.herb.clean,
                                                                inter.herb.plants.clean))
            
          }else{
            
            my.intra.plants <- year.plot.matrix[plant.positions,plant.positions]
            
            intra.plants.overlap <- interaction_overlap(abs(my.intra.plants))
            intra.plants.clean <- subset(intra.plants.overlap, sp1 != sp2 & !is.nan(overlap) & !is.na(overlap))
            intra.plants.unique <-  intra.plants.clean[!duplicated(t(apply(intra.plants.clean, 1, sort))),]
            
            my.intra.fv <- year.plot.matrix[fv.positions,fv.positions]
            
            intra.fv.overlap <- interaction_overlap(abs(my.intra.fv))
            intra.fv.clean <- subset(intra.fv.overlap, sp1 != sp2 & !is.nan(overlap) & !is.na(overlap))
            intra.fv.unique <-  intra.fv.clean[!duplicated(t(apply(intra.fv.clean, 1, sort))),]
            
            my.intra.herb <- year.plot.matrix[herb.positions,herb.positions]
            
            intra.herb.overlap <- interaction_overlap(abs(my.intra.herb))
            intra.herb.clean <- subset(intra.herb.overlap, sp1 != sp2 & !is.nan(overlap) & !is.na(overlap))
            intra.herb.unique <-  intra.herb.clean[!duplicated(t(apply(intra.herb.clean, 1, sort))),]
            
            my.inter.plants.fv <- year.plot.matrix[plant.positions,fv.positions]
            
            inter.plants.fv.overlap <- interaction_overlap(abs(my.inter.plants.fv))
            inter.plants.fv.clean <- subset(inter.plants.fv.overlap, sp1 != sp2 & !is.nan(overlap) & !is.na(overlap))
            inter.plants.fv.unique <-  inter.plants.fv.clean[!duplicated(t(apply(inter.plants.fv.clean, 1, sort))),]
            
            my.inter.fv.plants <- year.plot.matrix[fv.positions,plant.positions]
            
            inter.fv.plants.overlap <- interaction_overlap(abs(my.inter.fv.plants))
            inter.fv.plants.clean <- subset(inter.fv.plants.overlap, sp1 != sp2 & !is.nan(overlap) & !is.na(overlap))
            inter.fv.plants.unique <-  inter.fv.plants.clean[!duplicated(t(apply(inter.fv.plants.clean, 1, sort))),]
            
            my.inter.plants.herb <- year.plot.matrix[plant.positions,herb.positions]
            
            inter.plants.herb.overlap <- interaction_overlap(abs(my.inter.plants.herb))
            inter.plants.herb.clean <- subset(inter.plants.herb.overlap, sp1 != sp2 & !is.nan(overlap) & !is.na(overlap))
            inter.plants.herb.unique <-  inter.plants.herb.clean[!duplicated(t(apply(inter.plants.herb.clean, 1, sort))),]
            
            my.inter.herb.plants <- year.plot.matrix[herb.positions,plant.positions]
            
            inter.herb.plants.overlap <- interaction_overlap(abs(my.inter.herb.plants))
            inter.herb.plants.clean <- subset(inter.herb.plants.overlap, sp1 != sp2 & !is.nan(overlap) & !is.na(overlap))
            inter.herb.plants.unique <-  inter.herb.plants.clean[!duplicated(t(apply(inter.herb.plants.clean, 1, sort))),]
            
            my.avg.intra.overlap <- avg_total_overlap(bind_rows(intra.plants.clean,
                                                                intra.fv.clean,
                                                                intra.herb.clean))
            my.avg.inter.overlap <- avg_total_overlap(bind_rows(inter.plants.herb.clean,
                                                                inter.herb.plants.clean,
                                                                inter.plants.fv.clean,
                                                                inter.fv.plants.clean))
            
          }
          
          # connectance -------------------------------------------------------------
          # for now, binary
          
          sz.all <- get_structural_zeros(A = my.matrix,
                                         plants = plants, 
                                         fv = fv, 
                                         herb = herb,
                                         dfz = sz.combined)
          
          my.connectance <- connectance(interaction.matrix = my.matrix,
                                        quant = FALSE,
                                        structural_zeros = sz.all)
          
          # degree distribution -----------------------------------------------------
          
          my.degrees <- igraph::degree(my.network)
          # first take: gini index
          # my.deg.dist <- ineq::ineq(my.degrees,type = "Gini")
          
          # second take - shannon entropy, as referenced in Goñi et al. 2013
          my.deg.dist <- vegan::diversity(my.degrees, index = "shannon")
          
          # May's complexity --------------------------------------------------------
          
          # connectance explicitly deals with structural zeros
          # so sigma must, as well.
          
          all.interactions <- expand.grid(row = 1:nrow(my.matrix),
                                          col = 1:ncol(my.matrix))
          all.interactions$value <- c(my.matrix)
          
          # return all interactions wihout a match in sz
          valid.interactions <- anti_join(all.interactions,sz.all)
          
          # sd of valid interaction strengths
          my.sigma <- sd(valid.interactions$value)
          
          my.complexity <- MultitrophicFun::May_complexity(S = nrow(my.matrix),
                                                           C = my.connectance,
                                                           sigma = my.sigma)
          
          # modularity --------------------------------------------------------------
          
          if(nrow(my.matrix) >= min.node.modularity){
            
            bin.matrix <- my.matrix
            bin.matrix[bin.matrix != 0] <- 1
            linkrank_binary_modularity <- infomap_modularity(A = bin.matrix,
                                                             allowed.interactions.mat = matrix(1,nrow = nrow(my.matrix),
                                                                                               ncol = ncol(my.matrix)),
                                                             damping = lr.damping,
                                                             pr.algo = lr.alg)
            
          }else{ # less than 5 sp
            linkrank_binary_modularity <- NA_real_
          }
          # store -------------------------------------------------------------------
          
          pos <- which(community.metrics$intraguild.type == intraguild.types[i.type] &
                         community.metrics$year == years[i.year] &
                         community.metrics$plot == i.plot &
                         community.metrics$guilds == guild.combinations[i.guild])
          
          community.metrics$value[pos[which(community.metrics$metric[pos] == "richness")]] <- my.richness
          community.metrics$value[pos[which(community.metrics$metric[pos] == "connectance")]] <- my.connectance
          community.metrics$value[pos[which(community.metrics$metric[pos] == "diagonally_dominant_sp")]] <- my.diag.dom.sp
          community.metrics$value[pos[which(community.metrics$metric[pos] == "avg_diagonal_dominance")]] <- my.avg.diag.dom
          community.metrics$value[pos[which(community.metrics$metric[pos] == "avg_intraguild_niche_overlap")]] <- my.avg.intra.overlap
          community.metrics$value[pos[which(community.metrics$metric[pos] == "avg_interguild_niche_overlap")]] <- my.avg.inter.overlap
          community.metrics$value[pos[which(community.metrics$metric[pos] == "degree_distribution")]] <- my.deg.dist
          community.metrics$value[pos[which(community.metrics$metric[pos] == "complexity")]] <- my.complexity
          community.metrics$value[pos[which(community.metrics$metric[pos] == "qual_modularity")]] <- linkrank_binary_modularity
          
        }# for i.guild
      }# for i.plot
    }# for i.year
  }# for i.type
  
  # store metrics -----------------------------------------------------------
  
  # year - plot - guilds - metric - value
  # store a different file for observed and null networks
  
  write.csv2(community.metrics,file = paste("results/community_metrics",vers,".csv",sep=""),
             row.names = FALSE)
}# if include.observed

# repeat for null matrices ------------------------------------------------

if(include.null){
  
  # -------------------------------------------------------------------------
  library(foreach)
  library(doParallel)
  workers <- 4
  cl <- makeCluster(workers)
  # register the cluster for using foreach
  registerDoParallel(cl)
  
  # results dataframe -------------------------------------------------------
  
  community.null.metrics <- list()
  
  id <- expand.grid(intraguild.types,years,plots,guild.combinations,null.models,sprintf("%04d",init.replicates:end.replicates))
  id.char <- paste(id[,1],"_",id[,2],"_",id[,3],"_",id[,4],"_",id[,5],"_",id[,6],sep="")

  community.null.metrics <- foreach(i.id = 1:length(id.char),
                                    # .combine=comb.fun, 
                                    .packages = c("tidyverse",
                                                  "infomapecology",
                                                  "robustbase","igraph",
                                                  "MultitrophicFun")) %dopar% {
                                                    
                                                    source("/home/david/Work/Projects/EBD/multitrophic_feasibility/R/aux_linkRankModularity.R")
                                                    source("/home/david/Work/Projects/EBD/multitrophic_feasibility/R/aux_run_infomap_monolayer3.R")
                                                    source("/home/david/Work/Projects/EBD/multitrophic_feasibility/R/infomap_modularity.R")
                                                    source("/home/david/Work/Projects/EBD/multitrophic_feasibility/R/reshuffle_matrix.R")
                                                    source("/home/david/Work/Projects/EBD/multitrophic_feasibility/R/interaction_overlap.R")
                                                    
                                                    # first, recover the parameters of each matrix
                                                    
                                                    my.type <- NA
                                                    my.year <- NA
                                                    my.plot <- NA
                                                    my.guild <- NA
                                                    my.null <- NA
                                                    my.rep <- NA
                                                    
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
                                                    
                                                    # null model
                                                    if(grepl("topology",id.char[i.id])){
                                                      my.null <- "topology"
                                                    }else{
                                                      my.null <- "diagonal dominance"
                                                    }
                                                    
                                                    # replicate
                                                    my.rep <- substr(id.char[i.id],nchar(id.char[i.id])-3,nchar(id.char[i.id]))
                                                    
                                                    # plot
                                                    my.plot1 <- stringr::str_remove(id.char[i.id],my.type)
                                                    my.plot2 <- stringr::str_remove(my.plot1,my.guild)
                                                    my.plot3 <- stringr::str_remove(my.plot2,my.year)
                                                    my.plot4 <- stringr::str_remove(my.plot3,my.null)
                                                    my.plot5 <- stringr::str_remove(my.plot4,my.rep)
                                                    my.plot <- as.numeric(stringr::str_remove_all(my.plot5,"_"))
                                                    
                                                    year.plot.matrix <- 
                                                      community_matrices[[my.type]][[as.character(my.year)]][[my.plot]]
                                                    
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
                                                    
                                                      # this reshuffles interactions completely within any given block.
                                                      # it breaks topology and diagonal dominance. It is not possible 
                                                      # to keep diagonal dominance while breaking topology, because
                                                      # breaking topology in a block (say plant-fv) but keeping
                                                      # rowsums, necessarily translates in keeping colsums in the 
                                                      # transposed block (fv-plant), thus breaking diagonal dominance
                                                      
                                                      # therefore, it makes more sense to apply an incremental approach:
                                                      # in the previous null, topology is maintained, and here, 
                                                      # it is not.
                                                      
                                                      my.matrix.null <- my.matrix
                                                      
                                                      # this is to make sure that no interactions are assigned
                                                      # to pollinators-herbivores blocks. In fact,
                                                      # interactions are reshuffled only within their block.
                                                      
                                                      if(my.guild %in% c("plants",
                                                                         "floral visitors",
                                                                         "herbivores")){
                                                        
                                                        my.matrix.null <- reshuffle_matrix(my.matrix)
                                                        # for(i.row in 1:nrow(my.matrix.null)){
                                                        #   my.matrix.null[i.row,] <- reshuffle.non.diag(v = my.matrix.null[i.row,],
                                                        #                                                d = diag(my.matrix.null)[i.row])
                                                        # }# for i.row
                                                        
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
                                                    
                                                    # convert to igraph object ------------------------------------------------
                                                    
                                                    my.network <- graph.adjacency(my.matrix.null,
                                                                                  mode="directed",
                                                                                  weighted=TRUE)
                                                    
                                                    # richness ----------------------------------------------------------------
                                                    
                                                    my.richness <- length(unique(c(rownames(my.matrix.null),colnames(my.matrix.null))))
                                                    
                                                    # diagonal dominance ------------------------------------------------------
                                                    
                                                    # diagonal elements
                                                    d <- diag(my.matrix.null)
                                                    # rowsums
                                                    nd <- rowSums(my.matrix.null)
                                                    # but without the diagonal
                                                    nd <- nd - d
                                                    
                                                    # number of diagonally dominant sp
                                                    my.diag.dom.sp <- sum(abs(d) - abs(nd) > 0)
                                                    
                                                    # average degree of diagonal dominance
                                                    my.avg.diag.dom <- mean(abs(d) - abs(nd))
                                                    
                                                    # niche overlap -----------------------------------------------------------
                                                    # differentiate intra and inter-guild overlap
                                                    
                                                    if(my.guild %in% c("plants","floral visitors","herbivores")){
                                                      
                                                      # in these cases, only intra-guild overlap
                                                      pair.intra.overlap <- interaction_overlap(abs(my.matrix.null))
                                                      ov.clean <- subset(pair.intra.overlap, sp1 != sp2 & !is.nan(overlap) & !is.na(overlap))
                                                      # from https://stackoverflow.com/questions/32669609/identifying-unique-pairs-of-values-from-two-columns-in-a-dataframe
                                                      ov.unique <-  ov.clean[!duplicated(t(apply(ov.clean, 1, sort))),]
                                                      
                                                      my.avg.intra.overlap <- avg_total_overlap(ov.clean)
                                                      my.avg.inter.overlap <- NA
                                                      
                                                      # my.avg.intra.overlap <- mean(ov.unique$overlap,na.rm = T)
                                                      # my.avg.inter.overlap <- NA
                                                      
                                                    }else if(my.guild == "plants-floral visitors"){
                                                      
                                                      my.intra.plants <- my.matrix.null[my.plant.positions,my.plant.positions]
                                                      
                                                      intra.plants.overlap <- interaction_overlap(abs(my.intra.plants))
                                                      intra.plants.clean <- subset(intra.plants.overlap, sp1 != sp2 & !is.nan(overlap) & !is.na(overlap))
                                                      intra.plants.unique <-  intra.plants.clean[!duplicated(t(apply(intra.plants.clean, 1, sort))),]
                                                      
                                                      my.intra.fv <- my.matrix.null[my.fv.positions,my.fv.positions]
                                                      
                                                      intra.fv.overlap <- interaction_overlap(abs(my.intra.fv))
                                                      intra.fv.clean <- subset(intra.fv.overlap, sp1 != sp2 & !is.nan(overlap) & !is.na(overlap))
                                                      intra.fv.unique <-  intra.fv.clean[!duplicated(t(apply(intra.fv.clean, 1, sort))),]
                                                      
                                                      my.inter.plants.fv <- my.matrix.null[my.plant.positions,my.fv.positions]
                                                      
                                                      inter.plants.fv.overlap <- interaction_overlap(abs(my.inter.plants.fv))
                                                      inter.plants.fv.clean <- subset(inter.plants.fv.overlap, sp1 != sp2 & !is.nan(overlap) & !is.na(overlap))
                                                      inter.plants.fv.unique <-  inter.plants.fv.clean[!duplicated(t(apply(inter.plants.fv.clean, 1, sort))),]
                                                      
                                                      my.inter.fv.plants <- my.matrix.null[my.fv.positions,my.plant.positions]
                                                      
                                                      inter.fv.plants.overlap <- interaction_overlap(abs(my.inter.fv.plants))
                                                      inter.fv.plants.clean <- subset(inter.fv.plants.overlap, sp1 != sp2 & !is.nan(overlap) & !is.na(overlap))
                                                      inter.fv.plants.unique <-  inter.fv.plants.clean[!duplicated(t(apply(inter.fv.plants.clean, 1, sort))),]
                                                      
                                                      my.avg.intra.overlap <- avg_total_overlap(bind_rows(intra.plants.clean,
                                                                                                          intra.fv.clean))
                                                      my.avg.inter.overlap <- avg_total_overlap(bind_rows(inter.plants.fv.clean,
                                                                                                          inter.fv.plants.clean))
                                                      
                                                      # my.avg.intra.overlap <- mean(c(intra.plants.unique$overlap,
                                                      #                                intra.fv.unique$overlap))
                                                      # my.avg.inter.overlap <- mean(c(inter.plants.fv.unique$overlap,
                                                      #                                inter.fv.plants.unique$overlap))
                                                      
                                                    }else if(my.guild == "plants-herbivores"){
                                                      
                                                      my.intra.plants <- my.matrix.null[my.plant.positions,my.plant.positions]
                                                      
                                                      intra.plants.overlap <- interaction_overlap(abs(my.intra.plants))
                                                      intra.plants.clean <- subset(intra.plants.overlap, sp1 != sp2 & !is.nan(overlap) & !is.na(overlap))
                                                      intra.plants.unique <-  intra.plants.clean[!duplicated(t(apply(intra.plants.clean, 1, sort))),]
                                                      
                                                      my.intra.herb <- my.matrix.null[my.herb.positions,my.herb.positions]
                                                      
                                                      intra.herb.overlap <- interaction_overlap(abs(my.intra.herb))
                                                      intra.herb.clean <- subset(intra.herb.overlap, sp1 != sp2 & !is.nan(overlap) & !is.na(overlap))
                                                      intra.herb.unique <-  intra.herb.clean[!duplicated(t(apply(intra.herb.clean, 1, sort))),]
                                                      
                                                      my.inter.plants.herb <- my.matrix.null[my.plant.positions,my.herb.positions]
                                                      
                                                      inter.plants.herb.overlap <- interaction_overlap(abs(my.inter.plants.herb))
                                                      inter.plants.herb.clean <- subset(inter.plants.herb.overlap, sp1 != sp2 & !is.nan(overlap) & !is.na(overlap))
                                                      inter.plants.herb.unique <-  inter.plants.herb.clean[!duplicated(t(apply(inter.plants.herb.clean, 1, sort))),]
                                                      
                                                      my.inter.herb.plants <- my.matrix.null[my.herb.positions,my.plant.positions]
                                                      
                                                      inter.herb.plants.overlap <- interaction_overlap(abs(my.inter.herb.plants))
                                                      inter.herb.plants.clean <- subset(inter.herb.plants.overlap, sp1 != sp2 & !is.nan(overlap) & !is.na(overlap))
                                                      inter.herb.plants.unique <-  inter.herb.plants.clean[!duplicated(t(apply(inter.herb.plants.clean, 1, sort))),]
                                                      
                                                      my.avg.intra.overlap <- avg_total_overlap(bind_rows(intra.plants.clean,
                                                                                                          intra.herb.clean))
                                                      my.avg.inter.overlap <- avg_total_overlap(bind_rows(inter.plants.herb.clean,
                                                                                                          inter.herb.plants.clean))
                                                      
                                                    }else{
                                                      
                                                      my.intra.plants <- my.matrix.null[my.plant.positions,my.plant.positions]
                                                      
                                                      intra.plants.overlap <- interaction_overlap(abs(my.intra.plants))
                                                      intra.plants.clean <- subset(intra.plants.overlap, sp1 != sp2 & !is.nan(overlap) & !is.na(overlap))
                                                      intra.plants.unique <-  intra.plants.clean[!duplicated(t(apply(intra.plants.clean, 1, sort))),]
                                                      
                                                      my.intra.fv <- my.matrix.null[my.fv.positions,my.fv.positions]
                                                      
                                                      intra.fv.overlap <- interaction_overlap(abs(my.intra.fv))
                                                      intra.fv.clean <- subset(intra.fv.overlap, sp1 != sp2 & !is.nan(overlap) & !is.na(overlap))
                                                      intra.fv.unique <-  intra.fv.clean[!duplicated(t(apply(intra.fv.clean, 1, sort))),]
                                                      
                                                      my.intra.herb <- my.matrix.null[my.herb.positions,my.herb.positions]
                                                      
                                                      intra.herb.overlap <- interaction_overlap(abs(my.intra.herb))
                                                      intra.herb.clean <- subset(intra.herb.overlap, sp1 != sp2 & !is.nan(overlap) & !is.na(overlap))
                                                      intra.herb.unique <-  intra.herb.clean[!duplicated(t(apply(intra.herb.clean, 1, sort))),]
                                                      
                                                      my.inter.plants.fv <- my.matrix.null[my.plant.positions,my.fv.positions]
                                                      
                                                      inter.plants.fv.overlap <- interaction_overlap(abs(my.inter.plants.fv))
                                                      inter.plants.fv.clean <- subset(inter.plants.fv.overlap, sp1 != sp2 & !is.nan(overlap) & !is.na(overlap))
                                                      inter.plants.fv.unique <-  inter.plants.fv.clean[!duplicated(t(apply(inter.plants.fv.clean, 1, sort))),]
                                                      
                                                      my.inter.fv.plants <- my.matrix.null[my.fv.positions,my.plant.positions]
                                                      
                                                      inter.fv.plants.overlap <- interaction_overlap(abs(my.inter.fv.plants))
                                                      inter.fv.plants.clean <- subset(inter.fv.plants.overlap, sp1 != sp2 & !is.nan(overlap) & !is.na(overlap))
                                                      inter.fv.plants.unique <-  inter.fv.plants.clean[!duplicated(t(apply(inter.fv.plants.clean, 1, sort))),]
                                                      
                                                      my.inter.plants.herb <- my.matrix.null[my.plant.positions,my.herb.positions]
                                                      
                                                      inter.plants.herb.overlap <- interaction_overlap(abs(my.inter.plants.herb))
                                                      inter.plants.herb.clean <- subset(inter.plants.herb.overlap, sp1 != sp2 & !is.nan(overlap) & !is.na(overlap))
                                                      inter.plants.herb.unique <-  inter.plants.herb.clean[!duplicated(t(apply(inter.plants.herb.clean, 1, sort))),]
                                                      
                                                      my.inter.herb.plants <- my.matrix.null[my.herb.positions,my.plant.positions]
                                                      
                                                      inter.herb.plants.overlap <- interaction_overlap(abs(my.inter.herb.plants))
                                                      inter.herb.plants.clean <- subset(inter.herb.plants.overlap, sp1 != sp2 & !is.nan(overlap) & !is.na(overlap))
                                                      inter.herb.plants.unique <-  inter.herb.plants.clean[!duplicated(t(apply(inter.herb.plants.clean, 1, sort))),]
                                                      
                                                      my.avg.intra.overlap <- avg_total_overlap(bind_rows(intra.plants.clean,
                                                                                                          intra.fv.clean,
                                                                                                          intra.herb.clean))
                                                      my.avg.inter.overlap <- avg_total_overlap(bind_rows(inter.plants.herb.clean,
                                                                                                          inter.herb.plants.clean,
                                                                                                          inter.plants.fv.clean,
                                                                                                          inter.fv.plants.clean))
                                                      
                                                    }
                                                    
                                                    # connectance -------------------------------------------------------------
                                                    # for now, binary
                                                    
                                                    sz.all <- get_structural_zeros(A = my.matrix.null,
                                                                                   plants = plants, 
                                                                                   fv = fv, 
                                                                                   herb = herb,
                                                                                   dfz = sz.combined)
                                                    
                                                    my.connectance <- connectance(interaction.matrix = my.matrix.null,
                                                                                  quant = FALSE,
                                                                                  structural_zeros = sz.all)
                                                    
                                                    # degree distribution -----------------------------------------------------
                                                    
                                                    my.degrees <- igraph::degree(my.network)
                                                    # first take: gini index
                                                    # my.deg.dist <- ineq::ineq(my.degrees,type = "Gini")
                                                    
                                                    # second take - shannon entropy, as referenced in Goñi et al. 2013
                                                    my.deg.dist <- vegan::diversity(my.degrees, index = "shannon")
                                                    
                                                    # May's complexity --------------------------------------------------------
                                                    
                                                    # connectance explicitly deals with structural zeros
                                                    # so sigma must, as well.
                                                    
                                                    all.interactions <- expand.grid(row = 1:nrow(my.matrix.null),
                                                                                    col = 1:ncol(my.matrix.null))
                                                    all.interactions$value <- c(my.matrix.null)
                                                    
                                                    # return all interactions wihout a match in sz
                                                    valid.interactions <- anti_join(all.interactions,sz.all)
                                                    
                                                    # sd of valid interaction strengths
                                                    my.sigma <- sd(valid.interactions$value)
                                                    
                                                    my.complexity <- MultitrophicFun::May_complexity(S = nrow(my.matrix.null),
                                                                                                     C = my.connectance,
                                                                                                     sigma = my.sigma)
                                                    
                                                    # modularity --------------------------------------------------------------
                                                    
                                                    if(nrow(my.matrix.null) >= min.node.modularity){
                                                      bin.matrix.null <- my.matrix.null
                                                      bin.matrix.null[bin.matrix.null != 0] <- 1
                                                      linkrank_binary_modularity <- infomap_modularity(A = bin.matrix.null,
                                                                                                       allowed.interactions.mat = matrix(1,nrow = nrow(my.matrix.null),
                                                                                                                                         ncol = ncol(my.matrix.null)),
                                                                                                       temp.dir = paste("R/temp/d",i.id,sep=""),
                                                                                                       damping = lr.damping,
                                                                                                       pr.algo = lr.alg)
                                                      
                                                    }else{ # less than 5 sp
                                                      linkrank_binary_modularity <- NA_real_
                                                    }
                                                    # store -------------------------------------------------------------------
                                                    
                                                    cm <- data.frame(intraguild.type = my.type,
                                                                     year = my.year,
                                                                     plot = my.plot,
                                                                     guilds = my.guild,
                                                                     null.model = my.null,
                                                                     replicate = my.rep,
                                                                     richness = my.richness,
                                                                     connectance = my.connectance,
                                                                     diagonally_dominant_sp = my.diag.dom.sp,
                                                                     avg_diagonal_dominance = my.avg.diag.dom,
                                                                     avg_intraguild_niche_overlap = my.avg.intra.overlap,
                                                                     avg_interguild_niche_overlap = my.avg.inter.overlap,
                                                                     degree_distribution = my.deg.dist,
                                                                     qual_modularity = linkrank_binary_modularity,
                                                                     complexity = my.complexity)
                                                    
                                                    # return to the list
                                                    cm
                                                  }# foreach
  # store null --------------------------------------------------------------
  cm.null.df <- bind_rows(community.null.metrics)
  
  write.csv2(cm.null.df,file = paste("results/community_metrics_null_",end.replicates,".csv",sep=""),row.names = F)
  
  stopCluster(cl)
  
  null.files <- list.files("results/",pattern = "community_metrics_null_",full.names = T)
  null.complete <- null.files %>% map_dfr(read_csv2)
  write.csv2(null.complete,file = "results/community_metrics_null.csv",row.names = F)
  
}# if include.null

