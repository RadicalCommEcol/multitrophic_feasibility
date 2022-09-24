
# obtain species-level metrics from the observed communities. There is currently
# no code implemented for the randomized communities.

# INPUT
# community matrices: "results/community_matrices.RData"
# community names: "results/community_names.RData"

# OUTPUT
# dataframe of species-level metrics: "results/species_level_metrics.csv"

# -------------------------------------------------------------------------

library(tidyverse)

load(file = paste("results/community_matrices.RData",sep=""))
source("R/interaction_overlap.R")

# load sp names as well
load(file = paste("results/community_names.RData",sep=""))
plants <- sort(unique(c(sp.names[["2019"]][["plants"]],sp.names[["2020"]][["plants"]])))
fv <- sort(unique(c(sp.names[["2019"]][["floral.visitors"]],sp.names[["2020"]][["floral.visitors"]])))
herb <- sort(unique(c(sp.names[["2019"]][["herbivores"]],sp.names[["2020"]][["herbivores"]])))

years <- c(2019,2020)
plots <- 1:9

guild.combinations <- c("plants",
                        "floral visitors",
                        "herbivores",
                        "plants-floral visitors",
                        "plants-herbivores",
                        "all")

intraguild.types <- names(community_matrices)

# -------------------------------------------------------------------------
sp.metrics.list <- list()

# test
# i.type <- i.year <- i.plot <- 2
# i.guild <- 4

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
        
        present.sp <- rownames(my.matrix)
        
        df.sp <- expand.grid(year = years[i.year],
                             plot = i.plot,
                             sp.guild = NA,
                             guilds = guild.combinations[i.guild],
                             intraguild.type = intraguild.types[i.type],
                             species = present.sp,
                             diagonal_dominance = NA,
                             avg_intraguild_niche_overlap = NA,
                             avg_interguild_niche_overlap = NA,
                             total_intraguild_niche_overlap = NA,
                             total_interguild_niche_overlap = NA,
                             in_degree = NA,
                             out_degree = NA,stringsAsFactors = FALSE)
        
        # diagonal dominance does not make sense for bipartite matrices
        if(guild.combinations[i.guild] %in% c("plants","floral visitors","herbivores","all")){
          # diagonal dominance
          # diagonal elements
          d <- diag(my.matrix)
          # rowsums
          nd <- rowSums(my.matrix)
          # but without the diagonal
          nd <- nd - d
          
          df.sp$diagonal_dominance <- abs(d) - abs(nd)
        }
        
        
        # -------------------------------------------------------------------------
        # niche overlap for all species
        # differentiate intra and inter-guild overlap
        
        if(guild.combinations[i.guild] %in% c("plants","floral visitors","herbivores")){
          
          # in these cases, only intra-guild overlap
          pair.intra.overlap <- interaction_overlap(abs(my.matrix))
          ov.clean <- subset(pair.intra.overlap, sp1 != sp2 & !is.nan(overlap) & !is.na(overlap))
          # from https://stackoverflow.com/questions/32669609/identifying-unique-pairs-of-values-from-two-columns-in-a-dataframe
          ov.unique <-  ov.clean[!duplicated(t(apply(ov.clean, 1, sort))),]
          
          my.intra.overlap <- ov.unique
          my.inter.overlap <- NA
          
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
          
          my.intra.overlap <- bind_rows(intra.plants.unique,
                                         intra.fv.unique)
          my.inter.overlap <- bind_rows(inter.plants.fv.unique,
                                         inter.fv.plants.unique)
          
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
          
          my.intra.overlap <- bind_rows(intra.plants.unique,
                                        intra.herb.unique)
          my.inter.overlap <- bind_rows(inter.plants.herb.unique,
                                        inter.herb.plants.unique)
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
          
          my.intra.overlap <- bind_rows(intra.plants.unique,
                                        intra.fv.unique,
                                        intra.herb.unique)
          my.inter.overlap <- bind_rows(inter.plants.herb.unique,
                                        inter.herb.plants.unique,
                                        inter.plants.fv.unique,
                                        inter.fv.plants.unique)
        }
        
        # -------------------------------------------------------------------------
        # to harmonize uni and bipartite nets, better to go one sp at a time
        for(i.sp in 1:nrow(df.sp)){
          
          # first, assign guild
          my.sp <- df.sp$species[i.sp]
          
          if(my.sp %in% plants){
            df.sp$sp.guild[i.sp] <- "plants"
          }else if(my.sp %in% fv){
            df.sp$sp.guild[i.sp] <- "floral visitors"
          }else if(my.sp %in% herb){
            df.sp$sp.guild[i.sp] <- "herbivores"
          }
          
          # exclude the diagonal
          df.sp$in_degree[i.sp] <- sum(my.matrix[my.sp,] != 0) - 1
          df.sp$out_degree[i.sp] <- sum(my.matrix[,my.sp] != 0) - 1
          
          # niche overlap -----------------------------------------------------------
          intra.my.sp <- my.intra.overlap %>% filter(sp1 == my.sp | sp2 == my.sp)
          df.sp$avg_intraguild_niche_overlap[i.sp] <- mean(intra.my.sp$overlap,na.rm = TRUE)
          df.sp$total_intraguild_niche_overlap[i.sp] <- sum(intra.my.sp$overlap,na.rm = TRUE)
          
          if(is.data.frame(my.inter.overlap)){
            inter.my.sp <- my.inter.overlap %>% filter(sp1 == my.sp | sp2 == my.sp)
            df.sp$avg_interguild_niche_overlap[i.sp] <- mean(inter.my.sp$overlap,na.rm = TRUE)  
            df.sp$total_interguild_niche_overlap[i.sp] <- sum(inter.my.sp$overlap,na.rm = TRUE)  
          }else{
            inter.my.sp <- NA
          }
          
        }# for each sp
        
        sp.metrics.list[[length(sp.metrics.list)+1]] <- df.sp
        
      }# i.guild
    }# i.plot
  }# i.year
}# i.type

sp.metrics <- bind_rows(sp.metrics.list)

# -------------------------------------------------------------------------
write.csv2(sp.metrics,"results/species_level_metrics.csv",row.names = F)


