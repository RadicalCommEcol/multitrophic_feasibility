


# NOTE: we assume that there is no overlap between the members of herb and fv guilds.

#' mask for allowed interactions in a unipartite matrix
#' 
#' returns a matrix whose dimensions and (row, col) names are equal to those of my.matrix.
#' If an interaction (i, j) is forbidden, the element (i,j) is set to zero, otherwise it is
#' set to 1.
#'
#' @param my.matrix original interaction matrix
#' @param year.plot.matrix overall ommunity matrix, for getting sp names
#' @param guild.i which community are we looking at
#' @param sp.names names of all sp in the community
#' @param i.year year to filter interactions
#' @param forbidden_info_file a file with the structure of "results/potential_interactions.csv"
#'
#' @return matrix with the same dimensions as my.matrix
#' @export
allowed.interactions <- function(my.matrix,
                                 year.plot.matrix,
                                 guild.i,
                                 sp.names, 
                                 i.year,
                                 forbidden_info_file = "data/potential_interactions.csv"){
  
  allowed_interactions <- my.matrix
  allowed_interactions[,] <- 1
  
  if (guild.i %in% c("plants","floral visitors","herbivores")){
    
    return(allowed_interactions)
    
  }else{
    
    forbidden_plant_animal_interactions <- read.csv2(forbidden_info_file) %>%
      filter(year == years[i.year], structural.zero == T) %>% 
      dplyr::select(plant,animal,type)
    
    all.names <- rownames(year.plot.matrix)
    
    plant.names <- all.names[which(rownames(year.plot.matrix) %in% 
                                     sp.names[[i.year]][["plants"]])]
    herb.names <- all.names[which(rownames(year.plot.matrix) %in% 
                                    sp.names[[i.year]][["herbivores"]])]
    fv.names <- all.names[which(rownames(year.plot.matrix) %in% 
                                  sp.names[[i.year]][["floral.visitors"]])]
    
    # Remove forbidden 
    for (i.plant in 1:length(plant.names)){
      
      plant.index <- which(rownames(allowed_interactions) == plant.names[i.plant])
      forbidden.mates.names <-  forbidden_plant_animal_interactions %>%
        filter(plant == plant.names[i.plant]) %>% dplyr::select(animal) %>% unique() %>% pull()
      forbidden.mates.index <- which(rownames(allowed_interactions) %in% forbidden.mates.names)
      
      # Set forbidden interactions to zero
      allowed_interactions[plant.index, forbidden.mates.index] <- 0
      allowed_interactions[forbidden.mates.index, plant.index] <- 0
    }
    
    # Remove fv and herb. interactions when both guilds are present
    
    if(guild.i=="all"){
      
      for (i.fv in 1:length(fv.names)){
        
        fv.index <- which(rownames(allowed_interactions) == fv.names[i.fv])
        forbidden.mates.index <- which(rownames(allowed_interactions) %in% herb.names)
        
        # Set forbidden interactions to zero
        allowed_interactions[fv.index, forbidden.mates.index] <- 0
        allowed_interactions[forbidden.mates.index, fv.index] <- 0
      }
      
    }
    
    return(allowed_interactions)
    
  }

}
