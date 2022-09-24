
infomap_modularity <- function(A,
                allowed.interactions.mat,
                damping, 
                pr.algo,
                temp.dir = NULL){
  g  <- graph.adjacency(A,weighted=TRUE)
  
  # To run infomap it seems that weights should be non-negative
  my.edge.list <- get.data.frame(g) %>% mutate(weight=abs(weight))
  
  nodes <- my.edge.list$from %>% unique()
  
  nodes.ID <- tibble(node_id=as.numeric(1:length(nodes)),species=nodes)
  
  # Preparing edge.list and nodes.ID to run infomap
  
  my.edge.list.ID <- my.edge.list %>% rename(species=from) %>%
    left_join(nodes.ID,by="species") %>% select(-species) %>%
    rename(from=node_id,species=to) %>%
    left_join(nodes.ID,by="species") %>%
    select(-species) %>% rename(to=node_id) %>% select(from,to,weight)
  
  nodes.ID2 <- nodes.ID %>% rename(node_name=species)
  
  infomap_mono <- create_monolayer_object(x = my.edge.list,
                                          directed = T,
                                          bipartite = F,
                                          node_metadata = nodes.ID2[,c(2,1)])
  
  # infomap_mono2 <- infomap_mono
  # infomap_mono2$nodes$node_name <- as.numeric(infomap_mono2$nodes$node_name)
  
  # run Infomap
  if(!is.null(temp.dir)){
    # temp_dir <- paste("R/temp/d",i.id,sep="")
    # temp_dir <- "R/temp/d1"
    dir.create(temp.dir)
    modules_relax_rate <- run_infomap_monolayer3(x = infomap_mono,
                                                 flow_model = 'directed',
                                                 # infomap_executable = paste(infomap_dir,"Infomap",sep=""),
                                                 infomap_executable = "Infomap",
                                                 temp_dir = temp.dir,
                                                 silent=T,
                                                 trials=1000,
                                                 two_level=T,
                                                 seed=200952)
    unlink(temp.dir, recursive = TRUE)
  }else{
    modules_relax_rate <- run_infomap_monolayer2(x = infomap_mono,
                                                 flow_model = 'directed',
                                                 # infomap_executable = paste(infomap_dir,"Infomap",sep=""),
                                                 infomap_executable = "Infomap",
                                                 # temp_dir = temp_dir,
                                                 silent=T,
                                                 trials=1000,
                                                 two_level=T,
                                                 seed=200952)
  }

  
  my.modularity <- modules_relax_rate[["L"]] # modularity in bits
  
  # Extract module information
  modules <- modules_relax_rate$modules %>%
    dplyr::select(node_id,module_level1) %>%
    rename(module=module_level1) %>%
    left_join(nodes.ID,by="node_id")
  
  
  # modules.aux <- tibble(year = my.year,
  #                       plot = my.plot,
  #                       guild = my.guild,
  #                       species = modules$species,
  #                       module = modules$module)
  # 
  # module_members <- bind_rows(module_members,modules.aux)
  
  # If we dont want to use bits, I guess that we could translate the previous partition
  # to other units by uning the alternative definition of modularity
  # E. A. Leicht and M. E. J. Newman, Phys. Rev. Lett. 100, 118703 (2008).
  
  # Note 2: infomap do not optimise this generalized modularity function. Thus,
  # I dont expect high values.
  
  # According to results for plot 1, year 2019 and guild == "floral visitors",
  # it seems that it seems that isolated nodes module is NA (see species "Apoidea")
  # For each one of those nodes, we will add a new module. We denote such partition
  # as "corrected partition"
  
  # CORRECTED partition from infomap
  module_max <- max(modules$module,na.rm = T)
  
  for(i in 1:nrow(modules)){
    
    if(is.na(modules$module[i])){
      module_max <- module_max + 1
      modules$module[i] <- module_max
    }
  }
  
  # linkrank modularity -----------------------------------------------------
  # Note 1: I will use non-negative weights to be consistent with
  # the inputs that we used to feed the infomap algorithm.
  
  g  <- graph.adjacency(abs(A),weighted=TRUE)
  
  # this returns a "mask", a matrix of the same dimensions as my.matrix
  # with zeros if an interaction is forbidden, 1 if it is allowed.
  # allowed.interactions.mat <- allowed.interactions(my.matrix.null, year.plot.matrix,
  #                                                  my.guild,
  #                                                  sp.names, my.year)
  
  linkrank_modularity <- linkRankModularity(g,
                                            partition=modules$module,
                                            allowed.interactions.mat = allowed.interactions.mat,
                                            damping = damping,
                                            pr.algo = pr.algo)
  return(linkrank_modularity)
}

