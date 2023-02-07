
# this script is used to generate a spatial layout for plotting the 
# networks of Figure 1

# it is not possible to reproduce the figure perfectly because it involves
# manual positioning of the nodes, but it is -repeatable-

# note that it involves installing the add-in SNAhelper (see below)

# -------------------------------------------------------------------------

# generate a layout for all interactions in the multilayer net

library(tidygraph)
library(tidyverse)
library(graphlayouts)
library(ggraph)
library(igraph)

# build matrix simply from all species combinations
load(file = "results/community_names.RData")

all.names <- unique(unlist(sp.names))

# the interaction strenghts are irrelevant
my.matrix <- matrix(data = runif(length(all.names)^2,0,1),
                    nrow = length(all.names),
                    ncol = length(all.names),
                    dimnames = list(all.names,all.names))

# convert to tidygraph object, will be easier to plot

# 1 - nodes/edges dataframes
g  <- graph.adjacency(my.matrix,weighted=TRUE)
my.edge.list <- get.data.frame(g)

# 1.2 include type of nodes
id.list <- data.frame(ID = 1:length(unique(my.edge.list$from)),
                      node.name = sort(unique(my.edge.list$from)),
                      guild = NA_character_,stringsAsFactors = FALSE)

for(i.node in 1:nrow(id.list)){
  if(id.list$node.name[i.node] %in% sp.names[["2019"]][["plants"]]){
    id.list$guild[i.node] <- "plants"
  }else if(id.list$node.name[i.node] %in% sp.names[["2019"]][["floral.visitors"]]){
    id.list$guild[i.node] <- "floral.visitors"
  }else if(id.list$node.name[i.node] %in% sp.names[["2019"]][["herbivores"]]){
    id.list$guild[i.node] <- "herbivores"
  }
}# for i.node

# 1.3 edge list
my.edge.list$from.id <- id.list$ID[match(my.edge.list$from,id.list$node.name)]
my.edge.list$to.id <- id.list$ID[match(my.edge.list$to,id.list$node.name)]

# add edge type, just in case
for(i.edge in 1:nrow(my.edge.list)){
  if(my.edge.list$from[i.edge] %in% sp.names[["2019"]][["plants"]] & 
     my.edge.list$to[i.edge] %in% sp.names[["2019"]][["plants"]]){
    my.edge.list$edge.type[i.edge] <- "intraguild plants"
  }else if(my.edge.list$from[i.edge] %in% sp.names[["2019"]][["floral.visitors"]] & 
           my.edge.list$to[i.edge] %in% sp.names[["2019"]][["floral.visitors"]]){
    my.edge.list$edge.type[i.edge] <- "intraguild floral visitors"
  }else if(my.edge.list$from[i.edge] %in% sp.names[["2019"]][["herbivores"]] & 
           my.edge.list$to[i.edge] %in% sp.names[["2019"]][["herbivores"]]){
    my.edge.list$edge.type[i.edge] <- "intraguild herbivores"
  }else if((my.edge.list$from[i.edge] %in% sp.names[["2019"]][["plants"]] & 
           my.edge.list$to[i.edge] %in% sp.names[["2019"]][["floral.visitors"]]) |
           (my.edge.list$from[i.edge] %in% sp.names[["2019"]][["floral.visitors"]] & 
            my.edge.list$to[i.edge] %in% sp.names[["2019"]][["plants"]])){
    my.edge.list$edge.type[i.edge] <- "plant-floral visitors"
  }else if((my.edge.list$from[i.edge] %in% sp.names[["2019"]][["plants"]] & 
            my.edge.list$to[i.edge] %in% sp.names[["2019"]][["herbivores"]]) |
           (my.edge.list$from[i.edge] %in% sp.names[["2019"]][["herbivores"]] & 
            my.edge.list$to[i.edge] %in% sp.names[["2019"]][["plants"]])){
    my.edge.list$edge.type[i.edge] <- "plant-herbivores"
  # just in case
  }else if((my.edge.list$from[i.edge] %in% sp.names[["2019"]][["herbivores"]] & 
            my.edge.list$to[i.edge] %in% sp.names[["2019"]][["floral.visitors"]]) |
           (my.edge.list$from[i.edge] %in% sp.names[["2019"]][["floral.visitors"]] & 
            my.edge.list$to[i.edge] %in% sp.names[["2019"]][["herbivores"]])){
    my.edge.list$edge.type[i.edge] <- "herbivores-floral visitors"
  }
}# for i.node

my.edge.list <- my.edge.list[,c("from.id","to.id","weight","edge.type")]
names(my.edge.list)[1:2] <- c("from","to")
my.edge.list <- subset(my.edge.list,edge.type != "herbivores-floral visitors")
my.net <- tbl_graph(nodes = id.list,edges = my.edge.list, directed = TRUE)


# -------------------------------------------------------------------------
# NOTE: 
# the layout involves visually placing nodes in the 2-d canvas
# at this stage i opened the SNA add-in, and manually tweaked the layout
# https://github.com/schochastics/snahelper
# this will generate X and Y coordinates for each node of the network object
# -------------------------------------------------------------------------

# when happy, save the layout
# NOTE: it will not work without the previous step, as by default the network
# object does not contain coordinates
layout.metaweb <- data.frame(ID = as.numeric(V(my.net)),x = V(my.net)$x,y = V(my.net)$y)
layout.metaweb$species <- id.list$node.name[match(layout.metaweb$ID,id.list$ID)]

write.csv2(x = layout.metaweb[,c("ID","species","x","y")],file = "results/network_plot_layout.csv",row.names = FALSE)

# template for plotting
# library(colorblindr)
# ggraph(my.net, layout = "manual", x = layout.metaweb$x, y = layout.metaweb$y) + 
#    geom_edge_link0(edge_colour = "#A8A8A8", 
# 	                 edge_width = 0.1, 
# 	                 edge_alpha = 1) + 
# 	 geom_node_point(aes(fill = guild), 
# 	                 colour = "#000000", 
# 	                 size = 8, 
# 	                 shape = 21, 
# 	                 stroke = 0.3) + 
#     scale_fill_OkabeIto() + 
# 	 theme_graph() + 
# 	 # theme(legend.position = "none") +
#   NULL





