
# Calculate rarefaction curves of species interactions 
# for each type of interaction observed (plant-plant, plant-herb, plant-pol)

# INPUT
# interaction matrices (/data/*_matrices.RData)

# OUTPUT
# rarefaction curves for every interaction and local community

# -------------------------------------------------------------------------

library(tidyverse)
library(vegan)
library(iNEXT)
library(patchwork)

# -------------------------------------------------------------------------
# add a version suffix?
vers.out <- ""
vers <- ""

# -------------------------------------------------------------------------

years <- c(2019,2020)
plots <- 1:9

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
# rarefaction/extrapolation

# first, convert matrices to 1d vectors, usable by iNEXT, and store them in lists
pp.interaction.list <- pp.all.years
ph.interaction.list <- ph.all.years
pfv.interaction.list <- pfv.all.years

for(i.year in 1:length(years)){
  for(i.plot in plots){
    
    # plant-plant interactions
    my.pp.vec <- c(t(pp.all.years[[i.year]][[i.plot]]))
    my.pp.names.grid <- expand.grid(colnames(pp.all.years[[i.year]][[i.plot]]), 
                                    rownames(pp.all.years[[i.year]][[i.plot]]))
    names(my.pp.vec) <- sprintf('%s_%s', my.pp.names.grid[,1], my.pp.names.grid[,2])
    
    pp.interaction.list[[i.year]][[i.plot]] <- my.pp.vec
    
    # plant-herb interactions
    my.ph.vec <- c(t(ph.all.years[[i.year]][[i.plot]]))
    my.ph.names.grid <- expand.grid(colnames(ph.all.years[[i.year]][[i.plot]]), 
                                    rownames(ph.all.years[[i.year]][[i.plot]]))
    names(my.ph.vec) <- sprintf('%s_%s', my.ph.names.grid[,1], my.ph.names.grid[,2])
    
    ph.interaction.list[[i.year]][[i.plot]] <- my.ph.vec
    
    # plant-pol interactions
    my.pfv.vec <- c(t(pfv.all.years[[i.year]][[i.plot]]))
    my.pfv.names.grid <- expand.grid(colnames(pfv.all.years[[i.year]][[i.plot]]), 
                                    rownames(pfv.all.years[[i.year]][[i.plot]]))
    names(my.pfv.vec) <- sprintf('%s_%s', my.pfv.names.grid[,1], my.pfv.names.grid[,2])
    
    pfv.interaction.list[[i.year]][[i.plot]] <- my.pfv.vec
    
  }# for i.plot
}# for i.year

# second, generate iNEXT curves

# these lists will hold the iNEXT plots
pp.plot.list <- list()
ph.plot.list <- list()
pfv.plot.list <- list()

community.sc <- expand.grid(year = years,plot = plots, 
                            pp = NA, pfv = NA, ph = NA)

for(i.year in 1:length(years)){
  for(i.plot in plots){
    
    # plant-plant interactions
    pp.endpoint <- round(sum(pp.interaction.list[[i.year]][[i.plot]]) + .25 * sum(pp.interaction.list[[i.year]][[i.plot]]))
    pp.out <- iNEXT(unlist(pp.interaction.list[[i.year]][[i.plot]]), q=c(0), datatype="abundance", endpoint=pp.endpoint)
    pp.plot <- ggiNEXT(pp.out, type=2,color.var = "Order.q")+
      scale_shape_manual(values = c(19))+
      theme_bw(base_size = 12) +
      theme(legend.position="none")+
      labs(x="",y="") +
      ggtitle(paste("plant-plant: ",years[i.year]," plot ",i.plot,sep=""))
    pp.plot.list[[length(pp.plot.list)+1]] <- pp.plot

    # plant-herb interactions
    ph.endpoint <- round(sum(ph.interaction.list[[i.year]][[i.plot]]) + .25 * sum(ph.interaction.list[[i.year]][[i.plot]]))
    ph.out <- iNEXT(unlist(ph.interaction.list[[i.year]][[i.plot]]), q=c(0), datatype="abundance", endpoint=ph.endpoint)
    ph.plot <- ggiNEXT(ph.out, type=2,color.var = "Order.q")+
      scale_shape_manual(values = c(19))+
      theme_bw(base_size = 12) +
      theme(legend.position="none")+
      labs(x="",y="") +
      ggtitle(paste("plant-herbivore: ",years[i.year]," plot ",i.plot,sep=""))
    ph.plot.list[[length(ph.plot.list)+1]] <- ph.plot   
    
    # plant-pol interactions
    pfv.endpoint <- round(sum(pfv.interaction.list[[i.year]][[i.plot]]) + .25 * sum(pfv.interaction.list[[i.year]][[i.plot]]))
    pfv.out <- iNEXT(unlist(pfv.interaction.list[[i.year]][[i.plot]]), q=c(0), datatype="abundance", endpoint=pfv.endpoint)
    pfv.plot <- ggiNEXT(pfv.out, type=2,color.var = "Order.q")+
      scale_shape_manual(values = c(19))+
      theme_bw(base_size = 12) +
      theme(legend.position="none")+
      labs(x="",y="") +
      ggtitle(paste("plant-pollinator: ",years[i.year]," plot ",i.plot,sep=""))
    pfv.plot.list[[length(pfv.plot.list)+1]] <- pfv.plot  
    
    pos <- which(community.sc$year == years[i.year] & community.sc$plot == i.plot)
    community.sc$pp[pos] <- pp.out$DataInfo$SC
    community.sc$ph[pos] <- ph.out$DataInfo$SC
    community.sc$pfv[pos] <- pfv.out$DataInfo$SC
    
  }# for i.plot
}# for i.year

mean(community.sc$pp)
mean(community.sc$pfv)
mean(community.sc$ph)

sd(community.sc$pp)
sd(community.sc$pfv)
sd(community.sc$ph)

# all interactions together
mean(as.matrix(community.sc[,3:5]))

# -------------------------------------------------------------------------

pp.plot.full <- wrap_plots(pp.plot.list)
ph.plot.full <- wrap_plots(ph.plot.list)
pfv.plot.full <- wrap_plots(pfv.plot.list)

# -------------------------------------------------------------------------

# ggsave("results/images/Fig_plant_plant_rarefaction.pdf",pp.plot.full,
#        device = cairo_pdf,
#        width = 16,height = 12,dpi = 300)
# 
# ggsave("results/images/Fig_plant_herb_rarefaction.pdf",ph.plot.full,
#        device = cairo_pdf,
#        width = 16,height = 12,dpi = 300)
# 
# ggsave("results/images/Fig_plant_pol_rarefaction.pdf",pfv.plot.full,
#        device = cairo_pdf,
#        width = 16,height = 12,dpi = 300)
