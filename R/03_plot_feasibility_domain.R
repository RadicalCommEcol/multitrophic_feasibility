
# plot feasibility results

# INPUT
# feasibility domains: "results/feasibility_domain_observed" 

# OUTPUT
# Figure 2 of the paper

# -------------------------------------------------------------------------
library(tidyverse)
library(colorblindr)
library(patchwork)

# read data ---------------------------------------------------------------
# vers <- "_v2"
vers <- ""

file.name <- "results/feasibility_domain_observed"

file.name <- paste(file.name,vers,".csv",sep="")

fd <- read.csv2(file = file.name,
                         header = TRUE,
                         stringsAsFactors = FALSE)
names(fd)[names(fd) == "guild"] <- "guilds"

# clean data ---------------------------------------------------------------

fd$fd.average <- fd$omega_mean
fd$guilds <- factor(fd$guilds,levels = c("plants",
                                         "floral visitors",
                                         "herbivores",
                                         "plants-floral visitors",
                                         "plants-herbivores","all")) 
fd$guilds <- recode(fd$guilds, "plants" = "Plants", 
                    "floral visitors" = "Pollinators",
                    "herbivores" = "Herbivores",
                    "plants-floral visitors" = "Plants-Pollinators",
                    "plants-herbivores" = "Plants-Herbivores",
                    "all" = "All")

fd$plot <- as.factor(fd$plot)

fd$intraguild.type[fd$intraguild.type == "mean_field"] <- "Mean field"
fd$intraguild.type[fd$intraguild.type == "nesting_larvae"] <- "Resources"
fd$intraguild.type[fd$intraguild.type == "nesting_larvae_phenology"] <- "Resources and \nphenology"

fd.clean <- fd %>% select(fd.average,intraguild.type,guilds)

fd.all <- fd.clean

# plot --------------------------------------------------------------------

pj <- .20

standard.fd <- ggplot(fd.all, aes(x = intraguild.type, y = fd.average)) + 
  geom_point(aes(color = intraguild.type),
             # shape = 21, 
             position=position_jitter(width = pj),
             # size = .85, 
             alpha = .7) +
  geom_boxplot(aes(fill = intraguild.type),
               outlier.shape = NA,
               alpha = .4) +
  scale_fill_OkabeIto(darken = 0.2, name = "Intra-guild \ncompetition") +
  scale_color_OkabeIto(darken = 0.2, name = "Intra-guild \ncompetition") +
  # scale_shape_manual(values = c(21,23)) +
  facet_grid(.~guilds, drop = T)+
  labs(x = "",y = "Feasibility domain") +
  # xlim(0,0.26) +
  theme_bw()+
  # theme(axis.text.x  = element_text(angle = 60, hjust = 1))+
  theme(strip.background = element_blank())+
  # scale_y_discrete(breaks=NULL, limits=rev) +
  # scale_x_continuous(breaks=seq(0,0.5,by = 0.05), limits = c(0,0.5)) +
  # guides(color=guide_legend(override.aes = list(shape=21))) + 
  theme(axis.ticks.x = element_blank(), axis.text.x = element_blank()) +
  theme(panel.grid.minor.x = element_blank(),
        panel.grid.major.x = element_blank()) +
  NULL
# standard.fd

# ggsave(filename = paste("results/images/standard_fd",vers,".pdf",sep=""),plot = standard.fd,
#        device = cairo_pdf,
#        width = 10,height = 4,dpi = 300)

