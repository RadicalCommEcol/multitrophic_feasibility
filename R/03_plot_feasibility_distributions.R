
library(tidyverse)
library(cowplot)
library(colorblindr)
library(patchwork)

# read data ---------------------------------------------------------------
# vers <- "_v2"
vers <- ""

file.name <- "results/feasibility_domain_observed"
null.name <- "results/feasibility_domain_null"

file.name <- paste(file.name,vers,".csv",sep="")
null.name <- paste(null.name,vers,".csv",sep="")

fd <- read.csv2(file = file.name,
                header = TRUE,
                stringsAsFactors = FALSE)
names(fd)[names(fd) == "guild"] <- "guilds"

fd.null <- read.csv2(file = null.name,
                     header = TRUE,
                     stringsAsFactors = FALSE)

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

# -------------------------------------------------------------------------
# add null observations
fd.null.clean <- fd.null %>%
  select(omega_mean,year,plot,guild,
         intraguild.type,null.replicate,null.type)

fd.null.clean$fd.average <- fd.null.clean$omega_mean
fd.null.clean$guilds <- factor(fd.null.clean$guild,levels = c("plants",
                                                              "floral visitors",
                                                              "herbivores",
                                                              "plants-floral visitors",
                                                              "plants-herbivores","all"))
fd.null.clean$guilds <- recode(fd.null.clean$guilds, "plants" = "Plants",
                               "floral visitors" = "Pollinators",
                               "herbivores" = "Herbivores",
                               "plants-floral visitors" = "Plants-Pollinators",
                               "plants-herbivores" = "Plants-Herbivores",
                               "all" = "All")

fd.null.clean$plot <- as.factor(fd.null.clean$plot)

fd.null.clean$intraguild.type[fd.null.clean$intraguild.type == "mean_field"] <- "Mean field"
fd.null.clean$intraguild.type[fd.null.clean$intraguild.type == "nesting_larvae"] <- "Resources"
fd.null.clean$intraguild.type[fd.null.clean$intraguild.type == "nesting_larvae_phenology"] <- "Resources and \nphenology"

fd.null.clean$type <- "randomized"
fd$type <- "observed"
fd$null.replicate <- NA

fd.null.clean.2 <- fd.null.clean %>% select(fd.average,intraguild.type,guilds,type)

# -------------------------------------------------------------------------

fd.obs.all <- fd.clean %>% filter(guilds == "All")
fd.null.all <- fd.null.clean.2 %>% filter(guilds == "All")

# obtain moments and 97.5 percentile of the null distributions
fd.null.dist <- fd.null.all %>% 
  group_by(intraguild.type, type) %>%
  summarise(null.mean = mean(fd.average),
            null.sd = sd(fd.average),
            null.95.perc = qnorm(.975,mean(fd.average),sd(fd.average)))
fd.obs.dist <- fd.obs.all %>% left_join(fd.null.dist) %>%
  group_by(intraguild.type) %>%
  mutate(null.dist.perc = pnorm(fd.average,null.mean,null.sd))

# -------------------------------------------------------------------------
# feasibility domain

pfd <- ggplot(fd.null.all, aes(x = fd.average)) + 
  geom_density(aes(color = intraguild.type, fill = intraguild.type), alpha = .4) + 
  geom_vline(data = fd.null.dist,aes(xintercept = null.95.perc), color = "darkred") +
  geom_point(data = fd.obs.all, aes(x = fd.average, 
                                    y = 0, 
                                    color = intraguild.type),
             size = 2) +
  facet_grid(intraguild.type~.) +
  scale_color_OkabeIto(darken = 0.2) +
  scale_fill_OkabeIto(darken = 0.2) +
  theme_bw() +
  theme(strip.background = element_blank())+
  # theme(strip.text.y = element_blank()) +
  theme(legend.position="none")+
  labs(x = "feasibility domain",y = "density estimate") +
  NULL
# pfd

# ggsave(filename = paste("results/images/Fig_null_fd.pdf",sep=""),
#        plot = pfd,
#        device = cairo_pdf,
#        width = 6,height = 4,dpi = 300)
