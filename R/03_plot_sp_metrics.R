
# plot probabilities of exclusion

# INPUT
# probabilities of exclusion: "results/exclusion_probabilities_observed.csv" 
# community names: "results/community_names.RData"

# OUTPUT
# Figure 3 of the paper

# -------------------------------------------------------------------------
library(tidyverse)
library(patchwork)
library(colorblindr)

# -------------------------------------------------------------------------

sp.obs <- read.csv2("results/exclusion_probabilities_observed.csv")

# -------------------------------------------------------------------------

# obtain the final metric, 
# relative to a situation with equal probabilities of exclusion

sp.obs <- sp.obs %>% 
  group_by(year,plot,guild,intraguild.type) %>%
  mutate(prob_excl_rel = (prob_excl_mean/(1/n())))

# -------------------------------------------------------------------------
# subset to be consistent with the rest of analyses

sp.clean <- subset(sp.obs, guild == "all")

sp.clean$intraguild.type[sp.clean$intraguild.type == "mean_field"] <- "Mean field"
sp.clean$intraguild.type[sp.clean$intraguild.type == "nesting_larvae"] <- "Resources"
sp.clean$intraguild.type[sp.clean$intraguild.type == "nesting_larvae_phenology"] <- "Resources and \nphenology"

# add species guild
# load sp names as well
load(file = paste("results/community_names.RData",sep=""))
plants <- sort(unique(c(sp.names[["2019"]][["plants"]],sp.names[["2020"]][["plants"]])))
fv <- sort(unique(c(sp.names[["2019"]][["floral.visitors"]],sp.names[["2020"]][["floral.visitors"]])))
herb <- sort(unique(c(sp.names[["2019"]][["herbivores"]],sp.names[["2020"]][["herbivores"]])))

sp.clean$sp.guild <- "plants"
sp.clean$sp.guild[sp.clean$species %in% fv] <- "pollinators"
sp.clean$sp.guild[sp.clean$species %in% herb] <- "herbivores"
sp.clean$sp.guild <- factor(sp.clean$sp.guild, levels = c("plants","pollinators",
                                                          "herbivores"))
# -------------------------------------------------------------------------

exc.means <- sp.clean %>% 
  group_by(intraguild.type, sp.guild) %>%
  summarise(exc.avg = mean(log(prob_excl_rel),na.rm = TRUE),
            exc.sd = sd(log(prob_excl_rel),na.rm = TRUE),
            exc.se = sd(log(prob_excl_rel),na.rm = TRUE)/sqrt(n()))

# plot

exc.plot <- ggplot(sp.clean, aes(x = sp.guild, y = log(prob_excl_rel))) +
  geom_point(aes(color = intraguild.type),
             position=position_jitterdodge(jitter.width = .11),
             size = .95, alpha = .4) +
  geom_boxplot(aes(fill = intraguild.type),
               outlier.shape = NA,
               alpha = .4) +
  scale_color_OkabeIto(darken = .2,name = "Intra-guild \ninteractions") +
  scale_fill_OkabeIto(darken = .2,name = "Intra-guild \ninteractions") +
  labs(x = "", y = "log(exclusion ratios)") +
  theme_bw() +
  theme(strip.background = element_blank())+
  theme(axis.ticks.x = element_blank()) +
  theme(panel.grid.minor.x = element_blank(),
        panel.grid.major.x = element_blank()) +
  NULL
# exc.plot

# ggsave(filename = paste("results/images/Fig_prob_excl.pdf",sep=""),
#        plot = exc.plot,
#        device = cairo_pdf,
#        width = 7,height = 4,dpi = 300)






