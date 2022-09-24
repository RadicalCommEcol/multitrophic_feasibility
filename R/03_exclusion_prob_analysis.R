
# statistical analyses at the species level

# INPUT
# probabilities of exclusion: "results/exclusion_probabilities_observed.csv"
# "results/exclusion_probabilities_null.csv"
# species-level metrics: "results/species_level_metrics.csv"
# community-level metrics: "results/community_metrics.csv"

# OUTPUT
# Tables of the paper: 
#   Probabilities of exclusion - number of guilds: Table S3
#   Probabilities of exclusion - richness: Table S4
#   Null probabilities of exclusion: Table S5
#   Probabilities of exclusion - species metrics: Table S6

# Figures
#   Probabilities of exclusion - number of guilds: Figure S2
#   Probabilities of exclusion - richness: Figure S3
#   Null probabilities of exclusion: Figure S4

# -------------------------------------------------------------------------

library(tidyverse)
library(patchwork)
library(DHARMa)
library(broom.mixed)
library(colorblindr)
# source("R/feasibility_functions/anisotropy_index.R")

# -------------------------------------------------------------------------

exc.obs <- read.csv2("results/exclusion_probabilities_observed.csv")
exc.null <- read.csv2("results/exclusion_probabilities_null.csv")
names(exc.obs)[names(exc.obs) == "guild"] <- "guilds"
sp.metrics <- read.csv2("results/species_level_metrics.csv")
comm.metrics <- read.csv2("results/community_metrics.csv")

# -------------------------------------------------------------------------

# obtain the final metric, 
# relative to a situation with equal probabilities of exclusion

exc.obs <- exc.obs %>% 
  group_by(year,plot,guilds,intraguild.type) %>%
  mutate(prob_excl_rel = (prob_excl_mean/(1/n())))

exc.null <- exc.null %>% 
  group_by(year,plot,guilds,intraguild.type) %>%
  mutate(prob_excl_rel = (prob_excl_mean/(1/n())))

# -------------------------------------------------------------------------
# clean up names

exc.obs$guilds <- factor(exc.obs$guilds,levels = c("plants",
                                                 "floral visitors",
                                                 "herbivores",
                                                 "plants-floral visitors",
                                                 "plants-herbivores","all")) 
exc.obs$guilds <- recode(exc.obs$guilds, "plants" = "Plants", 
                        "floral visitors" = "Pollinators",
                        "herbivores" = "Herbivores",
                        "plants-floral visitors" = "Plants-Pollinators",
                        "plants-herbivores" = "Plants-Herbivores",
                        "all" = "All")

exc.obs$intraguild.type[exc.obs$intraguild.type == "mean_field"] <- "Mean field"
exc.obs$intraguild.type[exc.obs$intraguild.type == "nesting_larvae"] <- "Resources"
exc.obs$intraguild.type[exc.obs$intraguild.type == "nesting_larvae_phenology"] <- "Resources and \nphenology"

exc.null$guilds <- factor(exc.null$guilds,levels = c("plants",
                                                   "floral visitors",
                                                   "herbivores",
                                                   "plants-floral visitors",
                                                   "plants-herbivores","all")) 
exc.null$guilds <- recode(exc.null$guilds, "plants" = "Plants", 
                         "floral visitors" = "Pollinators",
                         "herbivores" = "Herbivores",
                         "plants-floral visitors" = "Plants-Pollinators",
                         "plants-herbivores" = "Plants-Herbivores",
                         "all" = "All")

exc.null$intraguild.type[exc.null$intraguild.type == "mean_field"] <- "Mean field"
exc.null$intraguild.type[exc.null$intraguild.type == "nesting_larvae"] <- "Resources"
exc.null$intraguild.type[exc.null$intraguild.type == "nesting_larvae_phenology"] <- "Resources and \nphenology"

comm.metrics$guilds <- factor(comm.metrics$guilds,levels = c("plants",
                                                           "floral visitors",
                                                           "herbivores",
                                                           "plants-floral visitors",
                                                           "plants-herbivores","all")) 
comm.metrics$guilds <- recode(comm.metrics$guilds, "plants" = "Plants", 
                             "floral visitors" = "Pollinators",
                             "herbivores" = "Herbivores",
                             "plants-floral visitors" = "Plants-Pollinators",
                             "plants-herbivores" = "Plants-Herbivores",
                             "all" = "All")

comm.metrics$intraguild.type[comm.metrics$intraguild.type == "mean_field"] <- "Mean field"
comm.metrics$intraguild.type[comm.metrics$intraguild.type == "nesting_larvae"] <- "Resources"
comm.metrics$intraguild.type[comm.metrics$intraguild.type == "nesting_larvae_phenology"] <- "Resources and \nphenology"

sp.metrics$guilds <- factor(sp.metrics$guilds,levels = c("plants",
                                                             "floral visitors",
                                                             "herbivores",
                                                             "plants-floral visitors",
                                                             "plants-herbivores","all")) 
sp.metrics$guilds <- recode(sp.metrics$guilds, "plants" = "Plants", 
                              "floral visitors" = "Pollinators",
                              "herbivores" = "Herbivores",
                              "plants-floral visitors" = "Plants-Pollinators",
                              "plants-herbivores" = "Plants-Herbivores",
                              "all" = "All")

sp.metrics$intraguild.type[sp.metrics$intraguild.type == "mean_field"] <- "Mean field"
sp.metrics$intraguild.type[sp.metrics$intraguild.type == "nesting_larvae"] <- "Resources"
sp.metrics$intraguild.type[sp.metrics$intraguild.type == "nesting_larvae_phenology"] <- "Resources and \nphenology"

# -------------------------------------------------------------------------
# 1 - probabilities of exclusion and complexity

exc.obs$guild.num <- "one"
exc.obs$guild.num[exc.obs$guilds %in% c("Plants-Pollinators",
                                      "Plants-Herbivores")] <- "two"
exc.obs$guild.num[exc.obs$guilds == "All"] <- "three"
exc.obs$guild.num <- factor(exc.obs$guild.num,levels = c("one","two","three"))
exc.obs$sp.guild <- sp.metrics$sp.guild[match(exc.obs$species,sp.metrics$species)]
exc.obs$sp.guild <- factor(exc.obs$sp.guild, levels = c("plants","herbivores","floral visitors"))  

richness.data <- comm.metrics %>% 
  filter(metric == "richness") %>% 
  pivot_wider(names_from = "metric",values_from = "value")

r1.data <- left_join(exc.obs,richness.data) %>%
  filter(!is.na(sp.guild))
r1.data$richness.scaled <- scales::rescale(r1.data$richness)

# -------------------------------------------------------------------------
exc.dif <- lm(log(prob_excl_rel) ~ guild.num * intraguild.type * sp.guild, data = r1.data)
exc.aov3 <- car::Anova(exc.dif,type = 3)

# print(xtable::xtable(as.data.frame(broom::tidy(exc.aov3)),
#                      floating=FALSE,
#                      digits = 3,
#                      latex.environments=NULL,
#                      booktabs=FALSE))

# -------------------------------------------------------------------------
exc.dif.p1 <- ggplot(r1.data, aes(x = intraguild.type, 
                                  y = log(prob_excl_rel))) +
  geom_boxplot() +
  labs(x = "Intra-guild type", y = "log(exclusion ratio)") +
  theme_bw() +
  ggtitle("a)") +
NULL
# exc.dif.p1

exc.dif.p2 <- ggplot(r1.data, aes(x = guild.num, 
                                  y = log(prob_excl_rel))) +
  geom_boxplot() +
  labs(x = "Number of guilds", y = "") +
  theme_bw() +
  ggtitle("b)") +
  NULL
# exc.dif.p2

exc.dif.p3 <- ggplot(r1.data, aes(x = sp.guild, 
                    y = log(prob_excl_rel))) +
  geom_boxplot() +
  labs(x = "Species guild", y = "") +
  theme_bw() +
  ggtitle("c)") +
  NULL
# exc.dif.p3

exc.factors.plot <- exc.dif.p1 + exc.dif.p2 + exc.dif.p3

# ggsave(filename = paste("results/images/Fig_prob_excl_factors.pdf",sep=""),plot = exc.factors.plot,
#        device = cairo_pdf,
#        width = 9,height = 4,dpi = 300)

# -------------------------------------------------------------------------
exc.rich <- lm(log(prob_excl_rel) ~ richness.scaled * intraguild.type * sp.guild, data = r1.data)
# summary(exc.rich)
# broom::tidy(exc.rich)
# DHARMa::testResiduals(exc.rich)
# 
# print(xtable::xtable(as.data.frame(broom::tidy(exc.rich)),
#                      floating=FALSE,
#                      digits = 3,
#                      latex.environments=NULL,
#                      booktabs=FALSE))

# -------------------------------------------------------------------------

# null distributions of probabilities of exclusion

# clean up null distributions for plotting
exc.null$sp.guild <- sp.metrics$sp.guild[match(exc.null$species,sp.metrics$species)]
exc.null.clean <- exc.null %>%
  filter(guilds == "All") %>%
  select(sp.guild, intraguild.type, prob_excl_rel) %>%
  mutate(obs.type = "null")

exc.obs.clean <- exc.obs %>%
  filter(!is.na(sp.guild)) %>%
  select(sp.guild, intraguild.type, prob_excl_rel) %>%
  mutate(obs.type = "observed")
  
exc.combined <- bind_rows(exc.obs.clean,exc.null.clean)

exc.combined.dist <- exc.combined %>%
  group_by(sp.guild,intraguild.type,obs.type) %>%
  summarise(mean.prob.excl = mean(prob_excl_rel),
            sd.prob.excl = sd(prob_excl_rel),
            min.prob.excl = min(prob_excl_rel),
            max.prob.excl = max(prob_excl_rel))

exc.combined.plot <- ggplot(exc.combined, aes(x = obs.type,
                                              y = log(prob_excl_rel))) +
                                              # y = prob_excl_rel)) +
  # geom_boxplot(aes(fill = intraguild.type)) + 
  stat_summary(fun = mean, fun.min = min, fun.max = max, aes(color = intraguild.type)) +
  
  facet_grid(intraguild.type~sp.guild, scales = "free") +
  scale_fill_OkabeIto(darken = .2) +
  scale_color_OkabeIto(darken = .2) +
  labs(x = "",y = "log(exclusion ratio)") +
  theme_bw() +
  theme(strip.background = element_blank())+
  # theme(strip.text = element_blank()) +
  theme(legend.position="none")+
  NULL
# exc.combined.plot

# print(xtable::xtable(exc.combined.dist,
#                      floating=FALSE,
#                      digits = 5,
#                      latex.environments=NULL,
#                      booktabs=FALSE))

# ggsave(filename = paste("results/images/Fig_prob_excl_null.pdf",sep=""),plot = exc.combined.plot,
#        device = cairo_pdf,
#        width = 7,height = 4,dpi = 300)

# -------------------------------------------------------------------------
pd <- .4
exc.rich.plot <- ggplot(r1.data, aes(x = richness,y = log(prob_excl_rel))) +
  geom_point(aes(color = intraguild.type), 
             position = position_jitter(pd), 
             alpha = .4) +
  geom_smooth(aes(group = intraguild.type), 
              colour = "black",size = 2, method = "lm") +
  geom_smooth(aes(color = intraguild.type), size = 1, method = "lm") +
  scale_color_OkabeIto(darken = .2, name = "Intra-guild type") +
  labs(x = "richness", y = "log(exclusion ratio)") +
  theme_bw() +
  NULL
# exc.rich.plot

# ggsave(filename = paste("results/images/Fig_prob_excl_richness.pdf",sep=""),plot = exc.rich.plot,
#        device = cairo_pdf,
#        width = 7,height = 3,dpi = 300)

# -------------------------------------------------------------------------
# subset to be consistent with the community-level analysis

sp.obs <- left_join(exc.obs,sp.metrics) %>%
  filter(intraguild.type == "Resources and \nphenology" &
           guilds == "All") %>%
  select(species,sp.guild,year,plot,prob_excl_rel,diagonal_dominance,
         in_degree,out_degree,avg_intraguild_niche_overlap,avg_interguild_niche_overlap,
         total_intraguild_niche_overlap,total_interguild_niche_overlap)

# scale variables
sp.obs.scaled <- sp.obs %>% select(plot,sp.guild,prob_excl_rel,diagonal_dominance,
                                   total_intraguild_niche_overlap,total_interguild_niche_overlap)
sp.obs.scaled$diagonal_dominance <- scales::rescale(sp.obs.scaled$diagonal_dominance)
sp.obs.scaled$total_intraguild_niche_overlap <- scales::rescale(sp.obs.scaled$total_intraguild_niche_overlap)
sp.obs.scaled$total_interguild_niche_overlap <- scales::rescale(sp.obs.scaled$total_interguild_niche_overlap)
sp.obs.scaled$sp.guild <- factor(sp.obs.scaled$sp.guild, levels = c("plants", "floral visitors", "herbivores"))
# -------------------------------------------------------------------------
# visualization
# hist(sp.obs$diagonal_dominance)
# hist(sp.obs$in_degree)
# hist(sp.obs$out_degree)
# plot(sp.obs$diagonal_dominance,sp.obs$prob_excl_mean)
# plot(sp.obs$diagonal_dominance,log(sp.obs$prob_excl_mean))

# very skewed response
hist(sp.obs.scaled$prob_excl_rel)
sum(sp.obs.scaled$prob_excl_rel == 0)

# log-transformed response
ms <- lmerTest::lmer(log(prob_excl_rel) ~ sp.guild + diagonal_dominance + 
                        total_intraguild_niche_overlap +
                        total_interguild_niche_overlap + (1|plot),data = sp.obs.scaled)
# car::vif(ms)
# summary(ms)
# DHARMa::testResiduals(ms)
# broom.mixed::tidy(ms)
# 
# print(xtable::xtable(as.data.frame(broom.mixed::tidy(ms)),
#                      floating=FALSE,
#                      digits = 3,
#                      latex.environments=NULL,
#                      booktabs=FALSE))
