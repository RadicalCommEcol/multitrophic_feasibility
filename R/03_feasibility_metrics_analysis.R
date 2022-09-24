
# statistical analyses at the community level

# INPUT
# feasibility domains: "results/feasibility_domain_observed.csv"
# community-level metrics: "results/community_metrics.csv", 
# "results/community_metrics_null.csv"

# OUTPUT
# Tables of the paper: 
#   feasibility domain - number of guilds: inline results in main text
#   feasibility domain - richness: Table S2
#   feasibility domain - metrics: Table 1

# Figures
#   feasiblity domain - metrics: Figure 5
#   feasibility domain rankings: Figure S1

# -------------------------------------------------------------------------

library(tidyverse)
library(broom)
library(broom.mixed)
library(cowplot)
library(corrplot)
library(ggfortify)
library(DHARMa)
library(colorblindr)
library(effects)
library(patchwork)
library(lmerTest)

# -------------------------------------------------------------------------
# read data

vers <- ""
# vers <- "_v2"

# feasibility domain of the observed data
fd.obs <- read.csv2("results/feasibility_domain_observed.csv")
# forgot an "s"
names(fd.obs)[which(names(fd.obs) == "guild")] <- "guilds"

# observed network metrics
metrics.obs <- read.csv2("results/community_metrics.csv")
metrics.null <- read_csv2("results/community_metrics_null.csv")

# -------------------------------------------------------------------------
# clean up names

fd.obs$guilds <- factor(fd.obs$guilds,levels = c("plants",
                                         "floral visitors",
                                         "herbivores",
                                         "plants-floral visitors",
                                         "plants-herbivores","all")) 
fd.obs$guilds <- dplyr::recode(fd.obs$guilds, "plants" = "Plants", 
                    "floral visitors" = "Pollinators",
                    "herbivores" = "Herbivores",
                    "plants-floral visitors" = "Plants-Pollinators",
                    "plants-herbivores" = "Plants-Herbivores",
                    "all" = "All")

# fd.obs$plot <- as.factor(fd.obs$plot)

fd.obs$intraguild.type[fd.obs$intraguild.type == "mean_field"] <- "Mean field"
fd.obs$intraguild.type[fd.obs$intraguild.type == "nesting_larvae"] <- "Resources"
fd.obs$intraguild.type[fd.obs$intraguild.type == "nesting_larvae_phenology"] <- "Resources and \nphenology"

# -------------------------------------------------------------------------
metrics.obs$guilds <- factor(metrics.obs$guilds,levels = c("plants",
                                                 "floral visitors",
                                                 "herbivores",
                                                 "plants-floral visitors",
                                                 "plants-herbivores","all")) 
metrics.obs$guilds <- dplyr::recode(metrics.obs$guilds, "plants" = "Plants", 
                        "floral visitors" = "Pollinators",
                        "herbivores" = "Herbivores",
                        "plants-floral visitors" = "Plants-Pollinators",
                        "plants-herbivores" = "Plants-Herbivores",
                        "all" = "All")

# metrics.obs$plot <- as.factor(metrics.obs$plot)

metrics.obs$intraguild.type[metrics.obs$intraguild.type == "mean_field"] <- "Mean field"
metrics.obs$intraguild.type[metrics.obs$intraguild.type == "nesting_larvae"] <- "Resources"
metrics.obs$intraguild.type[metrics.obs$intraguild.type == "nesting_larvae_phenology"] <- "Resources and \nphenology"

metrics.null$guilds <- factor(metrics.null$guilds,levels = c("plants",
                                                           "floral visitors",
                                                           "herbivores",
                                                           "plants-floral visitors",
                                                           "plants-herbivores","all")) 
metrics.null$guilds <- dplyr::recode(metrics.null$guilds, "plants" = "Plants", 
                             "floral visitors" = "Pollinators",
                             "herbivores" = "Herbivores",
                             "plants-floral visitors" = "Plants-Pollinators",
                             "plants-herbivores" = "Plants-Herbivores",
                             "all" = "All")

# metrics.null$plot <- as.factor(metrics.null$plot)

metrics.null$intraguild.type[metrics.null$intraguild.type == "mean_field"] <- "Mean field"
metrics.null$intraguild.type[metrics.null$intraguild.type == "nesting_larvae"] <- "Resources"
metrics.null$intraguild.type[metrics.null$intraguild.type == "nesting_larvae_phenology"] <- "Resources and \nphenology"

# -------------------------------------------------------------------------
# first question: fd ~ number of guilds/richness

fd.obs$guild.num <- "one"
fd.obs$guild.num[fd.obs$guilds %in% c("Plants-Pollinators",
                                      "Plants-Herbivores")] <- "two"
fd.obs$guild.num[fd.obs$guilds == "All"] <- "three"
fd.obs$guild.num <- factor(fd.obs$guild.num,levels = c("one","two","three"))

richness.data <- metrics.obs %>% 
  filter(metric == "richness") %>% 
  pivot_wider(names_from = "metric",values_from = "value")
richness.data$richness.scaled <- scales::rescale(richness.data$richness)

r1.data <- left_join(fd.obs,richness.data)

# -------------------------------------------------------------------------
fd.dif <- lm(omega_mean ~ guild.num * intraguild.type, data = r1.data)
fd.aov3 <- car::Anova(fd.dif,type = 3)

# -------------------------------------------------------------------------
fd.rich <- lm(omega_mean ~ richness.scaled * intraguild.type, data = r1.data)
# summary(fd.rich)
# DHARMa::testResiduals(fd.rich)

# print(xtable::xtable(as.data.frame(tidy(fd.rich)),floating=FALSE, 
#       digits = 5,
#       latex.environments=NULL,
#       booktabs=FALSE))

# -------------------------------------------------------------------------
# are rankings maintained? i.e. is the ranking in feasibility domains 
# consistent across parameterizations?

fd.rank <- fd.obs %>% 
  select(year,plot,guilds,intraguild.type,omega_mean) %>%
  mutate(community_id = as.character(paste(plot,"_",year))) %>%
  group_by(intraguild.type,guilds) %>%
  mutate(rank = factor(rank(omega_mean))) %>%
  filter(guilds == "All" & intraguild.type %in% c("Mean field",
                                                  "Resources and \nphenology"))
fd.rank$community_id <- as.factor(fd.rank$community_id)

pj <- .2
fd.rank.plot <- ggplot(fd.rank,aes(x = community_id, y = rank, group = intraguild.type)) + 
  geom_line(aes(color = intraguild.type)) +
  geom_point(aes(fill = intraguild.type), 
             # position = position_jitter(pj), 
             shape = 21) +
  # facet_grid(year~.) +
  scale_fill_manual(values = c("#51A48C","#CBA34D"), name = "Intraguild \ncompetition") +
  scale_color_manual(values = c("#51A48C","#CBA34D"), name = "Intraguild \ncompetition") +
  ylab("Feasibility domain ranking")+
  xlab("Community id (plot and year)") +
  theme_bw() +
  theme(strip.background = element_blank())+
  NULL
# fd.rank.plot

# ggsave(filename = "results/images/Fig_plot_ranking.pdf",plot = fd.rank.plot,
#        device = cairo_pdf,
#        width = 12,height = 4,dpi = 300)

cor.data <- fd.rank %>% 
  select(-omega_mean,-community_id) %>%
  mutate(rank = as.numeric(rank)) %>%
  pivot_wider(names_from = intraguild.type,values_from = rank) 
names(cor.data)[4:5] <- c("mf","res_ph")

# cor.test(cor.data$mf,cor.data$res_ph,method = "spearman")

# cor.results <- cor.data %>%
#   # group_by(year) %>%
#   summarise(correlation = cor(mf,res_ph,method = "spearman"))

# -------------------------------------------------------------------------
# tidy metrics for the models
# note that I also obtain z-scores, this is a legacy from preliminary analyses
# in the final version, I use raw metrics

metrics.null.long <- metrics.null %>% 
  # number of diag. dom. species is not normally distributed
  select(-diagonally_dominant_sp,-complexity,-skewness,-intra_inter_ratio) %>%
  pivot_longer(richness:quant_modularity,names_to = "metric",values_to = "value") %>%
  group_by(intraguild.type,guilds,null.model,metric) %>%
  summarise(mean.value = mean(value, na.rm = T),
            sd.value = sd(value, na.rm = T)) 
# fuck it
metrics.null.long$null.model[metrics.null.long$null.model == "diagonal dominance"] <- "diag.dom"

metrics.null.wide <- metrics.null.long %>% 
  select(metric,intraguild.type,guilds,null.model,mean.value,sd.value) %>%
  pivot_wider(names_from = null.model,values_from =c(mean.value,sd.value))

metrics.obs.z <- left_join(metrics.obs,metrics.null.wide) %>%
  na.omit() %>% mutate(dd = (value-mean.value_diag.dom)/sd.value_diag.dom,
                       topo = (value-mean.value_topology)/sd.value_topology) %>%
  select(-(mean.value_diag.dom:sd.value_topology)) %>%
  pivot_longer(cols = c(dd,topo),names_to = "null_model",values_to = "z_score") %>%
  mutate(z_score = replace(z_score, z_score > 1e3 | z_score < -1e3, NA))

# -------------------------------------------------------------------------

# final subset
metrics.z <- subset(metrics.obs.z,metric %in% c("avg_diagonal_dominance",
                                                "avg_intraguild_niche_overlap",
                                                "avg_interguild_niche_overlap") &
                      null_model == "topo" &
                      guilds == "All" &
                      intraguild.type == "Resources and \nphenology")

metrics.raw.fd <- left_join(metrics.z,fd.obs) %>%
  select(year,plot,omega_mean,omega_isotropic,isotropy_index_mean,metric,value) %>%
  pivot_wider(names_from = metric,values_from = value) %>%
  group_by(year,plot) %>%
  mutate(net = cur_group_id())

# -------------------------------------------------------------------------
# models

# scale variables
metrics.scaled <- metrics.raw.fd %>% select(plot,omega_mean,avg_diagonal_dominance:avg_interguild_niche_overlap)
metrics.scaled$avg_diagonal_dominance <- scales::rescale(metrics.scaled$avg_diagonal_dominance)
metrics.scaled$avg_intraguild_niche_overlap <- scales::rescale(metrics.scaled$avg_intraguild_niche_overlap)
metrics.scaled$avg_interguild_niche_overlap <- scales::rescale(metrics.scaled$avg_interguild_niche_overlap)

m2 <- lmerTest::lmer(omega_mean ~ avg_diagonal_dominance + avg_intraguild_niche_overlap + avg_interguild_niche_overlap + (1|plot),
                     data = metrics.scaled)
# car::vif(m2)
# broom.mixed::tidy(m2)
## the residuals look good
# DHARMa::testResiduals(m2)

# print(xtable::xtable(as.data.frame(broom.mixed::tidy(m2)),
#                      floating=FALSE, 
#                      digits = 5,
#                      latex.environments=NULL,
#                      booktabs=FALSE))
# -------------------------------------------------------------------------

edd <- effects::effect("avg_diagonal_dominance",m2,xlevels = 20)
eintra <- effects::effect("avg_intraguild_niche_overlap",m2,xlevels = 20)
# plot(ec)
# plot(ed)

edd.df <- data.frame(edd)
eintra.df <- data.frame(eintra)

dd.plot <- ggplot(edd.df,aes(y = fit, x = avg_diagonal_dominance)) +
  geom_ribbon(aes(ymin = lower, ymax = upper), fill = "grey", alpha = .3) +
  geom_line()+
  geom_point(data = metrics.scaled,aes(x = avg_diagonal_dominance, y = omega_mean))+
  theme_bw() +
  labs(x="Average diagonal dominance (scaled)",y = "Feasibility domain") +
  ylim(c(0.27,0.41)) +
  ggtitle(label = "A)") +
  NULL

intra.plot <- ggplot(eintra.df,aes(y = fit, x = avg_intraguild_niche_overlap)) +
  geom_ribbon(aes(ymin = lower, ymax = upper), fill = "grey", alpha = .3) +
  geom_line()+
  geom_point(data = metrics.scaled,aes(x = avg_intraguild_niche_overlap, y = omega_mean))+
  theme_bw() +
  labs(x="Average intra-guild interaction overlap (scaled)",y = "") +
  ylim(c(0.27,0.4)) +
  ggtitle(label = "B)") +
  NULL

metrics.plot <- dd.plot + intra.plot

# ggsave(filename = "results/images/Fig_metrics.pdf",plot = metrics.plot,
#        device = cairo_pdf,
#        width = 8,height = 4,dpi = 300)

# -------------------------------------------------------------------------




