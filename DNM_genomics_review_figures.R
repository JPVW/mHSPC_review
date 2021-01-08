##### Molecular Characterization of Metastatic Castration-Sensitive Prostate Cancer #####
library(ggplot2)
library(tidyverse)
library(cowplot)
library(gtable)
library(gridExtra)
library(grid)
library(ggstatsplot)
setwd("path")

pathways <- read.delim("pathways.txt")
pathways <- pathways[,-3]

# Function to extract legend from the plot
g_legend <- function(a.gplot){
  tmp <- ggplot_gtable(ggplot_build(a.gplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)}




### Figure 2: Disease state
disease_state <- read.csv2("disease_state.csv", header = T, sep = ";")
disease_state$state <- factor(disease_state$state, levels = c("HSPC-Loc", "HSPC-M1", "CRPC-M1"))
names(disease_state)[names(disease_state) == "Gene"] <- "gene"
disease_state <- left_join(disease_state, pathways, by ='gene')

# adding of weighted means to boxplot
# calculate weight of every study, do by absolute numbers per group lowest absolute number = 1, calculate ratio
disease_state$weight_all <- disease_state$number_state/16
disease_state$weight_all <- as.numeric(disease_state$weight_all)

CI_disease_state <- disease_state %>%
  group_by(state, gene) %>%
  dplyr::summarise(mean = mean(AF, na.rm = TRUE),
                   sd = sd(AF, na.rm = TRUE),
                   n = n()) %>%
  mutate(se = sd / sqrt(n),
         lower.ci = mean - qt(1 - (0.05 / 2), n - 1) * se,
         upper.ci = mean + qt(1 - (0.05 / 2), n - 1) * se)

CI_disease_state_weight <- disease_state %>%
  group_by(state, gene) %>%
  dplyr::summarise(weighted.mean = weighted.mean(AF,weight_all, na.rm = TRUE),
                   sd = sd(AF, na.rm = TRUE),
                   n = n()) %>%
  mutate(se = sd / sqrt(n),
         lower.ci = weighted.mean - qt(1 - (0.05 / 2), n - 1) * se,
         upper.ci = weighted.mean + qt(1 - (0.05 / 2), n - 1) * se)


# Plot figure 2
CI_disease_state_weight <- left_join(CI_disease_state_weight, pathways, by ='gene')

ggplot(disease_state, aes(gene, AF)) +
  geom_boxplot(fatten=NULL, aes(color = state), outlier.shape = NA) +
  geom_point(aes(shape = sample_type, color = state, group = state), 
             alpha = 0.6, position = position_dodge(width = 0.75)) +
  geom_point(data = CI_disease_state_weight, aes(gene, weighted.mean, color = state, fill = state), 
             shape = 124, size = 2, show.legend = F, position = position_dodge(width = 0.75)) +
  facet_grid(pathway~., scales = "free", space = "free",
             labeller=label_wrap_gen(width=10), switch = "y") +
  scale_color_manual(values = c("steelblue", "springgreen4","red3")) +
  scale_size_continuous(range = c(1,4)) +
  scale_y_continuous(breaks = c(0,5,10,15,20,30,40,50,60,70,80)) +
  coord_flip() +
  theme_linedraw() +
  scale_shape_manual(values = c(15,16,17)) +
  labs(shape = "Sample type", color = "Disease state", size = "Patients included in the study") + 
  xlab("Genes") + ylab("Alteration frequency (%)") +
  theme(plot.title = element_text(size = 20, face = "bold")) +
  ggtitle("Disease state")






### Figure 3: Disease volume + disease timing
## Disease volume
disease_volume <- read.csv2("disease_volume.csv", header = T, sep = ";")
disease_volume <- left_join(disease_volume, pathways, by ='gene')

# Calculate weighted means
disease_volume$weight_all <- disease_volume$absolute_number/14


CI_disease_volume <- disease_volume %>%
  group_by(disease_volume, gene) %>%
  dplyr::summarise(mean = mean(frequency, na.rm = TRUE),
                   sd = sd(frequency, na.rm = TRUE),
                   n = n()) %>%
  mutate(se = sd / sqrt(n),
         lower.ci = mean - qt(1 - (0.05 / 2), n - 1) * se,
         upper.ci = mean + qt(1 - (0.05 / 2), n - 1) * se)

CI_disease_volume_weight <- disease_volume %>%
  group_by(disease_volume, gene) %>%
  dplyr::summarise(weighted.mean = weighted.mean(frequency,weight_all, na.rm = TRUE),
                   sd = sd(frequency, na.rm = TRUE),
                   n = n()) %>%
  mutate(se = sd / sqrt(n),
         lower.ci = weighted.mean - qt(1 - (0.05 / 2), n - 1) * se,
         upper.ci = weighted.mean + qt(1 - (0.05 / 2), n - 1) * se)

## Disease timing
disease_timing <- read.csv2(file = "disease_timing.csv", header = T, sep = ";")

# Calculate weighted mean
disease_timing$weight_all <- disease_timing$total/41

CI_disease_timing <- disease_timing %>%
  group_by(state, Gene) %>%
  dplyr::summarise(mean = mean(percentage, na.rm = TRUE),
                   sd = sd(percentage, na.rm = TRUE),
                   n = n()) %>%
  mutate(se = sd / sqrt(n),
         lower.ci = mean - qt(1 - (0.05 / 2), n - 1) * se,
         upper.ci = mean + qt(1 - (0.05 / 2), n - 1) * se)

CI_disease_timing_weight <- disease_timing %>%
  group_by(state, Gene) %>%
  dplyr::summarise(weighted.mean = weighted.mean(percentage,weight_all, na.rm = TRUE),
                   sd = sd(percentage, na.rm = TRUE),
                   n = n()) %>%
  mutate(se = sd / sqrt(n),
         lower.ci = weighted.mean - qt(1 - (0.05 / 2), n - 1) * se,
         upper.ci = weighted.mean + qt(1 - (0.05 / 2), n - 1) * se)






# Plot left panel of figure 3
f3a <- 
  ggplot(disease_volume,aes(disease_volume, frequency, color = disease_volume)) + 
  geom_boxplot(fatten = NULL) + 
  geom_point(data = CI_disease_volume_weight, 
             aes(disease_volume, weighted.mean, fill = disease_volume, 
                 color = disease_volume), shape = 124, size = 4, show.legend = F) +
  geom_point(aes(color = disease_volume, size = absolute_number, shape = sample_type), 
             alpha = 0.6) + 
  scale_shape_manual(values = c(15,16,17)) +
  scale_color_manual(name = 'Disease volume', values =c('red', 'blue')) + 
  theme_linedraw() + 
  facet_grid(gene~., scales = 'free', space = 'free', switch = "y") +
  theme(strip.text.y.left = element_text(angle=0, size = 15)) +
  theme(axis.text.y = element_blank(), axis.title.x = element_text(size = 15),axis.text.x = element_text(size = 15)) + xlab('') + ylab('Alteration frequency (%)') + 
  ylim(0,50) +
  coord_flip() + 
  theme(strip.text.y = element_text(angle = 180)) + 
  theme(axis.ticks.y = element_blank()) +
  labs(shape = "Sample type", size = "Patients included in study") +
  guides(color = guide_legend(reverse = T)) +
  ggtitle("Disease volume") + 
  theme(plot.title = element_text(size = 20, face = "bold")) +
  theme(legend.position = "none") # Don't include if legend needs to be extracted

leg3a <- g_legend(f3a)


# Plot right panel of figure 3
f3b <-
  ggplot(disease_timing,aes(state, percentage, color = state)) + 
  geom_boxplot(fatten=NULL) + 
  geom_point(aes(color = state, size = total, shape = sample_type), alpha = 0.6) + 
  geom_point(data = CI_disease_timing_weight, 
             aes(state, weighted.mean, fill = state, color = state), 
             shape = 124, size = 1, show.legend = F, stroke = 5) +
  scale_color_manual(name = 'Timing of metastasis', values =c('gold2', 'darkorchid4')) + 
  theme_linedraw() + 
  facet_grid(Gene~., scales = 'free', space = 'free', switch = "y") + 
  theme(strip.background = element_blank(), strip.text.x = element_blank()) +
  theme(axis.text.y = element_blank(), axis.title.x = element_text(size = 15), axis.text.x = element_text(size = 15)) + xlab('') + ylab('Alteration frequency (%)') + 
  ylim(0,50) +
  coord_flip() + 
  theme(strip.text.y = element_text(angle = 180)) + 
  theme(axis.ticks.y = element_blank()) +
  labs(shape = "Sample type", size = "Patients included in study") +
  guides(size = F, color = guide_legend(reverse = T)) +
  ggtitle("Timing of metastases") +
  theme(plot.title = element_text(size = 20, face = "bold")) +
  theme(legend.position = "none")


leg3b <- g_legend(f3b)

  # Combine figures
figure3 <- plot_grid(f3a, f3b, ncol = 2, nrow = 1)
legend3 <- plot_grid(leg3a, leg3b, ncol = 1, rel_heights = c(1.5,0.5)) 

