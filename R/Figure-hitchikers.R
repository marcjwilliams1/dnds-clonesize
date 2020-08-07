library(tidyverse)
library(jcolors)
library(argparse)
library(tidybayes)
library(bayesplot)
library(modelr)
library(cowplot)

dfsims <- read_csv(snakemake@input$simulationdata)

dfdata <- read_csv(snakemake@input$oesophagusdata)
donor <- readxl::read_xlsx(snakemake@input$oesophagusmetadata, skip = 1) %>%
  dplyr::rename(donor = PD)
dfdata <- left_join(dfdata, donor)

# Plot data per patch
message("Read in data from xlsx file")

dfall <- readxl::read_xlsx(snakemake@input$oesophagusdata_all, sheet = 1, skip = 16)

dfallcounts <- dfall %>%
  mutate(muttype = case_when(
    str_detect(gene, "NOTCH1") & impact != "Synonymous" ~ "NOTCH1",
    str_detect(gene, "TP53") & impact != "Synonymous" ~ "TP53",
    impact == "Synonymous" ~ "Synonymous"
  )) %>%
  filter(!is.na(muttype)) %>%
  group_by(sampleID, muttype) %>%
  summarise(vaf = mean(vaf), n = n()) %>%
  ungroup() %>%
  complete(sampleID, nesting(muttype), fill = list(vaf = 0.0, n = 0))

nmuts <- dfallcounts %>%
  dplyr::select(-vaf) %>%
  spread(muttype, n)

vafmuts <- dfallcounts %>%
  dplyr::select(-n) %>%
  spread(muttype, vaf)

patch <- full_join(nmuts, vafmuts, by = "sampleID", suffix = c("_n", "_vaf"))

(g1 <- patch %>%
  mutate(NOTCH1_n = ifelse(NOTCH1_n > 9, 10, NOTCH1_n),
         TP53_n = ifelse(TP53_n > 9, 10, TP53_n)) %>%
  ggplot(aes(x = factor(NOTCH1_n), y = Synonymous_n)) +
  geom_boxplot(width = 0.2, alpha = 0.8) +
  geom_violin(alpha = 0.5, fill = "firebrick4") +
  theme_cowplot() +
  xlab("# NOTCH1 Non-synonymous mutations") +
  ylab("# Synonymous mutations"))

(g2 <- patch %>%
    mutate(NOTCH1_n = ifelse(NOTCH1_n > 9, 10, NOTCH1_n),
           TP53_n = ifelse(TP53_n > 9, 10, TP53_n)) %>%
    ggplot(aes(x = factor(TP53_n), y = Synonymous_n)) +
    geom_boxplot(width = 0.2, alpha = 0.8) +
    geom_violin(alpha = 0.5, fill = "firebrick4") +
    theme_cowplot() +
    xlab("# TP53 Non-synonymous mutations") +
    ylab(""))

(g3 <- patch %>%
    mutate(NOTCH1_n = ifelse(NOTCH1_n > 9, 10, NOTCH1_n),
           TP53_n = ifelse(TP53_n > 9, 10, TP53_n)) %>%
    ggplot(aes(x = factor(NOTCH1_n), y = Synonymous_vaf)) +
    geom_boxplot(width = 0.2, alpha = 0.8) +
    geom_violin(alpha = 0.5, fill = "deepskyblue4") +
    theme_cowplot() +
    scale_y_log10() +
    xlab("# NOTCH1 Non-synonymous mutations") +
    ylab("VAF synonymous mutations"))

(g4 <- patch %>%
    mutate(NOTCH1_n = ifelse(NOTCH1_n > 9, 10, NOTCH1_n),
           TP53_n = ifelse(TP53_n > 9, 10, TP53_n)) %>%
    ggplot(aes(x = factor(TP53_n), y = Synonymous_vaf)) +
    geom_boxplot(width = 0.2, alpha = 0.8) +
    geom_violin(alpha = 0.5, fill = "deepskyblue4") +
    theme_cowplot() +
    scale_y_log10() +
    xlab("# TP53 Non-synonymous mutations") +
    ylab(""))

gall <- plot_grid(g1, g2, g3, g4, ncol = 2, align = "hv", labels = c("a", "b", "c", "d"))

message("Save figure")
save_plot(snakemake@output$suppfigures[1], gall, base_height = 8, base_width = 20)

mylm <- lm(Synonymous_n ~ NOTCH1_n, patch)
print(summary(mylm))

mylm <- lm(Synonymous_n ~ TP53_n, patch)
print(summary(mylm))

mylm <- lm(Synonymous_vaf ~ NOTCH1_n, patch)
print(summary(mylm))

mylm <- lm(Synonymous_vaf ~ TP53_n, patch)
print(summary(mylm))



dfallcounts <- dfall %>%
  mutate(muttype = case_when(
    impact != "Synonymous" ~ "NS",
    impact == "Synonymous" ~ "Synonymous"
  )) %>%
  filter(!is.na(muttype)) %>%
  group_by(sampleID, muttype) %>%
  summarise(vaf = mean(vaf), n = n()) %>%
  ungroup() %>%
  complete(sampleID, nesting(muttype), fill = list(vaf = 0.0, n = 0))

nmuts <- dfallcounts %>%
  dplyr::select(-vaf) %>%
  spread(muttype, n)

vafmuts <- dfallcounts %>%
  dplyr::select(-n) %>%
  spread(muttype, vaf)

patch <- full_join(nmuts, vafmuts, by = "sampleID", suffix = c("_n", "_vaf"))


(g1 <- patch %>%
    ggplot(aes(x = NS_n, y = Synonymous_vaf, group = NS_n)) +
    geom_boxplot(width = 0.2, alpha = 0.8) +
    geom_violin(alpha = 0.5, fill = "firebrick4") +
    theme_cowplot() +
    xlab("# NOTCH1 Non-synonymous mutations") +
    ylab("# Synonymous mutations"))

mylm <- lm(Synonymous_n ~ NS_n, patch)
print(summary(mylm))

mylm <- lm(Synonymous_vaf ~ NS_n, patch)
print(summary(mylm))
