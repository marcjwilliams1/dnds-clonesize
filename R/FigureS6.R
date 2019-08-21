#load packages, you will need to install these if you don't have them already
library(cowplot)
library(tidyverse)
theme_set(theme_cowplot())
library(jcolors)
library(forcats)
library(ggforce)
library(Hmisc)

library(argparse)

parser <- ArgumentParser(description = "Generate Final Figures")
parser$add_argument('--figure', type='character',
                    help="Outpute figure files")
parser$add_argument('--suppfigures', type='character',
                    help="Outpute figure files", nargs = "+")
parser$add_argument('--skinfitmissense', type='character',
                    help="Fits for missense mutations oesophagus")
parser$add_argument('--skinfitnonsense', type='character',
                    help="Fits for nonsense mutations oesophagus")
parser$add_argument('--skinfitmissensepergene', type='character',
                    help="Fits for missense mutations oesophagus per gene")
parser$add_argument('--skinfitnonsensepergene', type='character',
                    help="Fits for nonsense mutations oesophagus per gene")
parser$add_argument('--mutationcutoff', type='double',
                    help="Only plot genes with at least this number of mutations")
parser$add_argument('--rsqcutoff', type='double',
                    help="Only plot genes with an rsq fit greater than this")
args <- parser$parse_args()
print(args)

message("Generating Figure 3...")
message("\t Reading in data...")

dfnon <- read_csv(args$skinfitnonsense, col_types = cols())
dfmiss <- read_csv(args$skinfitmissense, col_types = cols())

dfnon.gene <- read_csv(args$skinfitnonsensepergene, col_types = cols())
dfmiss.gene <- read_csv(args$skinfitmissensepergene, col_types = cols())

message("Make plots...")

textdf <- dfnon %>%
    filter(	patient == "PD18003") %>%
    mutate(label = paste("list(Delta[fit] == ", round(deltafit, 3), ",R^{2}==",round(rsq, 3) ,")"),
          y1 = min(dnds), y2 = min(dndsfit),
          y = floor(min(y1, y2)) + 0.5) %>%
    distinct(y, deltafit, label)
g1 <- dfnon %>%
    filter(	patient == "PD18003") %>%
    mutate(p = paste0(patient, ", ", Age)) %>%
    ggplot(aes(x = A, y = dnds)) +
    geom_point(alpha = 0.9, col = "plum4", size = 1) +
    geom_line(aes(y = dndsfit), alpha = 0.5, size = 1.5, col = "firebrick") +
    geom_ribbon(aes(ymin = dndsfitlq, ymax = dndsfituq), fill = "firebrick", alpha = 0.3) +
    xlab("Clone area") +
    ylab("Interval dN/dS - Nonsense") +
    ggtitle("PD18003 - All genes") +
    geom_text(data = textdf, aes(label = label, y = y), x = 1.0, size = 5, parse = TRUE)

textdf <- dfmiss %>%
    filter(	patient == "PD18003") %>%
    mutate(label = paste("list(Delta[fit] == ", round(deltafit, 3), ",R^{2}==",round(rsq, 3) ,")"),
          y1 = min(dnds), y2 = min(dndsfit),
          y = floor(min(y1, y2)) + 0.5) %>%
    distinct(y, deltafit, label)
g2 <- dfmiss %>%
    filter(	patient == "PD18003") %>%
    mutate(p = paste0(patient, ", ", Age)) %>%
    ggplot(aes(x = A, y = dnds)) +
    geom_point(alpha = 0.9, col = "plum4", size = 1) +
    geom_line(aes(y = dndsfit), alpha = 0.5, size = 1.5, col = "firebrick") +
    geom_ribbon(aes(ymin = dndsfitlq, ymax = dndsfituq), fill = "firebrick", alpha = 0.3) +
    xlab("Clone area") +
    ylab("Interval dN/dS - Missense") +
    ggtitle("PD18003 - All genes") +
    geom_text(data = textdf, aes(label = label, y = y), x = 1.0, size = 5, parse = TRUE)

textdf <- dfnon %>%
    filter(	patient == "PD20399") %>%
    mutate(label = paste("list(Delta[fit] == ", round(deltafit, 3), ",R^{2}==",round(rsq, 3) ,")"),
          y1 = min(dnds), y2 = min(dndsfit),
          y = floor(min(y1, y2)) + 0.5) %>%
    distinct(y, deltafit, label)
g3 <- dfnon %>%
    filter(	patient == "PD20399") %>%
    mutate(p = paste0(patient, ", ", Age)) %>%
    ggplot(aes(x = A, y = dnds)) +
    geom_point(alpha = 0.9, col = "plum4", size = 1) +
    geom_line(aes(y = dndsfit), alpha = 0.5, size = 1.5, col = "firebrick") +
    geom_ribbon(aes(ymin = dndsfitlq, ymax = dndsfituq), fill = "firebrick", alpha = 0.3) +
    xlab("Clone area") +
    ylab("Interval dN/dS - Nonsense") +
    ggtitle("PD20399 - All genes") +
    geom_text(data = textdf, aes(label = label, y = y), x = 1.0, size = 5, parse = TRUE)


textdf <- dfmiss.gene %>%
    filter(	patient == "PD18003", gene == "TP53") %>%
    mutate(label = paste("list(Delta[fit] == ", round(deltafit, 3), ",R^{2}==",round(rsq, 3) ,")"),
          y1 = min(dnds), y2 = min(dndsfit),
          y = floor(min(y1, y2)) + 1) %>%
    distinct(y, deltafit, label)
g4 <- dfnon.gene %>%
    filter(patient == "PD20399", gene == "NOTCH1") %>%
    ggplot(aes(x = A, y = dnds)) +
    geom_point(alpha = 0.9, size = 1) +
    geom_line(aes(y = dndsfit), alpha = 0.5, size = 1.5, col = "plum4") +
    geom_ribbon(aes(ymin = dndsfitlq, ymax = dndsfituq), fill = "plum4", alpha = 0.3) +
    xlab("Clone area") +
    ylab("Interval dN/dS - Nonsense") +
    ggtitle("PD20399 - NOTCH1") +
    geom_text(data = textdf, aes(label = label, y = y), x = 1.0, size = 5, parse = TRUE)

textdf <- dfmiss.gene %>%
    filter(	patient == "PD18003", gene == "NOTCH1") %>%
    mutate(label = paste("list(Delta[fit] == ", round(deltafit, 3), ",R^{2}==",round(rsq, 3) ,")"),
          y1 = min(dnds), y2 = min(dndsfit),
          y = floor(min(y1, y2)) + 1) %>%
    distinct(y, deltafit, label)
g5 <- dfmiss.gene %>%
    filter(patient == "PD18003", gene == "NOTCH1") %>%
    ggplot(aes(x = A, y = dnds)) +
    geom_point(alpha = 0.9, size = 1) +
    geom_line(aes(y = dndsfit), alpha = 0.5, size = 1.5, col = "plum4") +
    geom_ribbon(aes(ymin = dndsfitlq, ymax = dndsfituq), fill = "plum4", alpha = 0.3) +
    xlab("Clone area") +
    ylab("Interval dN/dS - Nonsense") +
    ggtitle("PD18003 - NOTCH1") +
    geom_text(data = textdf, aes(label = label, y = y), x = 1.0, size = 5, parse = TRUE)

message("Combine data frames...")

dfcombined <- bind_rows(mutate(dfmiss.gene, mutationtype = "Missense"),
                        mutate(dfnon.gene, mutationtype = "Nonsense"))
print(distinct(dfcombined, patient, gene, nmutations))

DFresults <- dfcombined %>%
    group_by(gene, mutationtype, deltafit, nmutations, deltafitlq, deltafituq) %>%
    summarise(rsq = first(rsq))

print(DFresults)

DFresults <- DFresults %>%
    filter(rsq > args$rsqcutoff, nmutations > args$mutationcutoff) %>%
    group_by(gene, mutationtype) %>% mutate(n = n())

message("Plot DFE...")

DFEmiss <- DFresults %>%
    filter(mutationtype == "Missense") %>%
    mutate(genelab = ifelse(gene == "TP53" | gene == "NOTCH1" | gene == "NOTCH2", gene, "Other genes")) %>%
    mutate(genelab = factor(genelab, levels = c("TP53", "NOTCH1", "NOTCH2", "Other genes"))) %>%
    ggplot(aes(x = deltafit, fill = genelab)) +
    geom_histogram(bins = 25, position  = "stack", color = "white", size = 0.2) +
    ggtitle("Missense") +
    xlim(c(0, 0.1)) +
    xlab(expression(paste(Delta))) +
    ylab("Counts") +
    theme(legend.position = c(0.5, 0.8)) +
    scale_fill_jcolors(pal = "pal7",  name = "")

g <- plot_grid(g1, g2, g3, g4, g5, DFEmiss, ncol = 3, labels = c("a", "b", "c", "d", "e", "f"))

save_plot(args$figure, g, base_height = 7, base_width = 10)


message("Plot all fits...")


textdf <- dfcombined %>%
    filter(rsq > args$rsqcutoff, nmutations > args$mutationcutoff, mutationtype == "Missense") %>%
    mutate(label = paste("list(Delta[fit] == ", round(deltafit, 3), ",R^{2}==",round(rsq, 3) ,")"),
          y1 = min(dnds), y2 = min(dndsfit),
          y = floor(min(y1, y2)) + 1) %>%
    distinct(y, deltafit, label)
print(textdf)
gmiss <- dfcombined %>%
    filter(rsq > args$rsqcutoff, nmutations > args$mutationcutoff, mutationtype == "Missense") %>%
    ggplot(aes(x = A, y = dnds)) +
    geom_point(alpha = 0.9, size = 1) +
    geom_line(aes(y = dndsfit), alpha = 0.5, size = 1.0, col = "plum4") +
    geom_ribbon(aes(ymin = dndsfitlq, ymax = dndsfituq), fill = "plum4", alpha = 0.3) +
    xlab("Clone area") +
    ylab("Interval dN/dS - Missense") +
    facet_wrap(gene~patient, scale = "free_y")
    geom_text(data = textdf, aes(label = label, y = y), x = 3.0, size = 5, parse = TRUE)

#save_plot("FinalFigures/Figure3SX.pdf", g3, base_height = 12, base_width = 12)

textdf <- dfcombined %>%
    filter(rsq > args$rsqcutoff, nmutations > 4, mutationtype == "Nonsense") %>%
    mutate(label = paste("list(Delta[fit] == ", round(deltafit, 3), ",R^{2}==",round(rsq, 3) ,")"),
          y1 = min(dnds), y2 = min(dndsfit),
          y = floor(min(y1, y2)) + 1) %>%
    distinct(y, deltafit, label)
print(textdf)
gnon <- dfcombined %>%
    filter(rsq > args$rsqcutoff, nmutations > 4, mutationtype == "Nonsense") %>%
    ggplot(aes(x = A, y = dnds)) +
    geom_point(alpha = 0.9, size = 1) +
    geom_line(aes(y = dndsfit), alpha = 0.5, size = 1.0, col = "darkseagreen4") +
    geom_ribbon(aes(ymin = dndsfitlq, ymax = dndsfituq), fill = "darkseagreen4", alpha = 0.3) +
    xlab("Clone area") +
    ylab("Interval dN/dS - Nonsense") +
    facet_wrap(gene~patient, scale = "free_y")
    geom_text(data = textdf, aes(label = label, y = y), x = 3.0, size = 5, parse = TRUE)

g <- plot_grid(gnon, gmiss, ncol = 1)

save_plot(args$suppfigures[1], g, base_height = 24, base_width = 12)
