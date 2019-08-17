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
parser$add_argument('--oesophagusfitmissense', type='character',
                    help="Fits for missense mutations oesophagus")
parser$add_argument('--oesophagusfitnonsense', type='character',
                    help="Fits for nonsense mutations oesophagus")
parser$add_argument('--oesophagusfitmissensepergene', type='character',
                    help="Fits for missense mutations oesophagus per gene")
parser$add_argument('--oesophagusfitnonsensepergene', type='character',
                    help="Fits for nonsense mutations oesophagus per gene")
args <- parser$parse_args()

message("Generating Figure 3...")
message("\t Reading in data...")

dfnon <- read_csv(args$oesophagusfitnonsense, col_types = cols())
dfmiss <- read_csv(args$oesophagusfitmissense, col_types = cols())

dfnon.gene <- read_csv(args$oesophagusfitnonsensepergene, col_types = cols())
dfmiss.gene <- read_csv(args$oesophagusfitnonsensepergene, col_types = cols())

message("Combine data frames...")
dfcombined <- bind_rows(mutate(dfmiss.gene, mutationtype = "Missense"),
                        mutate(dfnon.gene, mutationtype = "Nonsense"))

library(ggforce)

DFresults <- dfcombined %>%
    group_by(gene, mutationtype, deltafit, nmutations, deltafitlq, deltafituq) %>%
    summarise(rsq = first(rsq)) %>%
    filter(rsq > 0.6, nmutations > 7) %>%
    group_by(gene, mutationtype) %>% mutate(n = n()) %>% filter(n > 1)

message("Summarise per gene values...")

g1del <- DFresults %>%
    filter(mutationtype == "Missense") %>%
    group_by(gene) %>%
    mutate(s = median(deltafit)) %>%
    ungroup() %>%
    ggplot(aes(x = fct_reorder(gene, s, .desc = TRUE), y = deltafit)) +
    geom_point(col = "firebrick") +
    geom_hline(yintercept = 0.0, lty = 2) +
    geom_boxplot(alpha = 0.1, width = 0.3, fill = "firebrick") + #scale_y_log10() +
    xlab("") +
    ylab(expression(paste(Delta))) +
    coord_flip() + ggtitle("Missense")

g2del <- DFresults %>%
    filter(mutationtype == "Nonsense") %>%
    group_by(gene) %>%
    mutate(s = median(deltafit)) %>%
    ungroup() %>%
    ggplot(aes(x = fct_reorder(gene, s, .desc = TRUE), y = deltafit)) +
    geom_point(col = "deepskyblue4") +
    geom_hline(yintercept = 0.0, lty = 2) +
    geom_boxplot(alpha = 0.1, width = 0.3, fill = "deepskyblue4") + #scale_y_log10() +
    xlab("") +
    ylab(expression(paste(Delta))) +
    coord_flip() + ggtitle("Nonsense")

geneplot <- plot_grid(g1del, g2del, align = T, labels = c("a", "b"))

message("Plot DFE...")


DFresults <- dfcombined %>%
    group_by(gene, mutationtype, deltafit, nmutations, deltafitlq, deltafituq) %>%
    summarise(rsq = first(rsq)) %>%
    filter(rsq > 0.6, nmutations > 7) %>%
    group_by(gene, mutationtype) %>% mutate(n = n())

DFEmiss <- DFresults %>%
    filter(mutationtype == "Missense") %>%
    mutate(genelab = ifelse(gene == "TP53" | gene == "NOTCH1", gene, "Other genes")) %>%
    mutate(genelab = factor(genelab, levels = c("TP53", "NOTCH1", "Other genes"))) %>%
    ggplot(aes(x = deltafit, fill = genelab)) +
    geom_histogram(bins = 25, position  = "stack", color = "white", size = 0.2) +
    ggtitle("Missense") +
    xlim(c(0, 0.1)) +
    xlab(expression(paste(Delta))) +
    ylab("Counts") +
    theme(legend.position = c(0.7, 0.8)) +
    scale_fill_jcolors(pal = "pal7",  name = "")

DFEnon <- DFresults %>%
    filter(mutationtype == "Nonsense") %>%
    mutate(genelab = ifelse(gene == "TP53" | gene == "NOTCH1", gene, "Other genes")) %>%
    mutate(genelab = factor(genelab, levels = c("TP53", "NOTCH1", "Other genes"))) %>%
    ggplot(aes(x = deltafit, group = genelab, fill = genelab, col = genelab)) +
    geom_histogram(bins = 25, position = "stack", color = "white", size = 0.2) +
    scale_fill_jcolors(pal = "pal7") +
    #scale_color_jcolors(pal = "pal7") +
    ggtitle("Nonsense") +
    xlim(c(0, 0.1)) +
    xlab(expression(paste(Delta))) +
    ylab("Counts") +
    theme(legend.position = "none")

DFE <- plot_grid(DFEmiss, DFEnon, labels = c("c", "d"))
message("Combine figures and save plot...")
figure3 <- plot_grid(geneplot, DFE, ncol = 1)
save_plot(args$figure, figure3, base_height = 7, base_width = 10)


message("Summary of per gene values...")
DFresults %>%
    filter(mutationtype == "Missense") %>%
    group_by(gene) %>%
    summarise(deltam = mean(deltafit),
              deltacilow = quantile(deltafit, 0.025),
              deltacihigh = quantile(deltafit, 0.975)) %>%
    arrange(desc(deltam)) %>% print()

DFresults %>%
    filter(mutationtype == "Nonsense") %>%
    group_by(gene) %>%
    summarise(deltam = mean(deltafit),
              deltacilow = quantile(deltafit, 0.025),
              deltacihigh = quantile(deltafit, 0.975)) %>%
    arrange(desc(deltam)) %>% print()

DFresults %>%
    ungroup() %>%
    filter(mutationtype == "Missense", gene %in% c("CREBBP", "FAT1", "NOTCH2", "PIK3CA")) %>%
    summarise(deltam = mean(deltafit),
              deltacilow = quantile(deltafit, 0.025),
              deltacihigh = quantile(deltafit, 0.975),
               deltamin = min(deltafit),
              deltamax = max(deltafit)) %>%
    arrange(desc(deltam)) %>% print()

DFresults %>%
    ungroup() %>%
    filter(mutationtype == "Nonsense", gene %in% c("CREBBP", "FAT1", "NOTCH2", "PIK3CA")) %>%
    summarise(deltam = mean(deltafit),
              deltacilow = quantile(deltafit, 0.025),
              deltacihigh = quantile(deltafit, 0.975),
               deltamin = min(deltafit),
              deltamax = max(deltafit)) %>%
    arrange(desc(deltam)) %>% print()

message("Plot all per patient fits")
textdfmiss <- dfmiss %>%
    mutate(p = paste0(patient, ", ", Age2),
    p = fct_reorder(p, Age2)) %>%
    mutate(label = paste("list(Delta == ", round(deltafit, 3), ",R^{2}==",round(rsq, 3) ,")")) %>%
    distinct(p, label)

gmiss <- dfmiss %>%
    mutate(p = paste0(patient, ", ", Age2),
    p = fct_reorder(p, Age2)) %>%
    ggplot(aes(x = A, y = dnds)) +
    geom_point(alpha = 0.5, col = "plum4") +
    geom_line(aes(y = dndsfit), alpha = 0.7, col = "firebrick") +
    geom_ribbon(aes(ymin = dndsfitlq, ymax = dndsfituq), fill = "firebrick", alpha = 0.3) +
    xlab("Clone area") +
    ylab("Interval dN/dS Missense") +
    geom_text(data = textdfmiss,
        mapping = aes(x = -Inf, y = -Inf, label = label),
        hjust   = -0.4,
        vjust   = -0.5, parse = T) +
    facet_wrap(~p, scales = "free", ncol = 3)

textdfnon <- dfnon %>%
    mutate(p = paste0(patient, ", ", Age2),
    p = fct_reorder(p, Age2)) %>%
    mutate(label = paste("list(Delta == ", round(deltafit, 3), ",R^{2}==",round(rsq, 3) ,")")) %>%
    distinct(p, label)

gnon <- dfnon %>%
    mutate(p = paste0(patient, ", ", Age2),
    p = fct_reorder(p, Age2)) %>%
    ggplot(aes(x = A, y = dnds)) +
    geom_point(alpha = 0.5, col = "plum4") +
    geom_line(aes(y = dndsfit), alpha = 0.7, col = "firebrick") +
    geom_ribbon(aes(ymin = dndsfitlq, ymax = dndsfituq), fill = "firebrick", alpha = 0.3) +
    xlab("Clone area") +
    ylab("Interval dN/dS Nonsense") +
    geom_text(data = textdfnon,
        mapping = aes(x = -Inf, y = -Inf, label = label),
        hjust   = -0.4,
        vjust   = -0.5, parse = T) +
    facet_wrap(~p, scales = "free", ncol = 3)

g <- plot_grid(gmiss, gnon, ncol = 2, labels = c("a", "b"))

save_plot(args$suppfigures[1], g, base_height = 8, base_width = 16)


message("summarise per patient global delta values")

nonsenseperpatient <- dfnon %>%
    distinct(patient, Age, Age2, deltafit, deltafitlq, deltafituq, rsq) %>%
    mutate(xlab = paste0(patient, " (", Age, ")"),
          barlab = paste0("R^{2}==",round(rsq, 3))) %>%
    filter(rsq > 0.6) %>%
    ggplot(aes(x = fct_reorder(xlab, Age2), y = deltafit)) +
    geom_bar(stat = "identity", fill = "darkslategrey", width = 0.5) +
    geom_text(aes(label = barlab), y = 0.04, parse = T) +
    geom_linerange(aes(ymin = deltafitlq, ymax = deltafituq)) +
    coord_flip() +
    ylim(c(0, 0.05)) +
    xlab("") +
    ylab(expression(Delta)) +
    ggtitle("Nonsense")

missenseperpatient <- dfmiss %>%
    distinct(patient, Age, Age2, deltafit, deltafitlq, deltafituq, rsq) %>%
    mutate(xlab = paste0(patient, " (", Age, ")"),
          barlab = paste0("R^{2}==",round(rsq, 3))) %>%
    filter(rsq > 0.6) %>%
    ggplot(aes(x = fct_reorder(xlab, Age2), y = deltafit)) +
    geom_bar(stat = "identity", fill = "darkolivegreen", width = 0.5) +
    geom_text(aes(label = barlab), y = 0.008, parse = T) +
    geom_linerange(aes(ymin = deltafitlq, ymax = deltafituq)) +
    coord_flip() +
    ylim(c(0, 0.01)) +
    xlab("") +
    ylab(expression(Delta)) +
    ggtitle("Missense")

patientplot <- plot_grid(missenseperpatient, nonsenseperpatient, align = T, labels = c("a", "b"))

message("summarise per patient global lambda values")
nonsenseperpatient <- dfnon %>%
    distinct(patient, Age, Age2, lambdarfit, lambdarfitlq, lambdarfituq, rsq) %>%
    mutate(xlab = paste0(patient, " (", Age, ")"),
          barlab = paste0("R^{2}==",round(rsq, 3))) %>%
    filter(rsq > 0.7) %>%
    ggplot(aes(x = fct_reorder(xlab, Age2), y = lambdarfit)) +
    geom_bar(stat = "identity", fill = "darkslategrey", width = 0.5) +
    geom_linerange(aes(ymin = lambdarfitlq, ymax = lambdarfituq)) +
    geom_text(aes(label = barlab), y = 12.0, parse = T) +
    coord_flip() +
    ylim(c(0, 15.0)) +
    xlab("") +
    ylab(expression(lambda)) +
    ggtitle("Nonsense")

missenseperpatient <- dfmiss %>%
    distinct(patient, Age, Age2, lambdarfit, lambdarfitlq, lambdarfituq, rsq) %>%
    mutate(xlab = paste0(patient, " (", Age, ")"),
          barlab = paste0("R^{2}==",round(rsq, 3))) %>%
    filter(rsq > 0.7) %>%
    ggplot(aes(x = fct_reorder(xlab, Age2), y = lambdarfit)) +
    geom_bar(stat = "identity", fill = "darkolivegreen", width = 0.5) +
    geom_linerange(aes(ymin = lambdarfitlq, ymax = lambdarfituq)) +
    geom_text(aes(label = barlab), y = 30.0, parse = T) +
    coord_flip() +
    ylim(c(0, 34.0)) +
    xlab("") +
    ylab(expression(lambda)) +
    ggtitle("Missense")

patientplot2 <- plot_grid(missenseperpatient, nonsenseperpatient, align = T, labels = c("c", "d"))

finalplot <- plot_grid(patientplot, patientplot2, ncol = 1)

save_plot(args$suppfigures[2], finalplot, base_height = 7, base_width = 10)


message("plot all per patient per gene fits")

textdf <- dfcombined %>%
    filter(rsq > 0.6, nmutations > 9, mutationtype == "Missense") %>%
    mutate(label = paste("list(Delta[fit] == ", round(deltafit, 3), ",R^{2}==",round(rsq, 3) ,")"),
          y1 = min(dnds), y2 = min(dndsfit),
          y = floor(min(y1, y2)) + 1) %>%
    distinct(y, deltafit, label)
gmiss <- dfcombined %>%
    filter(rsq > 0.6, nmutations > 7, mutationtype == "Missense") %>%
    ggplot(aes(x = A, y = dnds)) +
    geom_point(alpha = 0.9, size = 1) +
    geom_line(aes(y = dndsfit), alpha = 0.5, size = 1.0, col = "plum4") +
    geom_ribbon(aes(ymin = dndsfitlq, ymax = dndsfituq), fill = "plum4", alpha = 0.3) +
    xlab("Clone area") +
    ylab("Interval dN/dS - Missense") +
    facet_wrap(gene~patient, scale = "free_y")
    geom_text(data = textdf, aes(label = label, y = y), x = 3.0, size = 5, parse = TRUE)

textdf <- dfcombined %>%
    filter(rsq > 0.6, nmutations > 9, mutationtype == "Nonsense") %>%
    mutate(label = paste("list(Delta[fit] == ", round(deltafit, 3), ",R^{2}==",round(rsq, 3) ,")"),
          y1 = min(dnds), y2 = min(dndsfit),
          y = floor(min(y1, y2)) + 1) %>%
    distinct(y, deltafit, label)
gnon <- dfcombined %>%
    filter(rsq > 0.6, nmutations > 7, mutationtype == "Nonsense") %>%
    ggplot(aes(x = A, y = dnds)) +
    geom_point(alpha = 0.9, size = 1) +
    geom_line(aes(y = dndsfit), alpha = 0.5, size = 1.0, col = "darkseagreen4") +
    geom_ribbon(aes(ymin = dndsfitlq, ymax = dndsfituq), fill = "darkseagreen4", alpha = 0.3) +
    xlab("Clone area") +
    ylab("Interval dN/dS - Nonsense") +
    facet_wrap(gene~patient, scale = "free_y")
    geom_text(data = textdf, aes(label = label, y = y), x = 3.0, size = 5, parse = TRUE)

g <- plot_grid(gnon, gmiss, ncol = 1)

save_plot(args$suppfigures[3], g, base_height = 24, base_width = 12)

message("Delta and lambda summary")

glambda <- dfcombined %>%
    distinct(patient, Age, Age2, deltafit, deltafitlq, deltafituq, lambdarfit,
             lambdarfitlq, lambdarfituq, rsq, mutationtype, gene, nmutations) %>%
    mutate(plab = paste0(patient, " (", Age, ")"),
          barlab = paste0("R^{2}==",round(rsq, 3))) %>%
    mutate(plab = fct_reorder(plab, Age2)) %>%
    filter(rsq > 0.6, lambdarfit < 40.0, nmutations > 7) %>%
    ggplot(aes(x = gene, y = lambdarfit, fill = mutationtype, group = mutationtype)) +
    geom_bar(stat = "identity", position = position_dodge2(width = 0.9, preserve = "single"), width = 0.5) +
    geom_linerange(aes(ymin = lambdarfitlq, ymax = lambdarfituq),
                       position = position_dodge2(width = 0.5)) +
    coord_flip() +
    xlab("") +
    ylab(expression(lambda)) +
    facet_wrap(~plab, drop = TRUE, scales = "free") +
    scale_fill_jcolors(palette = "default") +
    theme(legend.position = "none") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))

gdelta <- dfcombined %>%
    distinct(patient, Age, Age2, deltafit, deltafitlq, deltafituq, lambdarfit,
             lambdarfitlq, lambdarfituq, rsq, mutationtype, gene, nmutations) %>%
    mutate(plab = paste0(patient, " (", Age, ")"),
          barlab = paste0("R^{2}==",round(rsq, 3))) %>%
    mutate(plab = fct_reorder(plab, Age2)) %>%
    filter(rsq > 0.6, nmutations > 7) %>%
    ggplot(aes(x = gene, y = deltafit, fill = mutationtype)) +
    geom_bar(stat = "identity", position = position_dodge2(width = 0.9, preserve = "single"), width = 0.5) +
    geom_linerange(aes(ymin = deltafitlq, ymax = deltafituq),
                  position = position_dodge2(width = 0.5)) +
    coord_flip() +
    xlab("") +
    ylab(expression(Delta)) +
    facet_wrap(~plab, drop = TRUE, scales = "free") +
    scale_fill_jcolors(palette = "default") +
    theme(legend.title = element_blank()) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))

save_plot(args$suppfigures[4], plot_grid(glambda, gdelta, ncol = 2), base_height = 10, base_width = 20)
