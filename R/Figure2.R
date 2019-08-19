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
parser$add_argument('--oesophagusfitmissense', type='character',
                    help="Fits for missense mutations oesophagus")
parser$add_argument('--oesophagusfitnonsense', type='character',
                    help="Fits for nonsense mutations oesophagus")
parser$add_argument('--skinfitmissense', type='character',
                    help="Fits for missense mutations oesophagus")
parser$add_argument('--skinfitnonsense', type='character',
                    help="Fits for nonsense mutations oesophagus")
parser$add_argument('--oesophagusfitmissensepergene', type='character',
                    help="Fits for missense mutations oesophagus per gene")
parser$add_argument('--oesophagusfitnonsensepergene', type='character',
                    help="Fits for nonsense mutations oesophagus per gene")
parser$add_argument('--oesophagusfitneutral', type='character',
                    help="Fits for neutral mutations oesophagus")
parser$add_argument('--stemcellpower', type='character',
                    help="Simulations for power to infer dynamics")
parser$add_argument('--stemcellexamplefit', type='character',
                    help="Example simulations fits")
args <- parser$parse_args()

message("Generating Figure 2...")
message("\t Reading in data...")

dfnon <- read_csv(args$oesophagusfitnonsense, col_types = cols())
dfmiss <- read_csv(args$oesophagusfitmissense, col_types = cols())

dfnon.gene <- read_csv(args$oesophagusfitnonsensepergene, col_types = cols())
dfmiss.gene <- read_csv(args$oesophagusfitnonsensepergene, col_types = cols())

dfneutral <- read_csv(args$oesophagusfitneutral, col_types = cols())

dfsim <- read.csv(args$stemcellexamplefit)
dfpower <- read_csv(args$stemcellpower, col_types = cols())


message("Summarise data...")

dfmiss %>%
    group_by(patient) %>%
    summarise(delta = first(deltafit), lambda = first(lambdarfit)) %>%
    summarise(deltam = mean(delta),
              deltacilow = quantile(delta, 0.025),
              deltacihigh = quantile(delta, 0.975),
              lambdam = mean(lambda),
              lambdacilow = quantile(lambda, 0.025),
              lambdacihigh = quantile(lambda, 0.975))

dfnon %>%
    group_by(patient) %>%
    summarise(delta = first(deltafit), lambda = first(lambdarfit)) %>%
    summarise(deltam = mean(delta),
              deltacilow = quantile(delta, 0.025),
              deltacihigh = quantile(delta, 0.975),
              lambdam = mean(lambda),
              lambdacilow = quantile(lambda, 0.025),
              lambdacihigh = quantile(lambda, 0.975),
              lambdamin = min(lambda),
              lambdamax = max(lambda))

message("Creating simulation figure...")
textdf <- dfsim %>%
    filter(deltatrue > 0.3, lambdartrue == 0.25) %>%
    group_by(deltatrue) %>%
    mutate(deltag = paste0(deltatrue)) %>%
     mutate(label = paste("list(Delta[input] == ", deltatrue,", Delta[fit] == ", round(deltafit, 3), ",R^{2}==",round(rsq, 3) ,")")) %>%
    distinct(deltatrue, deltag, label) %>%
    mutate(y = ifelse(deltatrue == 0.4, 6.5, 8.6))

gsim <- dfsim %>%
    filter(row_number() %% 5 == 0) %>%
    filter(deltatrue > 0.3, lambdartrue == 0.25) %>%
    filter(row_number() %% 5 == 0) %>%
    #filter(A < 15.0) %>%
    mutate(deltag = paste0(deltatrue)) %>%
    ggplot(aes(x = A)) +
    geom_point(aes(y = dnds,  group = deltag, fill = deltag, col = deltag), alpha = 0.9, size = 1) +
    geom_line(aes(y = dndsfit, group = deltag, fill = deltag, col = deltag), alpha = 0.3, size = 1) +
    geom_ribbon(aes(ymin = dndsfitlq, ymax = dndsfituq, group = deltag, fill = deltag, col = deltag), alpha = 0.1) +
    geom_text(data = textdf, aes(label = label, y = y, col = deltag), x = 15.0, size = 3, parse = TRUE) +
    xlab("Clone Area") +
    ylab("dN/dS") +
    ggtitle("Simulated cohort") +
    theme(legend.position = "none") +
    scale_color_jcolors(palette = "pal6")


library(Hmisc)

g1 <- dfpower %>%
    mutate(nsamples = factor(paste0(nsamples),
    levels = Cs(5, 8, 10, 25, 50, 100, 250))) %>%
    ggplot(aes(x = nsamples, y = delta)) +
    #geom_violin(fill = "steelblue4") +
    geom_sina(col = "steelblue4") +
    geom_boxplot(width = 0.3, alpha = 0.7, col = "steelblue4") +
    ylim(c(0, 1.0)) +
    geom_hline(yintercept = 0.25, lty = 2) +
    xlab("Number of mutations") +
    ylab(expression(paste("Inferred  ", Delta)))  +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))

g2 <- dfpower %>%
    mutate(nsamples = factor(paste0(nsamples),
    levels = Cs(5, 8, 10, 25, 50, 100, 250))) %>%
    ggplot(aes(x = nsamples, y = rsq)) +
    geom_sina(col = "steelblue4") +
    geom_boxplot(width = 0.3, alpha = 0.7, col = "steelblue4") +
    ylim(c(0, 1.0)) +
    xlab("Number of mutations") +
    ylab(expression(R^2)) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))

figure2top <- plot_grid(gsim, g1, g2, ncol = 3, labels = c("a", "b", "c"))

textdf <- dfnon.gene %>%
    filter(	patient == "PD30274", gene == "FAT1") %>%
    mutate(label = paste("list(Delta[fit] == ", round(deltafit, 3), ",R^{2}==",round(rsq, 3) ,")"),
          y1 = min(dnds), y2 = min(dndsfit),
          y = floor(min(y1, y2)) + 1) %>%
    distinct(y, deltafit, label)
g2 <- dfnon.gene %>%
    filter(patient == "PD30274", gene == "FAT1") %>%
    ggplot(aes(x = A, y = dnds)) +
    geom_point(alpha = 0.9, col = "black", size = 1) +
    geom_line(aes(y = dndsfit), alpha = 0.5, size = 1.0, col = "plum4") +
    geom_ribbon(aes(ymin = dndsfitlq, ymax = dndsfituq), fill = "plum4", alpha = 0.3) +
    xlab("Clone area") +
    ylab("Interval dN/dS - Nonsense") +
    ggtitle("PD30274 - FAT1") +
    geom_text(data = textdf, aes(label = label, y = y), x = 2.5, size = 5, parse = TRUE)

textdf <- dfneutral %>%
    mutate(label = paste("list(Delta[fit] == ", round(deltafit, 3),")"),
          y1 = min(dnds), y2 = min(dndsfit),
          y = floor(min(y1, y2)) + 0) %>%
    distinct(y, deltafit, label)
g1 <- dfneutral %>%
    ggplot(aes(x = A, y = dnds)) +
    geom_point(alpha = 0.9, col = "black", size = 1) +
    geom_line(aes(y = dndsfit), alpha = 0.5, size = 1.0, col = "cornsilk3") +
    geom_ribbon(aes(ymin = dndsfitlq, ymax = dndsfituq), fill = "cornsilk3", alpha = 0.3) +
    xlab("Clone area") +
    ylab("Interval dN/dS") +
    ggtitle("All Neutral genes") +
    geom_text(data = textdf, aes(label = label, y = y), x = 2.5, size = 5, parse = TRUE) +
    ylim(c(0, 2))

textdf <- dfmiss.gene %>%
    filter(	patient == "PD30988", gene == "NOTCH1") %>%
    mutate(label = paste("list(Delta[fit] == ", round(deltafit, 3), ",R^{2}==",round(rsq, 3) ,")"),
          y1 = min(dnds), y2 = min(dndsfit),
          y = floor(min(y1, y2)) + 1) %>%
    distinct(y, deltafit, label)
g2 <- dfmiss.gene %>%
    filter(patient == "PD30988", gene == "NOTCH1") %>%
    ggplot(aes(x = A, y = dnds)) +
    geom_point(alpha = 0.9, col = "black", size = 1) +
    geom_line(aes(y = dndsfit), alpha = 0.5, size = 1.0, col = "darkseagreen4") +
    geom_ribbon(aes(ymin = dndsfitlq, ymax = dndsfituq), fill = "darkseagreen4", alpha = 0.3) +
    xlab("Clone area") +
    ylab("Interval dN/dS - Missense") +
    ggtitle("PD31182 - NOTCH1") +
    geom_text(data = textdf, aes(label = label, y = y), x = 2.5, size = 5, parse = TRUE)


textdf <- dfmiss.gene %>%
    filter(	patient == "PD30273", gene == "TP53") %>%
    mutate(label = paste("list(Delta[fit] == ", round(deltafit, 3), ",R^{2}==",round(rsq, 3) ,")"),
          y1 = min(dnds), y2 = min(dndsfit),
          y = floor(min(y1, y2)) + 1) %>%
    distinct(y, deltafit, label)
g3 <- dfmiss.gene %>%
    filter(patient == "PD30273", gene == "TP53") %>%
    ggplot(aes(x = A, y = dnds)) +
    geom_point(alpha = 0.9, col = "black", size = 1) +
    geom_line(aes(y = dndsfit), alpha = 0.5, size = 1.0, col = "darkseagreen4") +
    geom_ribbon(aes(ymin = dndsfitlq, ymax = dndsfituq), fill = "darkseagreen4", alpha = 0.3) +
    xlab("Clone area") +
    ylab("Interval dN/dS - Missense") +
    ggtitle("PD30273 - TP53") +
    geom_text(data = textdf, aes(label = label, y = y), x = 2.5, size = 5, parse = TRUE)

textdf <- dfnon.gene %>%
    filter(	patient == "PD31182", gene == "NOTCH1") %>%
    mutate(label = paste("list(Delta[fit] == ", round(deltafit, 3), ",R^{2}==",round(rsq, 3) ,")"),
          y1 = min(dnds), y2 = min(dndsfit),
          y = floor(min(y1, y2)) + 1) %>%
    distinct(y, deltafit, label)
g4 <- dfnon.gene %>%
    filter(patient == "PD31182", gene == "NOTCH1") %>%
    ggplot(aes(x = A, y = dnds)) +
    geom_point(alpha = 0.9, col = "black", size = 1) +
    geom_line(aes(y = dndsfit), alpha = 0.5, size = 1.0, col = "plum4") +
    geom_ribbon(aes(ymin = dndsfitlq, ymax = dndsfituq), fill = "plum4", alpha = 0.3) +
    xlab("Clone area") +
    ylab("Interval dN/dS - Nonsense") +
    ggtitle("PD31182 - NOTCH1") +
    geom_text(data = textdf, aes(label = label, y = y), x = 2.5, size = 5, parse = TRUE)

figure2bottom <- plot_grid(g1, g2, g3, g4, labels = c("d", "e", "f", "g"), ncol = 4)

message("Combine plots together...")
figure2 <- plot_grid(figure2top, figure2bottom, ncol = 1)

message("Save plot to file")
save_plot(args$figure, figure2, base_height = 7, base_width = 12)
