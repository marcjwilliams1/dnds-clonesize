library(tidyverse)
library(jcolors)
library(argparse)
library(cowplot)

parser <- ArgumentParser(description = "Plot simulation fits output")
parser$add_argument('--suppfigure', type='character',
                    help="Output figure files", nargs = "+")
parser$add_argument('--exp', type='character',
                    help="Progenitor density")
parser$add_argument('--beta', type='character',
                    help="Binsize for fitting")
args <- parser$parse_args()

exp <- read_csv(args$exp)
beta <- read_csv(args$beta)


g1 <- exp %>%
  filter(rsq > 0.6, deltafit > 0.0, deltatrue == 0.05) %>%
  #sample_n(500) %>%
  ggplot(aes(x = deltafit)) +
  geom_histogram(bins = 50, alpha = 0.9, fill = "black") +
  geom_vline(xintercept = 0.05, lty = 2, col = "firebrick", size = 2) +
  xlab(expression("Inferred"~Delta)) +
  theme_cowplot() +
  ylab("Counts") +
  ggtitle(expression("True "~Delta~"= 0.05 (Mean of exponential)")) +
  xlim(c(0, 1.0))

g2 <- exp %>%
  #sample_n(500) %>%
  filter(rsq > 0.6, deltafit > 0.0, deltatrue == 0.1) %>%
  ggplot(aes(x = deltafit)) +
  geom_histogram(bins = 50, alpha = 0.9, fill = "black") +
  geom_vline(xintercept = 0.1, lty = 2, col = "firebrick", size = 2) +
  xlab(expression("Inferred"~Delta)) +
  theme_cowplot() +
  ylab("Counts") +
  ggtitle(expression("True "~Delta~"= 0.1 (Mean of exponential)")) +
  xlim(c(0, 1.0))

g2a <- exp %>%
  #sample_n(500) %>%
  filter(rsq > 0.6, deltafit > 0.0, deltatrue == 0.2) %>%
  ggplot(aes(x = deltafit)) +
  geom_histogram(bins = 50, alpha = 0.9, fill = "black") +
  geom_vline(xintercept = 0.2, lty = 2, col = "firebrick", size = 2) +
  xlab(expression("Inferred"~Delta)) +
  theme_cowplot() +
  ylab("Counts") +
  ggtitle(expression("True "~Delta~"= 0.2 (Mean of exponential)")) +
  xlim(c(0, 1.0))

g3 <- beta %>%
  filter(rsq > 0.6, deltafit > 0.0, deltatrue == 0.05) %>%
  ggplot(aes(x = deltafit)) +
  geom_histogram(bins = 50, alpha = 0.9, fill ="dodgerblue4") +
  geom_vline(xintercept = 0.05, lty = 2, col = "firebrick", size = 2) +
  xlab(expression("Inferred"~Delta)) +
  theme_cowplot() +
  ylab("Counts") +
  ggtitle(expression("True "~Delta~"= 0.05")) +
  xlim(c(0, 1.0))

g4 <- beta %>%
  filter(rsq > 0.6, deltafit > 0.0, deltatrue == 0.1) %>%
  ggplot(aes(x = deltafit)) +
  geom_histogram(bins = 50, alpha = 0.9, fill = "dodgerblue4") +
  geom_vline(xintercept = 0.1, lty = 2, col = "firebrick", size= 2) +
  xlab(expression("Inferred"~Delta)) +
  theme_cowplot() +
  ylab("Counts") +
  ggtitle(expression("True "~Delta~"= 0.1")) +
  xlim(c(0, 1.0))


g <- plot_grid(g1, g2, g2a, ncol = 3)
save_plot(args$suppfigure, g, base_width = 13, base_height = 4)

exp %>%
  filter(rsq > 0.6, deltafit > 0.0) %>%
  #sample_n(500) %>%
  group_by(deltatrue) %>%
  summarise(deltaMedian = median(deltafit), deltaMean = mean(deltafit),  deltafitl = quantile(deltafit, 0.025),
            deltafitu = quantile(deltafit, 0.975),
            n = n()) %>%
  print()


beta %>%
  filter(rsq > 0.6, deltafit > 0.0) %>%
  group_by(deltatrue) %>%
  summarise(deltaMedian = median(deltafit), deltaMean = mean(deltafit),  deltafitl = quantile(deltafit, 0.025),
	    deltafitu = quantile(deltafit, 0.975),
            n = n()) %>%
  print()

