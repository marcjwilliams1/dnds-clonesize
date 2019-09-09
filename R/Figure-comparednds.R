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
                    help="Output figure files")
parser$add_argument('--suppfigures', type='character',
                    help="Outpute figure files", nargs = "+", default = NULL)
parser$add_argument('--alldndscv', type='character',
                    help="dndscv files for all genes")
parser$add_argument('--SSB', type='character',
                    help="SSB file for all genes")
parser$add_argument('--dndscvmissense', type='character',
                    help="dndscv missense")
parser$add_argument('--dndscvnonsense', type='character',
                    help="dndscv nonsense")
parser$add_argument('--dndscvnonsensepergene', type='character',
                    help="dndscv nonsense per gene")
parser$add_argument('--dndscvmissensepergene', type='character',
                    help="dndscv missense per gene")
args <- parser$parse_args()

message("Read in files....")
dndscvall <- read_csv(args$alldndscv) %>%
    mutate(mutation_class = "All")
dndscvnonsense <- read_csv(args$dndscvnonsense) %>%
    mutate(mutation_class = "Nonsense")
dndscvmissense <- read_csv(args$dndscvmissense) %>%
    mutate(mutation_class = "Missense")
dndscv <- bind_rows(dndscvall, dndsvnonsense, dndscvmissense)
SSB <- read_csv(args$SSB)

df <- left_join(dndscv, SSB, by = c("patient", "Age", "Age2", "mutation_class"))


summarydf <- df %>%
  distinct(deltafit.x, deltafit.y, patient, Age2, mutation_class,
           nmutations.x, nmutations.y, rsq.x, rsq.y, gene)

g <- summarydf %>%
  filter(gene == "NOTCH1") %>%
  ggplot(aes(x = deltafit.x, y = deltafit.y)) +
  geom_point() +
  facet_wrap(~mutation_class, scales = "free") +
  geom_abline(lty = 2, col = "firebrick") +
  xlab(expression(Delta*" (dndscv)")) +
  ylab(expression(Delta*" (SSB)"))
save_plot(args$figure, g, base_height = 3.5, base_width = 10)



message("Plot all per patient fits for SSB")
textdfmiss <- SSB %>%
  filter(mutation_class == "Missense", gene == "All") %>%
  mutate(p = paste0(patient, ", ", Age2),
         p = fct_reorder(p, Age2)) %>%
  mutate(label = paste("list(Delta == ", round(deltafit, 3), ",R^{2}==",round(rsq, 3) ,")")) %>%
  distinct(p, label)

gmiss <- SSB %>%
  filter(mutation_class == "Missense", gene == "All") %>%
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

textdfnon <- SSB %>%
  filter(mutation_class == "Nonsense", gene == "All") %>%
  mutate(p = paste0(patient, ", ", Age2),
         p = fct_reorder(p, Age2)) %>%
  mutate(label = paste("list(Delta == ", round(deltafit, 3), ",R^{2}==",round(rsq, 3) ,")")) %>%
  distinct(p, label)

gnon <- SSB %>%
  filter(mutation_class == "Nonsense", gene == "All") %>%
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

save_plot("Figures/SSBfits.pdf", g, base_height = 8, base_width = 16)
