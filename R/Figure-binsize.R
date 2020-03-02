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
                    help="Output figure files", default = NULL)
parser$add_argument('--suppfigures', type='character',
                    help="Outpute figure files", nargs = "+", default = NULL)
parser$add_argument('--binsizesims', type='character',
                    help="Bin size simulations")
args <- parser$parse_args()
dfsim <- read_csv(args$binsizesims) %>%
  filter(stepsize < 0.1)

del <- 0.4

gsummary <- dfsim %>%
  filter(deltatrue == del) %>%
  distinct(stepsize, deltafit, deltafitlq, deltafituq, deltatrue) %>%
  ggplot(aes(x = factor(stepsize), y = deltafit, ymin = deltafitlq, ymax = deltafituq)) +
  geom_point() +
  geom_linerange() +
  ylim(c(0.4, 0.41)) +
  xlab("Bin size") +
  ylab(expression("Inferred " *Delta)) +
  ggtitle(expression("Input "*Delta*" = 0.4"))


textdf <- dfsim %>%
  filter(deltatrue == del) %>%
  group_by(deltatrue, stepsize) %>%
  mutate(deltag = paste0(deltatrue)) %>%
  mutate(label = paste("list(Delta[input] == ", deltatrue,", Delta[fit] == ", round(deltafit, 3), ",R^{2}==",round(rsq, 3) ,")")) %>%
  distinct(deltatrue, deltag, label) %>%
  mutate(y = 6.0)

gsim <- dfsim %>%
  filter(row_number() %% 5 == 0) %>%
  filter(deltatrue == del) %>%
  filter(row_number() %% 5 == 0) %>%
  #filter(A < 15.0) %>%
  mutate(deltag = paste0(deltatrue)) %>%
  ggplot(aes(x = A)) +
  geom_point(aes(y = dnds,  group = deltag, fill = deltag, col = deltag), alpha = 0.9, size = 1) +
  geom_line(aes(y = dndsfit, group = deltag, fill = deltag, col = deltag), alpha = 0.3, size = 1) +
  geom_ribbon(aes(ymin = dndsfitlq, ymax = dndsfituq, group = deltag, fill = deltag, col = deltag), alpha = 0.1) +
  geom_text(data = textdf, aes(label = label, y = y, col = deltag), x = 2.5, size = 3, parse = TRUE) +
  xlab("Clone Area") +
  ylab("dN/dS") +
  theme(legend.position = "none") +
  scale_color_jcolors(palette = "pal6") +
  facet_wrap(~stepsize) +
  xlim(c(0, 5))


gout <- plot_grid(gsim, gsummary, labels = c("a", "b"))
save_plot(args$suppfigures[1], gout, base_height = 5, base_width = 13)
