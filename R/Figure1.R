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
args <- parser$parse_args()

message("Generating Figure 1...")
f <- seq(0.0, 1.0, 0.01)
dnds <- 2 * f + 1
dfpos <- data_frame(f = f, dnds = dnds, s = "Positive selection")

dnds <- -2 * f + 1
dfneg <- data_frame(f = f, dnds = dnds, s = "Negative selection")

df <- bind_rows(dfpos, dfneg)

gdnds <- df %>%
    ggplot(aes(x = f, y = dnds, group = s, col = s)) +
    geom_line(size = 2) +
    theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank()) +
    geom_hline(yintercept = 1.0, lty = 2) +
    scale_y_continuous(breaks = c(1.0)) +
    xlab("Frequency") +
    ylab("dN/dS") +
    scale_color_manual(values = c( "dodgerblue4", "firebrick4")) +
    theme(legend.position=c(0.1, 0.8)) +
    theme(legend.title=element_blank())

dnds <- 0.5 * f + 1
dfpos1 <- data_frame(f = f, dnds = dnds, s = "Low s")
dnds <- 2.0 * f + 1
dfpos2 <- data_frame(f = f, dnds = dnds, s = "Medium s")
dnds <- 3.0 * f + 1
dfpos3 <- data_frame(f = f, dnds = dnds, s = "High s")

df <- bind_rows(dfpos1, dfpos2, dfpos3) %>%
    mutate(s = factor(s, levels = c("Low s", "Medium s", "High s")))

gdifferents <- df %>%
    ggplot(aes(x = f, y = dnds, group = s, col = s)) +
    geom_line(size = 2) +
    geom_text(data = filter(df, f == last(f)),
              aes(label = s, x = f - 0.05, y = 0.85 * dnds, color = s)) +
    theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank()) +
    geom_hline(yintercept = 1.0, lty = 2) +
    scale_y_continuous(breaks = c(1.0), limits = c(-1, 4)) +
    xlab("Frequency") +
    ylab("dN/dS") +
    scale_color_manual(values = c( "firebrick1", "firebrick3", "firebrick4")) +
    theme(legend.position="none")

#create figures leaving space to insert cartoons
figure1 <- plot_grid(NULL, gdnds, gdifferents, ncol = 3, rel_widths = c(1, 0.7, 0.7),
                         labels = c("a", "b", "c"))
save_plot(args$figure, figure1, base_height = 3.5, base_width = 15)
