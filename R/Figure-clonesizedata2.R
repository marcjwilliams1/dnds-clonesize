library(brms)
library(tidyverse)
library(jcolors)
library(argparse)
library(tidybayes)
library(bayesplot)
library(modelr)
library(cowplot)


parser <- ArgumentParser(description = "Plot simulation fits output")
parser$add_argument('--datafits', type='character',
                    help="Data fits")
parser$add_argument('--figure', type='character',
                    help="Output figure files", default = NULL)
parser$add_argument('--oesophagusmetadata', type='character',
                    help=" oesophagus meta data")
parser$add_argument('--suppfigures', type='character',
                    help="Output figure files", nargs = "+")
parser$add_argument('--rho', type='double',
                    help="Progenitor density", default = 5000.0)
parser$add_argument('--oesophagusfitmissense', type='character',
                    help="Fits for missense mutations oesophagus")
args <- parser$parse_args()



#args <- list(datafits = "~/Documents/apocrita/BCInew/marc/dnds/dnds-clonesize/results/dataforfigures/brmsfit.Rdata",
#             datamodelfits = "~/Documents/apocrita/BCInew/marc/dnds/dnds-clonesize/results/dataforfigures/data-clonesizefit-models.Rdata",
#             oesophaagusmetadata = "~/Documents/apocrita/BCInew/marc/dnds/dnds-clonesize/data/oesophagus/patient_info.xlsx",
#             rho = 5000)

modelfits <- readRDS(args$datafits)


message("Perform some model selection")
waicdf <- data.frame(waic = c(modelfits$frechet$waic$estimates[3],
                              modelfits$lognormal$waic$estimates[3],
                              modelfits$normal$waic$estimates[3]),
                     waicse = c(modelfits$frechet$waic$estimates[3,2],
                                modelfits$lognormal$waic$estimates[3,2],
                                modelfits$normal$waic$estimates[3,2]),
                     R2 = c(mean(modelfits$frechet$R2),
                            mean(modelfits$lognormal$R2),
                            mean(modelfits$normal$R2)),
                     R2min = c(quantile(modelfits$frechet$R2, 0.025),
                               quantile(modelfits$lognormal$R2, 0.025),
                               quantile(modelfits$normal$R2, 0.025)),
                     R2max = c(quantile(modelfits$frechet$R2, 0.975),
                               quantile(modelfits$lognormal$R2, 0.975),
                               quantile(modelfits$normal$R2, 0.975)),
                     model = c("Frechet", "Log Normal", "Normal"),
                     ord = c(1,2,3))

gwaic <- waicdf %>%
  filter(model != "Frechet") %>%
  mutate(ymin = waic - 2*waicse, ymax = waic + 2*waicse) %>%
  ggplot(aes(x = fct_reorder(model, ord), y = waic, ymin = ymin, ymax = ymax)) +
  geom_point() +
  geom_linerange() +
  xlab("") +
  ylab("WAIC") +
  theme_cowplot() +
  coord_flip()

(ppcheckfrech <- plot(pp_check(modelfits$frechet, nsamples = 50)) + 
    scale_x_log10(limits = c(1e-3, 1e1)) +
     theme_cowplot() +
    xlab("Area") +
    ylab("Density"))

(ppchecknormal <- plot(pp_check(modelfits$normal, nsamples = 50)) + 
    scale_x_log10(limits = c(1e-3, 1e1)) +
    theme_cowplot() +
    xlab("Area") +
    ylab("Density"))

(ppchecklognormal <- plot(pp_check(modelfits$lognormal, nsamples = 50)) + 
    scale_x_log10(limits = c(1e-3, 1e1)) +
    theme_cowplot() +
    xlab("Area") +
    ylab("Density"))

ppcheck <- ppchecklognormal

(g <- plot_grid(ppcheckfrech + ggtitle("Frechet"),
          ppchecklognormal + ggtitle("Log normal"),
          ppchecknormal + ggtitle("Normal"), ncol = 3))
save_plot(args$suppfigure[2], g, base_aspect_ratio = 4.0)


frechetCoef <- modelfits$lognormal %>%
  spread_draws(r_gene1[gene, var], b_Age2) %>%
  filter(var == "Age2") %>%
  median_qi(coef = r_gene1 + b_Age2) %>%
  arrange(desc(coef))

(gsummary <- frechetCoef %>%
  ggplot(aes(x = fct_reorder(gene, coef, .desc = TRUE), y = coef, 
             xmin = .lower, xmax = .upper)) +
  geom_pointinterval() +
  xlab("") +
  ylab("Regression coefficient") +
  geom_hline(yintercept = 0.0, lty = 2, col = "firebrick") +
  theme_cowplot() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)))


dfg <- read_csv(args$oesophagusfitmissense)
DFresults <- dfg %>%
  group_by(gene, deltafit, nmutations, deltafitlq, deltafituq) %>%
  summarise(rsq = first(rsq)) %>%
  filter(rsq > 0.6, nmutations > 7) %>%
  group_by(gene) %>% mutate(n = n()) %>%
  summarise(deltafit = mean(deltafit))

(gcompare <- left_join(DFresults, frechetCoef) %>%
  ggplot(aes(x = deltafit, y = coef)) +
  geom_text(aes(label = gene)) +
  scale_y_log10() +
  scale_x_log10() +
  geom_smooth(method = "lm", se = F) +
  xlab(expression(Delta~" fit - dN/dS")) +
  ylab("Regression coefficient") +
  theme_cowplot())

message("Generate final figure")
g <- plot_grid(gwaic, ppcheck + theme(legend.position = c(0.7, 0.9)), 
               gcompare, ncol = 3, labels = c("a", "b", "c"), align = "h")
(gall <- plot_grid(g, gsummary, ncol = 1, labels = c("", "d")))
save_plot(args$suppfigures[1], gall, base_height = 10, base_width = 20)




