library(tidyverse)
library(jcolors)
library(argparse)
library(tidybayes)
library(bayesplot)
library(modelr)
library(cowplot)

parser <- ArgumentParser(description = "Plot simulation fits output")
parser$add_argument('--simulationdata', type='character',
                    help="Simulation data for hitchikers")
parser$add_argument('--suppfigures', type='character',
                    help="Output figure files", nargs = "+")
parser$add_argument('--rho', type='double',
                    help="Progenitor density", default = 5000.0)
parser$add_argument('--binsize', type='double',
                    help="Binsize for fitting", default = 0.002)
parser$add_argument('--oesophagusdata', type='character',
                    help="Mutations and VAF for oesophagus data")
parser$add_argument('--oesophagusmetadata', type='character',
                    help=" oesophagus meta data")
args <- parser$parse_args()

args <- list(simulationdata = "~/Documents/apocrita/BCInew/marc/dnds/dnds-clonesize/results/simulations/clonesize_hitchikers.csv",
             datamodelfits = "~/Documents/apocrita/BCInew/marc/dnds/dnds-clonesize/results/dataforfigures/data-clonesizefit-models.Rdata",
             oesophaagusdata = "~/Documents/apocrita/BCInew/marc/dnds/dnds-clonesize/data/oesophagus/esophagus.csv",
             oesophaagusmetadata = "~/Documents/apocrita/BCInew/marc/dnds/dnds-clonesize/data/oesophagus/patient_info.xlsx",
             rho = 5000)

dfsims <- read_csv(args$simulationdata)

dfdata <- read_csv(args$oesophagusdata)
dfdata <- read_csv("~/Documents/apocrita/BCInew/marc/dnds/dnds-clonesize/data/oesophagus/esophagus.csv")
donor <- readxl::read_xlsx("~/Documents/apocrita/BCInew/marc/dnds/dnds-clonesize/data/oesophagus/patient_info.xlsx", skip = 1) %>%
  dplyr::rename(donor = PD)
dfdata <- left_join(dfdata, donor)


dfsims %>%
  #filter(muttype == "NonSyn") %>%
  filter(muNS == 0.001) %>%
  mutate(f = f / 5000) %>%
  ggplot(aes(x = f)) +
  #geom_bar(stat = "identity") +
  geom_histogram(bins = 100) +
  scale_x_log10() +
  facet_grid(muttype ~ tend, scales = "free_y")

