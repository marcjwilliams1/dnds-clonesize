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
parser$add_argument('--binsizesims', type='character',
                    help="Bin size simulations")
args <- parser$parse_args()

df <- read_csv(args$binsizesims)
