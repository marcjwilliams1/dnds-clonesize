#load packages, you will need to install these if you don't have them already
library(cowplot)
library(tidyverse)
theme_set(theme_cowplot())
library(jcolors)
library(brms)

library(argparse)

parser <- ArgumentParser(description = "Fit brms models")
parser$add_argument('--oesophagusdata', type='character',
                    help="Mutations and VAF for oesophagus data")
parser$add_argument('--oesophagusmetadata', type='character',
                    help=" oesophagus meta data")
parser$add_argument('--output', type='character',
                    help=" oesophagus meta data")
args <- parser$parse_args()


message("Read in data")
df <- read_csv(args$oesophagusdata)
donor <- readxl::read_xlsx(args$oesophagusmetadata, skip = 1) %>%
  dplyr::rename(donor = PD)
df <- left_join(df, donor)

message("Create data set")

dat <- df %>%
  mutate(area= 2*sumvaf) %>%
  filter(str_detect(impact, "Synonymous|Missense|Nonsense")) %>%
  mutate(gene1 = gene) %>%
  mutate(gene1 = ifelse(impact == "Synonymous", "Synonymous", gene1))

message("Define brms paramters")
nchains <- 4
its <- 50
formula <- bf(area ~ (1 + Age2|gene1))

###########################################################
# Log normal distribution
###########################################################
message("")
message("###########################################################")
message("Brms fit with lognormal distribution")
brms_lognormal <- brm(formula,
                dat = dat,
                family = lognormal(),
                chains = nchains,
                cores = nchains,
                iter = its)
print(brms_lognormal)

###########################################################
# Frechet distribution
###########################################################
message("")
message("###########################################################")
message("Brms fit with frechet distribution")
brms_frechet <- brm(formula,
                dat = dat,
                family = frechet(),
                chains = nchains,
                cores = nchains,
                iter = its)
print(brms_frechet)

###########################################################
# Normal distribution
###########################################################
message("")
message("###########################################################")
message("Brms fit with normal distribution")
brms_normal <- brm(formula,
                dat = dat,
                family = gaussian(),
                chains = nchains,
                cores = nchains,
                iter = its)
print(brms_normal)

out <- list(normal = brms_normal,
            frecher = brms_frechet,
            lognormal = brms_lognormal)

message("")
message("###########################################################")
message("Saving file")
saveRDS(out, args$output)
