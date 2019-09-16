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
  filter(impact != "Synonymous", aachange != ".")

message("Define brms paramters")
nchains <- 4
its <- 5000
formula <- bf(area ~ (1 + Age2|aachange))

genes <- c("NOTCH1", "TP53")
output <- list()

###########################################################
# Frechet distribution
###########################################################
message("")
message("###########################################################")
message("Brms fit with frechet distribution")
for (g in genes){
    message("###########################################################")
    message(paste0("Gene: ", g))
    dattemp <- dat %>%
                filter(gene == g)
    brms_frechet <- brm(formula,
                    dat = dat,
                    family = frechet(),
                    chains = nchains,
                    cores = nchains,
    		        control = list(adapt_delta = 0.8),
                    iter = its)
    brms_frechet <- add_criterion(brms_frechet,
                                c("loo", "waic", "R2"))
    print(brms_frechet)
    print(bayes_R2(brms_frechet))
    output[[g]] <- brms_frechet
    message("")
}

message("")
message("###########################################################")
message("Saving file")

saveRDS(output, args$output)
