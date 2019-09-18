#load packages, you will need to install these if you don't have them already
library(brms)
library(cowplot)
library(tidyverse)
theme_set(theme_cowplot())
library(jcolors)

library(argparse)

parser <- ArgumentParser(description = "Fit brms models")
parser$add_argument('--oesophagusdata', type='character',
                    help="Mutations and VAF for oesophagus data")
parser$add_argument('--oesophagusmetadata', type='character',
                    help=" oesophagus meta data")
parser$add_argument('--output', type='character',
                    help="oesophagus meta data")
parser$add_argument('--threads', type='integer',
                    help="Number of threads", default = 1)
args <- parser$parse_args()


message("Read in data")
df <- read_csv(args$oesophagusdata)
donor <- readxl::read_xlsx(args$oesophagusmetadata, skip = 1) %>%
  dplyr::rename(donor = PD)
df <- left_join(df, donor)

message("Create data set")

dat <- df %>%
  mutate(area= 2 * sumvaf) %>%
  filter(impact != "Synonymous", aachange != ".")

message("Define brms paramters")
nchains <- args$threads
its <- 10000
formula <- bf(area ~ (1 + Age2|aachange))

genes <- c("NOTCH1", "TP53")
out <- list()

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
                filter(gene == g) %>%
                group_by(aachange) %>%
                mutate(n = n()) %>%
                ungroup() %>%
                filter(n > 4)
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
    out[[g]] <- brms_frechet
}

message("")
message("###########################################################")
message("Saving file")

saveRDS(out, args$output)
