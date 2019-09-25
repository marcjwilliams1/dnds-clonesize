library(brms)
library(tidyverse)
library(jcolors)
library(argparse)

parser <- ArgumentParser(description = "Fit brms models")
parser$add_argument('--oesophagusdata', type='character',
                    help="Mutations and VAF for oesophagus data")
parser$add_argument('--oesophagusmetadata', type='character',
                    help=" oesophagus meta data")
parser$add_argument('--output', type='character',
                    help=" oesophagus meta data")
parser$add_argument('--threads', type='integer',
                    help="Number of threads", default = 1)
parser$add_argument('--rho', type='double',
                    help="Progenitor density", default = 5000.0)
parser$add_argument('--binsize', type='double',
                    help="Bin size for fitting", default = 0.005)
parser$add_argument('--its', type='integer',
                    help="Progenitor density", default = 5000)
parser$add_argument('--quantile', type='double',
                    help="quantile with which to filter data", default = 0.9)
parser$add_argument('--minvaf', type='double',
                    help="Progenitor density", default = 0.008)
args <- parser$parse_args()


message("Read in data")
df <- read_csv(args$oesophagusdata)
donor <- readxl::read_xlsx(args$oesophagusmetadata, skip = 1) %>%
  dplyr::rename(donor = PD)
df <- left_join(df, donor)

df <- df %>%
    filter(sumvaf > args$minvaf)

# Functions to calculate ground truth

midcut<-function(x,from,to,by){
  ## cut the data into bins...
  x=cut(x,seq(from,to,by),include.lowest=T, right = F)
  ## make a named vector of the midpoints, names=binnames
  vec=seq(from+by/2,to-by/2,by)
  #vec=seq(from,to,by)
  names(vec)=levels(x)
  ## use the vector to map the names of the bins to the midpoint values
  unname(vec[x])
}

message("")
message("###########################################################")
message("Fit model for age")

mydat <- df %>%
  filter(impact != "Synonymous") %>%
  group_by(Age) %>%
  mutate(maxvaf = quantile(sumvaf, args$quantile)) %>%
  ungroup() %>%
  mutate(nidx = midcut(sumvaf, args$minvaf, 2, args$binsize)) %>%
  filter(!is.na(nidx)) %>%
  group_by(Age, nidx, maxvaf) %>%
  summarise(C = n()) %>%
  ungroup() %>%
  rename(n = nidx) %>%
  complete(Age, nesting(n), fill = list(C = 0)) %>%
  #filter(C > 1) %>%
  filter(n < maxvaf)

prior1 <- prior(normal(5, 2), nlpar = "A", lb = 0.0001) +
  prior(normal(0, 5), nlpar = "B")
nchains <- args$threads

message("Fit model...")
fitage <- brm(bf(C ~ (A / n) * exp(-n / exp(B)),
               A ~ 1 + (1|Age),
               B ~ 1 + (1|Age),
               nl = TRUE),
            data = mydat,
            prior = prior1,
            family = gaussian,
            control = list(adapt_delta = 0.9),
            chains = nchains,
            iter = args$its)

fitage<- add_criterion(fitage, c("loo", "waic", "R2"))
print(fitage)
print(bayes_R2(fitage))



message("")
message("###########################################################")
message("Fit model for age synonymous")

mydat <- df %>%
  filter(impact == "Synonymous") %>%
  group_by(Age) %>%
  mutate(maxvaf = quantile(sumvaf, args$quantile)) %>%
  ungroup() %>%
  mutate(nidx = midcut(sumvaf, args$minvaf, 2, args$binsize)) %>%
  filter(!is.na(nidx)) %>%
  group_by(Age, nidx, maxvaf) %>%
  summarise(C = n()) %>%
  ungroup() %>%
  rename(n = nidx) %>%
  complete(Age, nesting(n), fill = list(C = 0)) %>%
  #filter(C > 1) %>%
  filter(n < maxvaf)

prior1 <- prior(normal(5, 2), nlpar = "A", lb = 0.0001) +
  prior(normal(0, 5), nlpar = "B")
nchains <- args$threads

message("Fit model...")
fitagesynon <- brm(bf(C ~ (A / n) * exp(-n / exp(B)),
               A ~ 1 + (1|Age),
               B ~ 1 + (1|Age),
               nl = TRUE),
            data = mydat,
            prior = prior1,
            family = gaussian,
            control = list(adapt_delta = 0.9),
            chains = nchains,
            iter = args$its)

fitagesynon <- add_criterion(fitagesynon, c("loo", "waic", "R2"))
print(fitagesynon)
print(bayes_R2(fitagesynon))



message("")
message("###########################################################")
message("Fit model pergene")

mydat <- df %>%
  filter(impact != "Synonymous") %>%
  mutate(gene = paste0(gene, "-", Age2)) %>%
  group_by(gene) %>%
  mutate(nmuts = n(),
        maxvaf = quantile(sumvaf, args$quantile)) %>%
  ungroup() %>%
  filter(nmuts > 9) %>%
  mutate(nidx = midcut(sumvaf, args$minvaf, 2, args$binsize)) %>%
  filter(!is.na(nidx)) %>%
  group_by(gene, nidx, maxvaf) %>%
  summarise(C = n()) %>%
  ungroup() %>%
  rename(n = nidx) %>%
  complete(gene, nesting(n), fill = list(C = 0)) %>%
  #filter(C > 0) %>%
  filter(!is.na(n)) %>%
  filter(n < maxvaf)

prior1 <- prior(normal(5, 2), nlpar = "A", lb = 0.0001) +
  prior(normal(0, 5), nlpar = "B")
nchains <- args$threads

message("Fit model...")
fitgene <- brm(bf(C ~ (A / n) * exp(-n / exp(B)),
               A ~ 1 + (1|gene),
               B ~ 1 + (1|gene),
               nl = TRUE),
            data = mydat,
            prior = prior1,
            family = gaussian,
            control = list(adapt_delta = 0.9),
            chains = nchains,
            cores = nchains,
            iter = args$its)

fitgene<- add_criterion(fitgene, c("loo", "waic", "R2"))
print(fitgene)
print(bayes_R2(fitgene))


message("")
message("###########################################################")
message("Fit model pergene and per age")

mydat <- df %>%
  filter(impact != "Synonymous") %>%
  filter(str_detect(gene, "NOTCH1|TP53")) %>%
  #group_by(gene, Age2) %>%
  #mutate(nmuts = n(),          maxvaf = quantile(sumvaf, args$quantile)) %>%
  #ungroup() %>%
  #filter(nmuts > 9) %>%
  mutate(nidx = midcut(sumvaf, args$minvaf, 2, args$binsize)) %>%
  filter(!is.na(nidx)) %>%
  group_by(gene, Age2, nidx, maxvaf) %>%
  summarise(C = n()) %>%
  ungroup() %>%
  rename(n = nidx) %>%
  complete(gene, Age2, nesting(n), fill = list(C = 0)) %>%
  #filter(C > 0) %>%
  filter(!is.na(n)) %>%
  filter(n < maxvaf)

prior1 <- prior(normal(5, 2), nlpar = "A", lb = 0.00001) +
  prior(normal(0, 5), nlpar = "B")
nchains <- args$threads

message("Fit model...")
fitgeneage <- brm(bf(C ~ (A / n) * exp(-n / exp(B)),
               A ~ 1 + (Age2|gene),
               B ~ 1 + (Age2|gene),
               nl = TRUE),
            data = mydat,
            prior = prior1,
            family = gaussian,
            control = list(adapt_delta = 0.99),
            chains = nchains,
            cores = nchains,
            iter = args$its)

fitgeneage <- add_criterion(fitgeneage, c("loo", "waic", "R2"))
print(fitgeneage)
print(bayes_R2(fitgeneage))

message("")
message("###########################################################")
message("Fit model pergene and per age for synonymous variants")

mydat <- df %>%
  filter(impact == "Synonymous") %>%
  filter(str_detect(gene, "NOTCH1|TP53")) %>%
  #group_by(gene, Age2) %>%
  #mutate(nmuts = n(),          maxvaf = quantile(sumvaf, args$quantile)) %>%
  #ungroup() %>%
  #filter(nmuts > 9) %>%
  mutate(nidx = midcut(sumvaf, args$minvaf, 2, args$binsize)) %>%
  filter(!is.na(nidx)) %>%
  group_by(gene, Age2, nidx, maxvaf) %>%
  summarise(C = n()) %>%
  ungroup() %>%
  rename(n = nidx) %>%
  complete(gene, Age2, nesting(n), fill = list(C = 0)) %>%
  #filter(C > 0) %>%
  filter(!is.na(n)) %>%
  filter(n < maxvaf)

prior1 <- prior(normal(5, 2), nlpar = "A", lb = 0.00001) +
  prior(normal(0, 5), nlpar = "B")
nchains <- args$threads

message("Fit model...")
fitgeneagesyn <- brm(bf(C ~ (A / n) * exp(-n / exp(B)),
               A ~ 1 + (Age2|gene),
               B ~ 1 + (Age2|gene),
               nl = TRUE),
            data = mydat,
            prior = prior1,
            family = gaussian,
            control = list(adapt_delta = 0.99),
            chains = nchains,
            cores = nchains,
            iter = args$its)

fitgeneagesyn <- add_criterion(fitgeneage, c("loo", "waic", "R2"))
print(fitgeneagesyn)
print(bayes_R2(fitgeneagesyn))

message("")
message("###########################################################")
message("Saving file")
out<- list(gene = fitgene, age = fitage,
           geneage = fitgeneage,
           geneagesyn = fitgeneagesyn,
            agesynon = fitagesynon)
saveRDS(out, args$output)
