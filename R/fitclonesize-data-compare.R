library(brms)
library(tidyverse)
library(jcolors)
library(argparse)


message("Read in data")
df <- read_csv(snakemake@input$oesophagusdata)
donor <- readxl::read_xlsx(snakemake@input$oesophaguspatientinfo, skip = 1) %>%
  dplyr::rename(donor = PD)
df <- left_join(df, donor)

df <- df %>%
    filter(sumvaf > snakemake@config$paramsdata$minvaf)

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
  filter(str_detect(impact, "Missense|Nonsense")) %>%
  group_by(Age) %>%
  mutate(maxvaf = quantile(sumvaf, snakemake@config$paramsdata$quantile)) %>%
  ungroup() %>%
  mutate(nidx = midcut(sumvaf, snakemake@config$paramsdata$minvaf, 2, snakemake@config$paramsdata$binsize)) %>%
  filter(!is.na(nidx)) %>%
  group_by(Age, nidx, maxvaf) %>%
  summarise(C = n()) %>%
  ungroup() %>%
  rename(n = nidx) %>%
  complete(Age, nesting(n), fill = list(C = 0)) %>%
  #filter(C > 1) %>%
  filter(n < maxvaf)

message("Fit model 1: 1/n with exponential...")
prior1 <- prior(normal(5, 2), nlpar = "A", lb = 0.0001) +
  prior(normal(0, 5), nlpar = "B")
nchains <- snakemake@threads
fitmodel1 <- brm(bf(C ~ (A / n) * exp(-n / exp(B)),
               A ~ 1 + (1|Age),
               B ~ 1 + (1|Age),
               nl = TRUE),
            data = mydat,
            prior = prior1,
            family = gaussian,
            control = list(adapt_delta = 0.9),
            chains = nchains,
            cores = nchains,
            save_all_pars = TRUE,
            iter = snakemake@config$paramsdata$its)

fitmodel1 <- add_criterion(fitmodel1, c("loo", "waic", "R2"))
print(fitmodel1)
print(bayes_R2(fitmodel1))

message("Fit model 2: exponential...")
prior2 <- prior(normal(5, 2), nlpar = "A", lb = 0.0001) +
  prior(normal(0, 5), nlpar = "B")
nchains <- snakemake@threads
fitmodel2 <- brm(bf(C ~ (A) * exp(-n / exp(B)),
               A ~ 1 + (1|Age),
               B ~ 1 + (1|Age),
               nl = TRUE),
            data = mydat,
            prior = prior2,
            family = gaussian,
            control = list(adapt_delta = 0.9),
            chains = nchains,
            cores = nchains,
            save_all_pars = TRUE,
            iter = snakemake@config$paramsdata$its)

fitmodel2 <- add_criterion(fitmodel2, c("loo", "waic", "R2"))
print(fitmodel2)
print(bayes_R2(fitmodel2))

message("Fit model 3: 1/n power law...")
prior3 <- prior(normal(5, 2), nlpar = "A", lb = 0.0001)
nchains <- snakemake@threads
fitmodel3 <- brm(bf(C ~ (A/n),
               A ~ 1 + (1|Age),
               nl = TRUE),
               data = mydat,
               prior = prior3,
               family = gaussian,
               control = list(adapt_delta = 0.9),
               chains = nchains,
               cores = nchains,
	       save_all_pars = TRUE,
               iter = snakemake@config$paramsdata$its)

fitmodel3 <- add_criterion(fitmodel3, c("loo", "waic", "R2"))
print(fitmodel3)
print(bayes_R2(fitmodel3))

message("Fit model 4: exponential + power law exponent")
prior4 <- prior(normal(5, 2), nlpar = "A", lb = 0.0001) +
  prior(normal(0, 5), nlpar = "B") +
  prior(normal(1, 1), nlpar = "C", lb = 0.1, ub = 3)
nchains <- snakemake@threads
fitmodel4 <- brm(bf(C ~ (A / (n ^ C)) * exp(-n / exp(B)),
               A ~ 1 + (1|Age),
               B ~ 1 + (1|Age),
               C ~ 1,
               nl = TRUE),
            data = mydat,
            prior = prior4,
            family = gaussian,
            control = list(adapt_delta = 0.99),
            chains = nchains,
            cores = nchains,
	    save_all_pars = TRUE,
            iter = snakemake@config$paramsdata$its)

fitmodel4 <- add_criterion(fitmodel4, c("loo", "waic", "R2"))
print(fitmodel4)
print(bayes_R2(fitmodel4))

message("")
message("###########################################################")
message("Compare models")

modellist <- list(fullmodel = fitmodel1,
            exponential = fitmodel2,
            powerlaw = fitmodel3,
            fullmodelplus = fitmodel4)

loo_compare_loo <- loo_compare(modellist$powerlaw,
                               modellist$fullmodel,
                               modellist$exponential,
                               modellist$fullmodelplus,
                               criterion = "loo")
loo_compare_waic <- loo_compare(modellist$powerlaw,
                               modellist$fullmodel,
                               modellist$exponential,
                               modellist$fullmodelplus,
                               criterion = "waic")
print(loo_compare_loo)
print(loo_compare_waic)

mweights <- model_weights(modellist$powerlaw, modellist$fullmodel, modellist$exponential, modellist$fullmodelplus, weights = "waic")
print(mweights)

mweights <- model_weights(modellist$powerlaw, modellist$fullmodel, modellist$exponential, weights = "waic")
print(mweights)

library(bayestestR)
comparison <- bayesfactor_models(modellist$powerlaw, modellist$fullmodel, modellist$exponential)
print(comparison)

message("")
message("###########################################################")
message("Saving file")
saveRDS(modellist, snakemake@output[[1]])
