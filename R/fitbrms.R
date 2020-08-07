#load packages, you will need to install these if you don't have them already
library(cowplot)
library(tidyverse)
library(readxl)
theme_set(theme_cowplot())
library(jcolors)
library(brms)

message("Read in data")
df <- read_csv(snakemake@input$oesophagusdata)
donor <- readxl::read_xlsx(snakemake@input$oesophaguspatientinfo, skip = 1) %>%
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
its <- 5000
formula <- bf(area ~ Age2 + (1 + Age2|gene1))

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
		        control = list(adapt_delta = 0.9),
                iter = its)
brms_lognormal <- add_criterion(brms_lognormal,
                            c("loo", "waic", "R2"))
print(brms_lognormal)
print(bayes_R2(brms_lognormal))

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
		        control = list(adapt_delta = 0.9),
                iter = its)
brms_normal <- add_criterion(brms_normal,
                            c("loo", "waic", "R2"))
print(brms_normal)
print(bayes_R2(brms_normal))

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
		            control = list(adapt_delta = 0.8),
                iter = its)
brms_frechet <- add_criterion(brms_frechet,
                            c("loo", "waic", "R2"))
print(brms_frechet)
print(bayes_R2(brms_frechet))

message("")
message("###########################################################")
message("Compare models")

loo_compare_loo <- loo_compare(brms_normal, brms_frechet,
                                brms_lognormal,
                            criterion = "loo")
loo_compare_waic <- loo_compare(brms_normal, brms_frechet,
                                brms_lognormal,
                            criterion = "waic")
print(loo_compare_loo)
print(loo_compare_waic)

message("")
message("###########################################################")
message("Saving file")

out <- list(normal = brms_normal,
            frechet = brms_frechet,
            lognormal = brms_lognormal,
            model_comparison_loo = loo_compare_loo,
            model_comparison_waic = loo_compare_waic)

saveRDS(out, snakemake@$output[[1]])


#library(lme4)
#formula <- bf(area ~ Age2 + (1 + Age2|gene1))

#x <- lmer(log(area) ~ Age2 + (1 + Age2|gene1), data = dat %>% filter(impact == "Missense"))
#as.data.frame(coef(x)$gene1) %>%
  #rownames_to_column(., "gene") %>%
  #arrange(desc(Age2))
