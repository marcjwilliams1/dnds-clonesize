#load packages, you will need to install these if you don't have them already
library(brms)
library(cowplot)
library(tidyverse)
theme_set(theme_cowplot())
library(jcolors)
library(tidybayes)

message("Read in data")
df <- read_csv(snakemake@input$oesophagusdata)
donor <- readxl::read_xlsx(snakemake@input$oesophaguspatientinfo, skip = 1) %>%
  dplyr::rename(donor = PD)
df <- left_join(df, donor)

message("Create data set")

dat <- df %>%
  mutate(area= 2 * sumvaf) %>%
  filter(impact != "Synonymous", aachange != ".", aachange != "-")

message(summary(dat))

message("Define brms paramters")
nchains <- snakemake@threads
its <- 5000
formula <- bf(area ~ Age2 + (1 + Age2|aachange))

genes <- c("NOTCH1", "TP53")
out <- list()
outcoef <- list()

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
                filter(n > 3)
    message(summary(dattemp))
    brms_frechet <- brm(formula,
                    dat = dattemp,
                    family = lognormal(),
                    chains = nchains,
                    cores = nchains,
    		        control = list(adapt_delta = 0.99),
                    iter = its)
    brms_frechet <- add_criterion(brms_frechet,
                                c("loo", "waic", "R2"))
    print(brms_frechet)
    print(bayes_R2(brms_frechet))
    out[[g]] <- brms_frechet

    frechetCoef <- brms_frechet %>%
      spread_draws(r_aachange[gene, var], b_Age2) %>%
      filter(var == "Age2") %>%
      median_qi(coef = r_aachange + b_Age2) %>%
      arrange(desc(coef))
    outcoef[[g]] <- frechetCoef
    message(head(outcoef[[g]]))
}

message("")
message("###########################################################")
message("Saving file")

saveRDS(out, snakemake@output$fits)
saveRDS(outcoef, snakemake@output$coefficients)
