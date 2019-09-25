library(brms)
library(tidyverse)
library(jcolors)
library(argparse)

parser <- ArgumentParser(description = "Fit brms models")
parser$add_argument('--simulationdata', type='character',
                    help="Simulation data for clone sizes")
parser$add_argument('--simulationdatahitchike', type='character',
                    help="Simulation data for clone sizes")
parser$add_argument('--output', type='character',
                    help="Output file")
parser$add_argument('--threads', type='integer',
                    help="Number of threads", default = 1)
parser$add_argument('--rho', type='double',
                    help="Progenitor density", default = 5000.0)
parser$add_argument('--binsize', type='double',
                    help="Binsize for fitting", default = 0.002)
args <- parser$parse_args()


message("Read in data")
sims <- read_csv(args$simulationdata, guess_max = 10^5) %>%
  mutate(condition = factor(group_indices(., gene))) %>%
  mutate(A = f / args$rho)

simshitchike <- read_csv(args$simulationdata, guess_max = 10^5) %>%
  filter(muttype == "Syn", muNS == 0.01) %>%
  mutate(condition = paste0("t_", tend)) %>%
  mutate(A = f / args$rho)

# Functions to bin data
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



mydat <- sims %>%
  mutate(fidx = midcut(A, args$binsize, 1, args$binsize)) %>%
  filter(!is.na(fidx)) %>%
  group_by(condition) %>%
  mutate(maxA = max(A)) %>%
  group_by(fidx, gene, delta, rlam, t, Nsims, mu, N0, condition, maxA) %>%
  summarise(C = n()) %>%
  ungroup() %>%
  rename(n = fidx) %>%
  complete(gene, nesting(n), fill = list(C = 0)) %>%
  fill(delta, rlam, t, Nsims, mu,N0, condition, .direction = "down")
  filter(n < maxA)

mydathitchike <- simshitchike %>%
  mutate(fidx = midcut(A, args$binsize, 1, args$binsize)) %>%
  filter(!is.na(fidx)) %>%
  group_by(condition) %>%
  mutate(maxA = max(A)) %>%
  group_by(fidx, N0, Nsims, muNS, muS, maxA, muttype, tend) %>%
  summarise(C = n()) %>%
  ungroup() %>%
  rename(n = fidx) %>%
  complete(condition, nesting(n), fill = list(C = 0)) %>%
  fill(Nsims, muNS, muS, N0, tend, .direction = "down") %>%
  filter(n < maxA)

params <- distinct(mydat, gene, condition, A, logB, B) %>%
  mutate(condition = paste0("condition", condition))

prior1 <- prior(normal(5, 2), nlpar = "A", lb = 0.0001) +
  prior(normal(0, 5), nlpar = "B")
nchains <- args$threads
its <- 5000

message("Fit model...")
fit1 <- brm(bf(C ~ (A / n) * exp(-n / exp(B)),
               A ~ 1 + (1|gene),
               B ~ 1 + (1|gene),
               nl = TRUE),
            data = mydat,
            prior = prior1,
            family = gaussian,
            control = list(adapt_delta = 0.9),
            chains = nchains,
            cores = nchains,
            iter = its)

fit1 <- add_criterion(fit1, c("loo", "waic", "R2"))
print(fit1)
print(bayes_R2(fit1))

message("Fit model...")
fit2 <- brm(bf(C ~ (A / n) * exp(-n / exp(B)),
               A ~ 1 + (1|condition),
               B ~ 1 + (1|condition),
               nl = TRUE),
            data = mydathitchike,
            prior = prior1,
            family = gaussian,
            control = list(adapt_delta = 0.9),
            chains = nchains,
            cores = nchains,
            iter = its)

fit2 <- add_criterion(fit2, c("loo", "waic", "R2"))
print(fit2)
print(bayes_R2(fit2))

message("")
message("###########################################################")
message("Saving file")
saveRDS(list(fit = fit1, hitchikers = fit2), args$output)
