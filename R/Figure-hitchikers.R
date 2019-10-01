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
parser$add_argument('--oesophagusdata_all', type='character',
                    help="Mutations and VAF for oesophagus data including sample info")
parser$add_argument('--oesophagusmetadata', type='character',
                    help=" oesophagus meta data")
args <- parser$parse_args()

# args <- list(simulationdata = "~/Documents/apocrita/BCInew/marc/dnds/dnds-clonesize/results/simulations/clonesize_hitchikers.csv",
#              datamodelfits = "~/Documents/apocrita/BCInew/marc/dnds/dnds-clonesize/results/dataforfigures/data-clonesizefit-models.Rdata",
#              oesophaagusdata = "~/Documents/apocrita/BCInew/marc/dnds/dnds-clonesize/data/oesophagus/esophagus.csv",
#              oesophaagusdata_all = "~/Documents/apocrita/BCInew/marc/dnds/dnds-clonesize/data/oesophagus/aau3879_TableS2.xlsx",
#              oesophaagusmetadata = "~/Documents/apocrita/BCInew/marc/dnds/dnds-clonesize/data/oesophagus/patient_info.xlsx",
#              rho = 5000,
#              datafits = "~/Documents/apocrita/BCInew/marc/dnds/dnds-clonesize/results/dataforfigures/data-clonesizefit.Rdata",
#              binsize = 0.005)

dfsims <- read_csv(args$simulationdata)

dfdata <- read_csv(args$oesophagusdata)
donor <- readxl::read_xlsx(args$oesophagusmetadata, skip = 1) %>%
  dplyr::rename(donor = PD)
dfdata <- left_join(dfdata, donor)


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

t <- 73.5

mydat <- dfsims %>%
  #filter(tend == t) %>%
  mutate(A = f / args$rho) %>%
  mutate(fidx = midcut(A, args$binsize, 1, args$binsize)) %>%
  filter(!is.na(fidx)) %>%
  group_by(muNS, muS, muttype, tend) %>%
  mutate(maxA = max(A)) %>%
  ungroup() %>%
  group_by(fidx, N0, Nsims, muNS, muS, maxA, muttype, tend) %>%
  summarise(C = n()) %>%
  ungroup() %>%
  rename(n = fidx) %>%
  mutate(gene = paste0(muttype, muNS)) %>%
  complete(gene, nesting(n), fill = list(C = 0)) %>%
  fill(Nsims, muNS, muS, N0, tend, .direction = "down") %>%
  filter(n < maxA) %>%
  mutate(C = C / (Nsims * args$binsize))

mydat %>%
  filter(muttype == "Syn", muNS == 0.01) %>%
  ggplot(aes(x = n, y = C)) +
  geom_line() +
  #geom_bar(stat="identity") +
  scale_x_log10() +
  scale_y_log10() +
  facet_grid(~tend, scales = "free_y")

dfsims %>%
  #filter(muttype == "NonSyn") %>%
  filter(muNS == 0.001) %>%
  mutate(f = f / 5000) %>%
  ggplot(aes(x = f)) +
  #geom_bar(stat = "identity") +
  geom_histogram(bins = 100) +
  scale_x_log10() +
  facet_grid(muttype ~ tend, scales = "free_y")


#theory from Watson et al.
hitchikeSFS <- function(N0, tau, t, s, n, muNS, muS){
  var1 <- tau * muNS * muS / s
  var2 <- (exp(s * t) - 1) / s
  var3 <- n
  var4 <- (n + 1/s) ^2
  return(var1 * ((var2 - var3) / var4))
}

nonhitchikeSFS <- function(N0, t, s, n, muS, rlam){
  var1 <- N0 * muS / n
  var2 <- 1 + rlam * t
  return(var1 * exp(- n /var2))
}

N0 <- mydat$N0[1]
rlam <- 0.5
tau <- rlam * t
s <- 0.05

mydat$th <- hitchikeSFS(N0, rlam * mydat$tend, mydat$tend, s, mydat$n * args$rho, 2 * mydat$muNS, 2 * mydat$muS) * args$rho
mydat$th2 <- nonhitchikeSFS(N0, mydat$tend, s, mydat$n * args$rho, 2 * mydat$muS, rlam) * args$rho


mydat %>%
  filter(muttype != "NonSyn", muNS == 0.01, muttype == "SynHitch") %>%
  ggplot(aes(x = n, y = C)) +
  geom_line() +
  geom_line(aes(y = th), col = "red") +
  #geom_bar(stat="identity") +
  scale_x_log10() +
  scale_y_log10() +
  facet_grid(muNS~muttype, scales = "free_y")

mydat %>%
  filter(muttype != "NonSyn", muNS == 0.01) %>%
  ggplot(aes(x = n, y = C)) +
  geom_line() +
  geom_line(aes(y = th), col = "red") +
  geom_line(aes(y = th2), col = "blue") +
  #geom_bar(stat="identity") +
  scale_x_log10() +
  scale_y_log10() +
  facet_grid(muttype~tend, scales = "free_y")

dfsims %>%
  mutate(A = f / args$rho) %>%
  filter(muttype == "Syn", muNS == 0.01) %>%
  ggplot(aes(x = A)) +
  geom_histogram(bins = 100) +
  facet_wrap(~tend) +
  scale_x_log10() +
  scale_y_log10()


# Plot data per patch
message("Read in data from xlsx file")

dfall <- readxl::read_xlsx(args$oesophaagusdata_all, sheet = 1, skip = 16)

dfallcounts <- dfall %>%
  mutate(muttype = case_when(
    str_detect(gene, "NOTCH1") & impact != "Synonymous" ~ "NOTCH1",
    str_detect(gene, "TP53") & impact != "Synonymous" ~ "TP53",
    impact == "Synonymous" ~ "Synonymous"
  )) %>%
  filter(!is.na(muttype)) %>%
  group_by(sampleID, muttype) %>%
  summarise(vaf = mean(vaf), n = n()) %>%
  ungroup() %>%
  complete(sampleID, nesting(muttype), fill = list(vaf = 0.0, n = 0))

nmuts <- dfallcounts %>%
  dplyr::select(-vaf) %>%
  spread(muttype, n)

vafmuts <- dfallcounts %>%
  dplyr::select(-n) %>%
  spread(muttype, vaf)

patch <- full_join(nmuts, vafmuts, by = "sampleID", suffix = c("_n", "_vaf"))

(g1 <- patch %>%
  mutate(NOTCH1_n = ifelse(NOTCH1_n > 9, 10, NOTCH1_n),
         TP53_n = ifelse(TP53_n > 9, 10, TP53_n)) %>%
  ggplot(aes(x = factor(NOTCH1_n), y = Synonymous_n)) +
  geom_boxplot(width = 0.2, alpha = 0.8) +
  geom_violin(alpha = 0.5, fill = "firebrick4") +
  theme_cowplot() +
  xlab("# NOTCH1 Non-synonymous mutations") +
  ylab("# Synonymous mutations"))

(g2 <- patch %>%
    mutate(NOTCH1_n = ifelse(NOTCH1_n > 9, 10, NOTCH1_n),
           TP53_n = ifelse(TP53_n > 9, 10, TP53_n)) %>%
    ggplot(aes(x = factor(TP53_n), y = Synonymous_n)) +
    geom_boxplot(width = 0.2, alpha = 0.8) +
    geom_violin(alpha = 0.5, fill = "firebrick4") +
    theme_cowplot() +
    xlab("# TP53 Non-synonymous mutations") +
    ylab(""))

(g3 <- patch %>%
    mutate(NOTCH1_n = ifelse(NOTCH1_n > 9, 10, NOTCH1_n),
           TP53_n = ifelse(TP53_n > 9, 10, TP53_n)) %>%
    ggplot(aes(x = factor(NOTCH1_n), y = Synonymous_vaf)) +
    geom_boxplot(width = 0.2, alpha = 0.8) +
    geom_violin(alpha = 0.5, fill = "deepskyblue4") +
    theme_cowplot() +
    scale_y_log10() +
    xlab("# NOTCH1 Non-synonymous mutations") +
    ylab("VAF synonymous mutations"))

(g4 <- patch %>%
    mutate(NOTCH1_n = ifelse(NOTCH1_n > 9, 10, NOTCH1_n),
           TP53_n = ifelse(TP53_n > 9, 10, TP53_n)) %>%
    ggplot(aes(x = factor(TP53_n), y = Synonymous_vaf)) +
    geom_boxplot(width = 0.2, alpha = 0.8) +
    geom_violin(alpha = 0.5, fill = "deepskyblue4") +
    theme_cowplot() +
    scale_y_log10() +
    xlab("# TP53 Non-synonymous mutations") +
    ylab(""))

gall <- plot_grid(g1, g2, g3, g4, ncol = 2, align = "hv", labels = c("a", "b", "c", "d"))

message("Save figure")
print(args)
save_plot(args$suppfigures[1], gall, base_height = 8, base_width = 20)

mylm <- lm(Synonymous_n ~ NOTCH1_n, patch)
print(summary(mylm))

mylm <- lm(Synonymous_n ~ TP53_n, patch)
print(summary(mylm))

mylm <- lm(Synonymous_vaf ~ NOTCH1_n, patch)
print(summary(mylm))

mylm <- lm(Synonymous_vaf ~ TP53_n, patch)
print(summary(mylm))



dfallcounts <- dfall %>%
  mutate(muttype = case_when(
    impact != "Synonymous" ~ "NS",
    impact == "Synonymous" ~ "Synonymous"
  )) %>%
  filter(!is.na(muttype)) %>%
  group_by(sampleID, muttype) %>%
  summarise(vaf = mean(vaf), n = n()) %>%
  ungroup() %>%
  complete(sampleID, nesting(muttype), fill = list(vaf = 0.0, n = 0))

nmuts <- dfallcounts %>%
  dplyr::select(-vaf) %>%
  spread(muttype, n)

vafmuts <- dfallcounts %>%
  dplyr::select(-n) %>%
  spread(muttype, vaf)

patch <- full_join(nmuts, vafmuts, by = "sampleID", suffix = c("_n", "_vaf"))


(g1 <- patch %>%
    ggplot(aes(x = NS_n, y = Synonymous_vaf, group = NS_n)) +
    geom_boxplot(width = 0.2, alpha = 0.8) +
    geom_violin(alpha = 0.5, fill = "firebrick4") +
    theme_cowplot() +
    xlab("# NOTCH1 Non-synonymous mutations") +
    ylab("# Synonymous mutations"))

mylm <- lm(Synonymous_n ~ NS_n, patch)
print(summary(mylm))

mylm <- lm(Synonymous_vaf ~ NS_n, patch)
print(summary(mylm))
