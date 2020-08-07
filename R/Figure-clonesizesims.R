library(brms)
library(tidyverse)
library(jcolors)
library(argparse)
library(tidybayes)
library(bayesplot)
library(modelr)
library(cowplot)

message("Read in data")
sims <- read_csv(snakemake@input$simulationdata, guess_max = 10^5) %>%
  mutate(condition = factor(group_indices(., gene))) %>%
  mutate(A = f / snakemake@config$params$rho)

fits <- readRDS(snakemake@input$simulationfits)
fit <- fits$fit
options(scipen = 999)

message("Transform data")
Nt0 <- function(rlam = 0.5, delta = 0.0, t = 10.0){
  Nt <- 1 + rlam * t
  return(Nt)
}

Nt <- function(rlam = 0.5, delta = 0.001, t = 10.0){
  top <- (1 + delta) * exp(2*rlam*delta*t) - (1 - delta)
  bottom <- 2 * delta
  Nt <- top / bottom
  return(Nt)
}

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
  mutate(fidx = midcut(A, snakemake@config$params$binsize, 1, snakemake@config$paramsbinsize)) %>%
  filter(!is.na(fidx)) %>%
  group_by(condition) %>%
  mutate(maxA = max(A)) %>%
  group_by(fidx, gene, delta, rlam, t, Nsims, mu, N0, condition, maxA) %>%
  summarise(C = n()) %>%
  ungroup() %>%
  rename(n = fidx) %>%
  complete(gene, nesting(n), fill = list(C = 0)) %>%
  fill(delta, rlam, t, Nsims, mu,N0, condition, .direction = "down") %>%
  mutate(C = C) %>%
  mutate(A = mu * N0 * Nsims * (snakemake@config$paramsbinsize),
         B = ifelse(delta == 0,
                    Nt0(rlam = rlam, delta = delta, t = t),
                    Nt(rlam = rlam, delta = delta, t = t)),
         B = B / snakemake@config$paramsrho,
         logB = log(B)) %>%
  #mutate(n = n) %>%
  mutate(Ctheory = (A / n) * (1 / (1 + delta)) * exp(-n / B)) %>%
  group_by(gene) %>%
  mutate(P = C / sum(C)) %>%
  mutate(Ptheory = (1 / n) * (1 / log(B)) * exp(-n / B)) %>%
  mutate(Ptheory = Ptheory / sum(Ptheory)) %>%
  ungroup() %>%
  filter(n < maxA)

params <- distinct(mydat, gene, condition, A, logB, B, t, delta, rlam, mu, Nsims, N0) %>%
  mutate(condition = gene) %>%
  mutate(yax = paste0(t, ", ", delta)) %>%
  mutate(yaxmu = paste0(round(mu, 6)))
Nsims <- params$Nsims[1]

message("Summarize parameter fits")

gA <- fit %>%
  spread_draws(r_gene__A[condition,], b_A_Intercept, sd_gene__A_Intercept) %>%
  mutate(condmean = r_gene__A + b_A_Intercept) %>%
  left_join(params) %>%
  mutate(condmean = condmean / (Nsims * (snakemake@config$paramsbinsize))) %>%
  mutate(yax = paste0(mu)) %>%
  ggplot(aes(y = yaxmu)) +
  scale_color_brewer() +
  stat_pointintervalh(aes(x = condmean), alpha = 0.7,
                      .width = c(.66, .95), position = position_nudge(y = 0.0)) +
  # data
  geom_point(aes(x = A / (Nsims * (snakemake@config$paramsbinsize))), data = params, col = "firebrick", fill = "white", shape = 21, size = 2) +
  xlab(~n[0]~mu/r~lambda~rho) +
  coord_flip() +
  theme_cowplot() +
  ylab(expression("Input "~mu)) +
  theme_cowplot() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

gB <- fit %>%
  spread_draws(r_gene__B[condition,], b_B_Intercept, sd_gene__B_Intercept) %>%
  mutate(condmean = r_gene__B + b_B_Intercept) %>%
  mutate(prediction = rnorm(n(), condmean, sd_gene__B_Intercept)) %>%
  left_join(params) %>%
  mutate(condmean = exp(condmean)) %>%
  mutate(yax = paste0(t, ", ", delta)) %>%
  ggplot(aes(y = yax)) +
  stat_pointintervalh(aes(x = condmean), alpha = 0.7,
                      .width = c(.66, .95), position = position_nudge(y = -0.0)) +
  geom_point(aes(x = B), data = params, col = "firebrick", fill = "white", shape = 21, size = 2) +
  theme_cowplot() +
  ggtitle("") +
  theme(legend.position = "none") +
  scale_x_log10() +
  xlab(expression("N(t)/"~rho)) +
  coord_flip() +
  ylab(expression("time (Years), "~Delta)) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

message("Plot fits")

gfits <- mydat %>%
    filter(delta == snakemake@config$params$delta) %>%
    data_grid(n = unique(mydat$n), gene) %>%
    add_predicted_draws(fit, n = 1000) %>%
    left_join(., mydat) %>%
    filter(n < maxA) %>%
    ggplot(aes(x = n, y = C)) +
    stat_lineribbon(aes(y = .prediction), .width = c(.95, .8, .5), color = "#08519C") +
    geom_point(data = mydat %>% filter(delta == 0.1)) +
    scale_fill_brewer() +
    facet_wrap(~ t, ncol = 7, scales = "free") +
    theme_cowplot() + scale_y_log10() +
    xlab("Area") +
    ylab("Counts") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))

message("Save plot fits")
g <- plot_grid(gA, gB, gfits, ncol = 1, labels = c("a", "b", "c"), align = T)
save_plot(snakemake@output$suppfigure[1], g, base_height = 12, base_width = 10)

val <- params %>%
  summarise(x = mean(mu) * mean(Nsims) * mean(N0)) %>%
  pull(x)
val <- val * snakemake@config$params$binsize

infval <- fit %>%
  spread_draws(b_A_Intercept) %>%
  median_qi()

message(paste0("Population intercept coefficient: ", round(val, 3)))
message(paste0("Inferred population intercept coefficient: ", round(infval$b_A_Intercept, 3),
               " (", round(infval$.lower, 3) , ", ", round(infval$.upper, 3), ")"))


sims %>%
  filter(delta == 0.0) %>%
  group_by(t) %>%
  summarise(x= n() / first(Nsims)) %>%
  ggplot(aes(x = t, y = x)) + geom_point() +
  scale_y_log10()#+geom_smooth(method = "lm")

f2 <- sims %>%
  filter(delta == 0.0) %>%
  group_by(t, mu) %>%
  summarise(x= n() / first(Nsims))
summary(lm(x~log(t), f2))


message("Hitchikers")

dfparams <- data.frame(tstr = unique(fits$hitchikers$data$condition), Nsims = 10^4, mu = 0.01 / 3, N0 = 10^3,
                       rho = snakemake@config$params$rho, rlam = 0.5) %>%
  separate(tstr, c(NA, "t"), "_") %>%
  mutate(t = as.numeric(t), tstr = paste0("t_", t)) %>%
  mutate(A = N0 * mu * Nsims * (snakemake@config$params$binsize),
         B = (1 + rlam * t)/5000) %>%
  mutate(yax = t)

gA <- fits$hitchikers %>%
  spread_draws(r_condition__A[condition, ], b_A_Intercept) %>%
  mutate(condmean = r_condition__A + b_A_Intercept) %>%
  #left_join(params) %>%
  mutate(condmean = condmean / (Nsims * (snakemake@config$params$binsize))) %>%
  mutate(yax = as.numeric(str_split(condition, "_")[[1]][2])) %>%
  ggplot(aes(y = yax)) +
  scale_color_brewer() +
  stat_pointintervalh(aes(x = condmean), alpha = 0.7,
                      .width = c(.66, .95), position = position_nudge(y = 0.0)) +
  # data
  geom_point(aes(x = A), data = dfparams, col = "firebrick", fill = "white", shape = 21, size = 2) +
  xlab(~n[0]~mu/~rho) +
  coord_flip() +
  theme_cowplot() +
  ylab(expression("Input "~mu)) +
  theme_cowplot() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

gB <- fits$hitchikers %>%
  spread_draws(r_condition__B[condition, ], b_B_Intercept) %>%
  mutate(condmean = r_condition__B + b_B_Intercept) %>%
  #left_join(params) %>%
  mutate(condmean = exp(condmean)) %>%
  mutate(yax = as.numeric(str_split(condition, "_")[[1]][2])) %>%
  ggplot(aes(y = yax)) +
  stat_pointintervalh(aes(x = condmean), alpha = 0.7,
                      .width = c(.66, .95), position = position_nudge(y = -0.0)) +
  geom_point(aes(x = B), data = dfparams, col = "firebrick", fill = "white", shape = 21, size = 2) +
  theme_cowplot() +
  ggtitle("") +
  theme(legend.position = "none") +
  #scale_x_log10() +
  xlab(expression("N(t)/"~rho)) +
  coord_flip() +
  ylab(expression("time (Years)")) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

mydathitch <- fits$hitchikers$data %>%
  group_by(condition) %>%
  mutate(maxA = max(n))

gfits <- fits$hitchikers$data %>%
  #filter(delta == 0.1) %>%
  data_grid(n = unique(fits$hitchikers$data$n), condition) %>%
  add_predicted_draws(fits$hitchikers, n = 1000) %>%
  left_join(., mydathitch) %>%
  filter(n < maxA) %>%
  ggplot(aes(x = n, y = C)) +
  stat_lineribbon(aes(y = .prediction), .width = c(.95, .8, .5), color = "#08519C") +
  geom_point(data = mydathitch) +
  scale_fill_brewer() +
  facet_wrap(~ condition, ncol = 3, scales = "free") +
  theme_cowplot() + scale_y_log10() +
  xlab("Area") +
  ylab("Counts") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

save_plot(snakemake@output$suppfigure[2], plot_grid(plot_grid(gA, gB, ncol = 1, labels = c("a", "b")),
                                        gfits, ncol = 2, labels = c("", "c")),
          base_height = 7, base_width = 20)
