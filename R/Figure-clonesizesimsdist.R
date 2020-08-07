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
  mutate(fidx = midcut(A, snakemake@config$params$binsize, 1, snakemake@config$params$binsize)) %>%
  filter(!is.na(fidx)) %>%
  group_by(condition) %>%
  mutate(maxA = max(A)) %>%
  group_by(fidx, gene, delta, deltasample, rlam, t, Nsims, mu, N0, condition, maxA) %>%
  summarise(C = n()) %>%
  ungroup() %>%
  rename(n = fidx) %>%
  complete(gene, nesting(n), fill = list(C = 0)) %>%
  fill(delta, deltasample, rlam, t, Nsims, mu,N0, condition, .direction = "down") %>%
  mutate(C = C) %>%
  mutate(A = mu * N0 * Nsims * (snakemake@config$params$binsize),
         B = ifelse(delta == 0,
                    Nt0(rlam = rlam, delta = deltasample, t = t),
                    Nt(rlam = rlam, delta = deltasample, t = t)),
         B = B / snakemake@config$params$rho,
         logB = log(B)) %>%
  #mutate(n = n) %>%
  mutate(Ctheory = (A / n) * (1 / (1 + deltasample)) * exp(-n / B)) %>%
  group_by(gene) %>%
  mutate(P = C / sum(C)) %>%
  mutate(Ptheory = (1 / n) * (1 / log(B)) * exp(-n / B)) %>%
  mutate(Ptheory = Ptheory / sum(Ptheory)) %>%
  ungroup() %>%
  filter(n < maxA)

params <- distinct(mydat, gene, condition, A, logB, B, t, delta, deltasample, rlam, mu, Nsims, N0) %>%
  mutate(condition = gene) %>%
  mutate(yax = paste0(t, ", ", deltasample)) %>%
  mutate(yaxmu = paste0(round(mu, 5)))
print("Summary of input parameters:")
print(params)
Nsims <- params$Nsims[1]

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


x1 <- rep(0.0, 3333)
x2 <- rexp(3333, 1/0.01)
x3 <- rexp(3333, 1/0.1)
x2 <- rep(0.01, 3333)
x3 <- rep(0.1, 3333)

x <- c(x1, x2, x3)
df <- data.frame(delta = x)
g1 <- df %>%
  ggplot(aes(x = delta)) +
  #geom_density(fill = "firebrick4", col = "white") +
  geom_histogram(bins = 100) +
  theme_cowplot() +
  xlab(~Delta) +
  ylab("Density") +
  xlim(c(-0.02, 0.25)) +
  theme(axis.text.y=element_blank(),
        axis.ticks.y=element_blank())



fitness <- fit %>%
  spread_draws(r_gene__B[condition,], b_B_Intercept, sd_gene__B_Intercept) %>%
  mutate(condmean = r_gene__B + b_B_Intercept) %>%
  mutate(prediction = rnorm(n(), condmean, sd_gene__B_Intercept)) %>%
  left_join(params) %>%
  mutate(condmean = exp(condmean)) %>%
  filter(delta > 0.0) %>%
  mutate(type = "Non-neutral")

neutral <- fit %>%
  spread_draws(r_gene__B[condition,], b_B_Intercept, sd_gene__B_Intercept) %>%
  mutate(condmean = r_gene__B + b_B_Intercept) %>%
  mutate(prediction = rnorm(n(), condmean, sd_gene__B_Intercept)) %>%
  left_join(params) %>%
  mutate(condmean = exp(condmean)) %>%
  filter(delta == 0.0) %>%
  mutate(type = "Neutral")

dfB <- bind_rows(fitness, neutral)

g2 <- dfB %>%
  ggplot(aes(x = t, y = condmean, col = type, group = type)) +
  stat_pointinterval(aes(y = condmean), alpha = 0.7,
                     .width = c(.66, .95), position = position_dodge(width = 0.5)) +
  ylab(expression("N(t)/"~rho)) +
  theme_cowplot() +
  theme_cowplot() +
  scale_y_log10() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  geom_smooth(method = "lm") +
  scale_colour_manual(values = c("brown4", "darkslategray")) +
  theme(legend.position = c(0.1, 0.9),
        legend.title = element_blank()) +
  xlab("Time (years)")


message("Summarize parameter fits")

gA <- fit %>%
  spread_draws(r_gene__A[condition,], b_A_Intercept, sd_gene__A_Intercept) %>%
  mutate(condmean = r_gene__A + b_A_Intercept) %>%
  left_join(params) %>%
  mutate(condmean = condmean / (Nsims * (snakemake@config$params$binsize))) %>%
  mutate(yax = paste0(mu)) %>%
  ggplot(aes(y = yaxmu)) +
  scale_color_brewer() +
  stat_pointintervalh(aes(x = condmean), alpha = 0.7,
                      .width = c(.66, .95), position = position_nudge(y = 0.0)) +
  # data
  geom_point(aes(x = A / (Nsims * (snakemake@config$params$binsize))), data = params, col = "firebrick", fill = "white", shape = 21, size = 2) +
  xlab(~n[0]~mu/~rho) +
  coord_flip() +
  theme_cowplot() +
  scale_y_log10() +
  ylab(expression("Input "~mu)) +
  theme_cowplot() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

gB <- fit %>%
  spread_draws(r_gene__B[condition,], b_B_Intercept, sd_gene__B_Intercept) %>%
  mutate(condmean = r_gene__B + b_B_Intercept) %>%
  mutate(prediction = rnorm(n(), condmean, sd_gene__B_Intercept)) %>%
  left_join(params) %>%
  filter(deltasample > 0.0) %>%
  mutate(condmean = exp(condmean)) %>%
  mutate(yax = paste0(t, ", ", deltasample)) %>%
  ggplot(aes(y = t)) +
  stat_pointintervalh(aes(x = condmean), alpha = 0.7,
                      .width = c(.66, .95), position = position_nudge(y = -0.0)) +
  geom_point(aes(x = B), data = params %>% filter(deltasample > 0.0), col = "firebrick", fill = "white", shape = 21, size = 2) +
  theme_cowplot() +
  ggtitle("") +
  theme(legend.position = "none") +
  xlab(expression("N(t)/"~rho)) +
  scale_x_log10() +
  coord_flip() +
  ylab("Time (Years)") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))



g <- plot_grid(g2, gB, ncol = 2, labels = c("a", "b"))
save_plot(snakemake@output$suppfigures[1], g, base_height = 4, base_width = 9)




average <- sims %>%
  group_by(t, delta) %>%
  summarise(meanclones = mean(f)) %>%
  mutate(Nt = Nt(delta = delta, t = t), theory = (Nt - 1) /log(Nt))
print(average)
