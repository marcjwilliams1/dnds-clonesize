library(brms)
library(tidyverse)
library(jcolors)
library(argparse)
library(tidybayes)
library(bayesplot)
library(modelr)
library(cowplot)

parser <- ArgumentParser(description = "Plot simulation fits output")
parser$add_argument('--datafits', type='character',
                    help="Data fits")
parser$add_argument('--datamodelfits', type='character',
                    help="Data model fits")
parser$add_argument('--figure', type='character',
                    help="Output figure files")
parser$add_argument('--oesophagusdata', type='character',
                    help="Mutations and VAF for oesophagus data")
parser$add_argument('--oesophagusmetadata', type='character',
                    help=" oesophagus meta data")
parser$add_argument('--suppfigures', type='character',
                    help="Output figure files", nargs = "+")
parser$add_argument('--rho', type='double',
                    help="Progenitor density", default = 5000.0)
parser$add_argument('--binsize', type='double',
                    help="Binsize for fitting", default = 0.002)
args <- parser$parse_args()

#args <- list(datafits = "~/Documents/apocrita/BCInew/marc/dnds/dnds-clonesize/results/dataforfigures/data-clonesizefit.Rdata",
#             datamodelfits = "~/Documents/apocrita/BCInew/marc/dnds/dnds-clonesize/results/dataforfigures/data-clonesizefit-models.Rdata",
#             oesophaagusmetadata = "~/Documents/apocrita/BCInew/marc/dnds/dnds-clonesize/data/oesophagus/patient_info.xlsx",
#             rho = 5000,
#             binsize = 0.002)

message("Read in data")

df <- read_csv(args$oesophagusdata)
donor <- readxl::read_xlsx(args$oesophagusmetadata, skip = 1) %>%
  dplyr::rename(donor = PD)
df <- left_join(df, donor)

fits <- readRDS(args$datafits)
modelfits <- readRDS(args$datamodelfits)

binsize <- fits$age$data$n[2] - fits$age$data$n[1]
totalarea <- mean(unique(df$Nsamples)) * 2

message("Summary of data")
dftemp <- df %>%
  #filter(gene == "NOTCH1") %>%
  mutate(sumvaf = 2*sumvaf) %>%
  mutate(colgene = gene, 
         colgene = ifelse(str_detect(colgene, "NOTCH[0-1+]|TP53") & impact != "Synonymous", colgene, "Other genes"),
         colgene = case_when(
           str_detect(gene, "NOTCH[0-1+]") & impact != "Synonymous" ~ "NOTCH1",
           str_detect(gene, "TP53") & impact != "Synonymous" ~ "TP53",
           !str_detect(gene, "NOTCH[0-1+]|TP53") & impact == "Synonymous" ~ "Synonymous",
           !str_detect(gene, "NOTCH[0-1+]|TP53") & impact != "Synonymous" ~ "Other genes"
         )) %>%
  mutate(colgene = factor(colgene, levels = c("Synonymous", "Other genes", "NOTCH1", "TP53")))

(gdotplot <- dftemp %>%
    filter(!is.na(colgene)) %>%
    mutate(colgene = factor(colgene, levels = c("TP53", "NOTCH1", "Other genes", "Synonymous"))) %>%
    ggplot(aes(x = Age2, y = sumvaf, group = Age2)) +
    #geom_point(alpha = 0.5) +
    geom_jitter(alpha = 0.5, size = 0.4, aes(col = colgene)) +
    scale_y_log10() +
    stat_summary(fill = "white", col = "black", alpha = 0.65,
                 fun.y = mean, geom = "point", size = 3, shape = 22) +
    scale_color_jcolors(pal = "pal7",  name = "") +facet_wrap(~colgene, ncol = 4) +
    background_grid(major = c("xy")) +
    theme_cowplot() +
    theme(legend.position = "none") +
    xlab("Age") +
    ylab("Area"))

message("Plot model selection results")

waicdf <- data.frame(waic = c(modelfits$fullmodel$waic$estimates[3],
                              modelfits$exponential$waic$estimates[3],
                              modelfits$powerlaw$waic$estimates[3]),
                     waicse = c(modelfits$fullmodel$waic$estimates[3,2],
                              modelfits$exponential$waic$estimates[3,2],
                              modelfits$powerlaw$waic$estimates[3,2]),
                     R2 = c(mean(modelfits$fullmodel$R2),
                            mean(modelfits$exponential$R2),
                            mean(modelfits$powerlaw$R2)),
                     R2min = c(quantile(modelfits$fullmodel$R2, 0.025),
                               quantile(modelfits$exponential$R2, 0.025),
                               quantile(modelfits$powerlaw$R2, 0.025)),
                     R2max = c(quantile(modelfits$fullmodel$R2, 0.975),
                               quantile(modelfits$exponential$R2, 0.975),
                               quantile(modelfits$powerlaw$R2, 0.975)),
                     model = c("Exponential \n + Power law", "Exponential", "Power law"),
                     ord = c(1,2,3))

gwaic <- waicdf %>%
  mutate(ymin = waic - 2*waicse, ymax = waic + 2*waicse) %>%
  ggplot(aes(x = fct_reorder(model, ord), y = waic, ymin = ymin, ymax = ymax)) +
  geom_point() +
  geom_linerange() +
  xlab("") +
  ylab("WAIC") +
  theme_cowplot() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


gR2 <- waicdf %>%
  ggplot(aes(x = fct_reorder(model, ord), y = R2, ymin = R2min, ymax = R2max)) +
  geom_point() +
  geom_linerange() +
  xlab("") +
  ylab(expression("Bayesian"~R^2)) +
  theme_cowplot() +
  coord_flip() +
  theme(axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank())

#########################################
# Patient fits
#########################################

nonsynA <- fits$age %>%
  spread_draws(r_Age__A[Age,], b_A_Intercept) %>%
  mutate(condmean = r_Age__A + b_A_Intercept) %>%
  left_join(., donor) %>%
  mutate(type = "Non-Synonymous")

synonA <- fits$agesynon %>%
  spread_draws(r_Age__A[Age,], b_A_Intercept) %>%
  mutate(condmean = r_Age__A + b_A_Intercept) %>%
  left_join(., donor) %>%
  mutate(type = "Synonymous")

dfA <- bind_rows(nonsynA, synonA) %>%
  mutate(condmean2 = condmean) %>%
  mutate(condmean = condmean / (binsize * totalarea))

message("Regression on summaries for A parameters")
testdf <- nonsynA %>%
  median_qi(x = condmean) %>%
  mutate(type = "Non-synonymous") %>%
  bind_rows(synonA %>%
              median_qi(x = condmean) %>%
              mutate(type = "Synonymous")) %>%
  left_join(., donor) %>%
  mutate(smoking = ifelse(Smokinghistory == "Never smoker", "No", "Yes"))
message("Synonymous and smoking")
print(wilcox.test(x ~ smoking, testdf %>% filter(type == "Synonymous")))
message("Non-Synonymous and smoking")
print(wilcox.test(x ~ smoking, testdf %>% filter(type == "Non-synonymous")))
message("Age vs parameter")
print(summary(lm(x ~ Age2*type, testdf)))

nonsynB <- fits$age %>%
  spread_draws(r_Age__B[Age,], b_B_Intercept) %>%
  mutate(condmean = r_Age__B + b_B_Intercept) %>%
  left_join(., donor) %>%
  mutate(condmean = exp(condmean)) %>%
  mutate(type = "Non-Synonymous")

synonB <- fits$agesynon %>%
  spread_draws(r_Age__B[Age,], b_B_Intercept) %>%
  mutate(condmean = r_Age__B + b_B_Intercept) %>%
  left_join(., donor) %>%
  mutate(condmean = exp(condmean)) %>%
  mutate(type = "Synonymous")

dfB <- bind_rows(nonsynB, synonB)


message("Regression on summaries for B parameters")
testdf <- nonsynB %>%
  median_qi(x = condmean) %>%
  mutate(type = "Non-synonymous") %>%
  bind_rows(synonB %>%
              median_qi(x = condmean) %>%
              mutate(type = "Synonymous")) %>%
  left_join(., donor) %>%
  mutate(smoking = ifelse(Smokinghistory == "Never smoker", "No", "Yes"))
message("Synonymous and smoking")
print(wilcox.test(x ~ smoking, testdf %>% filter(type == "Synonymous")))
message("Non-Synonymous and smoking")
print(wilcox.test(x ~ smoking, testdf %>% filter(type == "Non-synonymous")))
message("Age vs parameter")
print(summary(lm(x ~ Age2*type, testdf)))
print(summary(lm(x ~ Age2:type, testdf)))
message("Age vs parameter synonymous")
print(summary(lm(x ~ Age2, testdf %>% filter(type == "Synonymous"))))
message("Age vs parameter non-synonymous")
print(summary(lm(x ~ Age2, testdf %>% filter(type == "Non-synonymous"))))
print(summary(lm(log(x) ~ log(Age2), testdf %>% filter(type == "Non-synonymous"))))


g1 <- dfA %>% 
  ggplot(aes(x = Age, y = condmean, col = type, group = type)) +
  stat_pointinterval(aes(y = condmean), alpha = 0.7,
                     .width = c(.66, .95), position = position_dodge(width = 0.5)) +
  theme_cowplot() +
  #ylab(expression(rho~mu~" (Cells/mm2 X mutations per division)")) +
  ylab(expression(rho~mu~" ")) +
  theme_cowplot() +
  scale_y_log10() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_colour_manual(values = c("brown4", "darkslategray")) +
  theme(legend.position = "none")

g2 <- dfB %>% 
  ggplot(aes(x = Age, y = condmean, col = type, group = type)) +
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
        legend.title = element_blank())

mydatage <- fits$age$data %>%
  left_join(., donor)

gfitsage <- fits$age$data %>%
    data_grid(n = unique(fits$gene$data$n), Age) %>%
    add_predicted_draws(fits$age, n = 1000) %>%
    left_join(., mydatage) %>%
    filter(C > 0) %>%
    ggplot(aes(x = n, y = C)) +
    stat_lineribbon(aes(y = .prediction), .width = c(.95, .8, .5), color = "#08519C") +
    geom_point(data = mydatage) +
    scale_fill_brewer() +
    facet_wrap(~ Age, ncol = 3) +#, scales = "free_x") +
    theme_cowplot() + scale_y_log10() + scale_x_log10() +
    xlab("Area") +
    ylab("Counts") + 
    theme(axis.text.x = element_text(angle = 45, hjust = 1))


gleft <-plot_grid(plot_grid(gdotplot, labels = c("a")), 
                  plot_grid(gwaic, g1 + xlab(""), g2 + xlab("") + scale_y_log10(),
                            labels = c("b", "c", "d"), ncol = 3, rel_widths = c(0.7, 1, 1)), ncol = 1, align = "h")
gall <- plot_grid(gleft, gfitsage, labels = c("", "e"), rel_widths = c(1.0, 0.7))

save_plot(filename = args$figure, gall, base_height = 8, base_width = 18)


#########################################
# Per gene per patient fits
#########################################

genecoefs <- fits$gene %>%
  spread_draws(r_gene__B[condition, int], b_B_Intercept, sd_gene__B_Intercept) %>%
  median_qi(x = exp(r_gene__B + b_B_Intercept)) %>%
  arrange(desc(x)) %>%
  separate(condition, c("gene", "age"), "-") %>%
  mutate(age = as.numeric(age)) 

(genes <- genecoefs %>%
  rename(Age2 = age) %>%
  left_join(., donor) %>%
  ggplot(aes(x = Age, y = x, ymin = .lower, ymax = .upper)) +
  geom_point() +
  geom_linerange() +
  facet_wrap(~gene, scales = "free_y") + scale_y_log10() +
  xlab("Age") + 
  theme_cowplot() +
  ylab(expression("N(t)/"~rho)) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)))

plotlist <- list()
i <- 1
for (mygene in c("NOTCH1", "TP53")){
  print(mygene)
  mydat <- fits$gene$data %>%
    separate(gene, c("gene1", "age"), "-", remove = F) %>%
    mutate(age = as.numeric(age)) %>%
    ungroup() %>%
    filter(gene1 == mygene )
  
  gfits <- fits$gene$data %>%
      separate(gene, c("gene1", "age"), "-", remove = F) %>%
      mutate(age = as.numeric(age)) %>%
      ungroup() %>%
      filter(gene1 == mygene ) %>%
      data_grid(n = unique(fits$gene$data$n), gene) %>%
      add_predicted_draws(fits$gene, n = 1000) %>%
      left_join(., mydat) %>%
      filter(C > 0) %>%
      ggplot(aes(x = n, y = C)) +
      stat_lineribbon(aes(y = .prediction), .width = c(.95, .8, .5), color = "#08519C") +
      geom_point(data = mydat) +
      scale_fill_brewer() +
      facet_wrap(~ age, ncol = 9) +
      theme_cowplot() + scale_y_log10() + scale_x_log10() +
      xlab("Area") +
      ylab("Counts") + 
      theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
      ggtitle(mygene)
  plotlist[[i]] <- gfits
  i <- i + 1
}

gleft <- plot_grid(plotlist = plotlist, ncol = 1)
gall <- plot_grid(gleft, genes, ncol = 2, labels = c("a", "b"))

save_plot(args$suppfigure[1], gall, base_height = 10, base_width = 20)



testdf1 <- nonsynB %>%
  group_by(Age) %>%
  summarise(Nt = median(condmean), sd = sd(condmean),
            Ntlog = median(log(condmean)), sdlog = sd(condmean)) %>%
  ungroup() %>%
  left_join(., donor)

testdf2<- synonB %>%
  group_by(Age) %>%
  summarise(Nt = median(condmean), sd = sd(condmean),
            Ntlog = median(log(condmean)), sdlog = sd(condmean)) %>%
  ungroup() %>%
  left_join(., donor)

lm1 <- lm(log(Nt) ~ log(Age2),testdf1)
summary(lm1)
lm2 <- lm(Nt ~ Age2,testdf1)
summary(lm2)

lm1 <- lm(log(Nt) ~ log(Age2),testdf2)
summary(lm1)
lm2 <- lm(Nt ~ Age2,testdf2)
summary(lm2)