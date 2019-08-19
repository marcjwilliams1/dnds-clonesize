#load packages, you will need to install these if you don't have them already
library(cowplot)
library(tidyverse)
theme_set(theme_cowplot())
library(jcolors)
library(forcats)
library(ggforce)
library(Hmisc)

library(argparse)

parser <- ArgumentParser(description = "Generate Final Figures")
parser$add_argument('--figure', type='character',
                    help="Outpute figure files")
parser$add_argument('--suppfigures', type='character',
                    help="Outpute figure files", nargs = "+")
parser$add_argument('--baseline', type='character',
                    help="baseline dn/ds values")
parser$add_argument('--dndsclonality', type='character',
                    help="Mutation VAFs with clonality assignment")
parser$add_argument('--dndsclonality_percancertype', type='character',
                    help="dn/ds as a function of clonality per cancertype")
parser$add_argument('--tcgadata', type='character',
                    help="TCGA data file")
args <- parser$parse_args()

message("Generating Figure 4...")
message("\t Reading in data...")

dfsubc <- read_csv(args$tcgadata, col_types = cols()) %>%
    mutate(clonality = case_when(
                    MCN < 0.5 ~ "Subclonal",
                    MCN > 0.5 & MCN < 1.5 ~ "Clonal",
                    MCN > 1.5 ~ "Amplified"))
baseline <- read_csv(args$baseline, col_types = cols())
dfdnds <- read_csv(args$dndsclonality, col_types = cols())
dfdnds.cancertype <- read_csv(args$dndsclonality_percancertype, col_types = cols())

print(head(dfdnds))
print(head(baseline))

message("Generating figure")

getdndsplot <- function(dfdnds, baseline, mutname, muttype, plottitle = "",
                        ylabel = "dN/dS", ylims = c(0.8, 2.0), pointsize = 2, linesize = 1){

    dftemp <- dfdnds %>%
        left_join(., baseline, by = "name") %>%
        mutate(mle = mle - (dnds_bl - 1),
              cilow = cilow - (dnds_bl - 1),
              cihigh = cihigh - (dnds_bl - 1)) %>%
        filter(name == mutname, mutationtype ==  muttype)
    out <- dfdnds %>%
    filter(name == mutname, mutationtype ==  muttype) %>%
    ggplot(aes(x = clonality, y = mle, ymin = cilow, ymax = cihigh, col = clonality)) +
    geom_point(size = pointsize) +
    geom_linerange(size = linesize) +
    geom_hline(yintercept = 1.0, lty = 2) +
    background_grid(major = "xy", minor = "none") +
    ylab(ylabel) + xlab("") + ggtitle(plottitle) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    ylim(ylims)
    return(out)
}

gMCN <- dfsubc %>%
    na.omit() %>%
    mutate(clonality = factor(clonality,
                              levels = c("Subclonal", "Clonal", "Amplified"))) %>%
    ggplot(aes(x = MCN, fill = clonality)) +
    geom_histogram(bins = 100, position = "identity") +
    geom_vline(xintercept = 1.0, lty = 2, col = "firebrick", alpha = 0.5) +
    xlim(c(0, 3.0)) +
    theme(legend.position = c(0.5, 0.8), legend.title = element_blank()) +
    xlab("Mutation copy number") + ylab("Counts") +
    scale_fill_jcolors(palette = "pal9") #+scale_color_jcolors(palette = "pal9")

dfdnds <- dfdnds %>%
    mutate(clonality = factor(clonality,
                              levels = c("Subclonal", "Clonal", "Amplified")))

gdndsdriver <- getdndsplot(dfdnds, baseline, "wmis", "drivers", "Driver mutations", "dN/dS - Missense \n (Drivers)", ylims = c(0.5, 2.5)) +
                scale_color_jcolors(palette = "pal9") + theme(legend.position = "none") + ggtitle("")
gdndsdrivernon <- getdndsplot(dfdnds, baseline, "wnon", "drivers", "Driver mutations", "dN/dS - Nonsense \n (Drivers)", ylims = c(0.5, 6.0)) +
                scale_color_jcolors(palette = "pal9") + theme(legend.position = "none") + ggtitle("")

gclonality <- cowplot::plot_grid(gMCN, gdndsdriver, gdndsdrivernon, ncol = 3, labels = c("a", "b", "c"))

g1 <- dfdnds.cancertype %>%
    mutate(clonality = factor(clonality,
                                  levels = c("Subclonal", "Clonal", "Amplified"))) %>%
    filter(!cancertype %in% c("THCA", "PRAD", "KIRC")) %>%
    filter(mutationtype == "drivers", name == "wmis") %>%
    ggplot(aes(x = cancertype, y = mle, ymin = cilow, ymax = cihigh, col = clonality)) +
    geom_point(position = position_dodge(0.4)) +
    geom_linerange(position = position_dodge(0.4)) +
    scale_color_jcolors(palette = "pal9") +
    #scale_y_log10() +
    xlab("") +
    ylab("dN/dS - Missense") +
    geom_hline(yintercept = 1.0, lty = 2) +
    theme(axis.text.x = element_text(angle = 30, hjust = 1)) +
    theme(legend.position = "none") +
    theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())

g2 <- dfdnds.cancertype %>%
    mutate(clonality = factor(clonality,
                                  levels = c("Subclonal", "Clonal", "Amplified"))) %>%
    filter(!cancertype %in% c("THCA", "PRAD", "KIRC")) %>%
    filter(mutationtype == "drivers", name == "wnon") %>%
    ggplot(aes(x = cancertype, y = mle, ymin = cilow, ymax = cihigh, col = clonality)) +
    geom_point(position = position_dodge(0.4)) +
    geom_linerange(position = position_dodge(0.4)) +
    scale_color_jcolors(palette = "pal9") +
    #scale_y_log10() +
    xlab("") +
    ylab("dN/dS - Nonsense") +
    geom_hline(yintercept = 1.0, lty = 2) +
    theme(axis.text.x = element_text(angle = 30, hjust = 1)) +
    theme(legend.position = c(0.86, 0.8), legend.title = element_blank())

gpercancer <- plot_grid(g1, g2, ncol = 1, rel_heights = c(1.0, 1.4), align = "v", labels = c("d", "e"))

gclonality <- plot_grid(gclonality, gpercancer, ncol = 1)

save_plot(args$figure, gclonality, base_height = 10, base_width = 10)
gclonality

message("Generate supp figure...")


#only include cancer types with at least 100 samples
cancertypes <- dfsubc %>%
    distinct(sampleid, cancertype) %>%
    group_by(cancertype) %>%
    summarise(n = n()) %>%
    filter(n >= 100) %>%
    filter(!cancertype %in% c("THCA", "PRAD")) %>%
    pull(cancertype)

plotlist <- list()
i <- 1

for (ct in cancertypes){

    dfdndsct <- dfdnds.cancertype %>%
        filter(cancertype == ct) %>%
        mutate(clonality = factor(clonality,
                              levels = c("Subclonal", "Clonal", "Amplified")))

    gvafall <- dfsubc %>%
        filter(cancertype == ct) %>%
        mutate(clonality = factor(clonality,
                                  levels = c("Subclonal", "Clonal", "Amplified"))) %>%
        ggplot(aes(x = MCN, fill = clonality)) +
        geom_histogram(bins = 100, position = "identity") +
        geom_vline(xintercept = 1.0, lty = 2, col = "firebrick", alpha = 0.5) +
        xlim(c(0, 3.0)) +
        theme(legend.position = c(0.5, 0.8), legend.title = element_blank()) +
        xlab("Mutation copy number") + ylab("Counts") +
        scale_fill_jcolors(palette = "pal9") +
        ggtitle(ct)


    lim <- round(max(filter(dfdndsct, name == "wmis", mutationtype == "all")$cihigh)) + 0.5
    gdndsdriver <- getdndsplot(dfdndsct, baseline, "wmis", "all", "Driver mutations", "dN/dS - Missense", ylims = c(0.8, lim)) +
                    scale_color_jcolors(palette = "pal9") + theme(legend.position = "none") + ggtitle("")
    lim <- round(max(filter(dfdndsct, name == "wnon", mutationtype == "all")$cihigh)) + 0.5
    gdndsdrivernon <- getdndsplot(dfdndsct, baseline, "wnon", "all", "Driver mutations", "dN/dS - Nonsense", ylims = c(0.8, lim)) +
                    scale_color_jcolors(palette = "pal9") + theme(legend.position = "none") + ggtitle("")


    gclonality <- plot_grid(gvafall, gdndsdriver, gdndsdrivernon, ncol = 3)

    plotlist[[i]] <- gclonality
    i <- i + 1
}

#save_plot("FinalFigures/FigureS6.pdf", plot_grid(plotlist = plotlist, ncol = 2), base_height = 20,base_width = 20)

plotlist <- list()
i <- 1

for (ct in cancertypes){

    dfdndsct <- dfdnds.cancertype %>%
        filter(cancertype == ct) %>%
        mutate(clonality = factor(clonality,
                              levels = c("Subclonal", "Clonal", "Amplified")))

    gvafall <- dfsubc %>%
        filter(cancertype == ct) %>%
        mutate(clonality = factor(clonality,
                                  levels = c("Subclonal", "Clonal", "Amplified"))) %>%
        ggplot(aes(x = MCN, fill = clonality)) +
        geom_histogram(bins = 100, position = "identity") +
        geom_vline(xintercept = 1.0, lty = 2, col = "firebrick", alpha = 0.5) +
        xlim(c(0, 3.0)) +
        theme(legend.position = c(0.5, 0.8), legend.title = element_blank()) +
        xlab("Mutation copy number") + ylab("Counts") +
        scale_fill_jcolors(palette = "pal9") +
        ggtitle(ct)


    lim <- round(max(filter(dfdndsct, name == "wmis")$cihigh)) + 0.5
    gdndsdriver <- getdndsplot(dfdndsct, baseline, "wmis", "drivers", "Driver mutations", "dN/dS - Missense", ylims = c(0.0, lim)) +
                    scale_color_jcolors(palette = "pal9") + theme(legend.position = "none") + ggtitle("")
    lim <- round(max(filter(dfdndsct, name == "wnon")$cihigh)) + 0.5
    gdndsdrivernon <- getdndsplot(dfdndsct, baseline, "wnon", "drivers", "Driver mutations", "dN/dS - Nonsense", ylims = c(0.0, lim)) +
                    scale_color_jcolors(palette = "pal9") + theme(legend.position = "none") + ggtitle("")


    gclonality <- plot_grid(gvafall, gdndsdriver, gdndsdrivernon, ncol = 3)

    plotlist[[i]] <- gclonality
    i <- i + 1
}

save_plot(args$suppfigure[1], plot_grid(plotlist = plotlist, ncol = 2), base_height = 20, base_width = 20)
