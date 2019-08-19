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
                    help="Outpute figure files", nargs = "+", default = NULL)
parser$add_argument('--intervaldNdStcga', type='character',
                    help="Interval dn/ds TCGA values")
parser$add_argument('--intervaldNdSsim', type='character',
                    help="Interval dn/ds simulation values")
parser$add_argument('--intervaldNdSsimmu', type='character',
                    help="Interval dn/ds simulation values")
parser$add_argument('--intervaldNdSsimpower', type='character',
                    help="Interval dn/ds simulation values - power")
parser$add_argument('--nmutations_gene', type='character',
                    help="Number of mutations per gene")
parser$add_argument('--nmutations_gene_percancertype', type='character',
                    help="Number of mutations per gene per cancertype")
parser$add_argument('--drivergenelist', type='character',
                    help="driver gene list")
args <- parser$parse_args()

message("Generating Figure 5...")
message("\t Reading in data...")

dfdndsInt <- read_csv(args$intervaldNdStcga, col_types = cols())
dfdndsIntSim <- read_csv(args$intervaldNdSsim, col_types = cols())
dfdndsIntpowerSim <- read_csv(args$intervaldNdSsimpower, col_types = cols())
mutspergene <- read_csv(args$nmutations_gene, col_types = cols())
mutspergenepercancer <- read_csv(args$nmutations_gene_percancertype, col_types = cols())
drivers <- read_delim(args$drivergenelist, delim = "\t", col_types = cols(), col_names = F)
dfdndsIntmuSim <- read_csv(args$intervaldNdSsimmu, col_types = cols())

message("Make simulation figure...")

dfplot <- dfdndsIntpowerSim %>%
    group_by(simnum, sample) %>%
    dplyr::summarize(sfit = first(sfit),
    trues = first(trues),
    ndrivers = first(ndrivers),
    npassengers = first(npassengers)) %>%
    mutate(driverspersample = ndrivers / 1000, passengerspersample = npassengers/1000) %>%
    ungroup()

glow <- dfdndsIntpowerSim %>%
    filter(CCF < 0.5) %>%
    filter(simnum == 1, npassengers == 20) %>%
    mutate(trues = paste0(trues)) %>%
    ggplot2::ggplot(ggplot2::aes(x = CCF, y = dnds)) +
    ggplot2::geom_point(size = 0.7, alpha = 0.7) +
    ggplot2::geom_line(ggplot2::aes(y = dndsth), col = "plum4", size = 1.2, alpha = 0.7) +
    ggplot2::ylab("Interval \n dN/dS") + ggplot2::xlab(expression(f[max])) +
    geom_hline(yintercept = 1.0, lty = 2) +
    scale_fill_hue(l=40) +
    background_grid(major = "xy", minor = "none") +
    theme(legend.position="none") +
    ggtitle("# mutations = 20")

ghigh <- dfdndsIntpowerSim %>%
    filter(CCF < 0.5) %>%
    filter(simnum == 1, npassengers == 500) %>%
    mutate(trues = paste0(trues)) %>%
    ggplot2::ggplot(ggplot2::aes(x = CCF, y = dnds)) +
    ggplot2::geom_point(size = 0.7, alpha = 0.7) +
    ggplot2::geom_line(ggplot2::aes(y = dndsth), col = "plum4", size = 1.2, alpha = 0.7) +
    ggplot2::ylab("Interval \n dN/dS") + ggplot2::xlab(expression(f[max])) +
    geom_hline(yintercept = 1.0, lty = 2) +
    scale_fill_hue(l=40) +
    background_grid(major = "xy", minor = "none") +
    theme(legend.position="none") +
    ggtitle("# mutations = 500")

gsummary <- dfplot %>%
    mutate(sample = fct_reorder(as.factor(as.character(npassengers)), npassengers)) %>%
    ggplot(aes(y = sfit, x = sample)) +
    geom_boxplot(alpha = 0.5) +
    geom_sina(col = "plum4", alpha = 0.8, size = 0.7) +
    geom_hline(yintercept = 0.25, lty = 2, col = "firebrick") +
    #geom_hline(yintercept = 0.3, lty = 2, alpha = 0.5) +
    #geom_hline(yintercept = 0.2, lty = 2, alpha = 0.5) +
    xlab("Number of synonymous mutations") +
    ylab("Estimated s")

gpower <- plot_grid(glow, ghigh, gsummary, ncol = 3, rel_widths = c(1,1,1.5), labels = c("a", "b", "c"))

message("Make bottom half of figure...")

labels <- dfdndsIntSim %>%
    mutate(trues = paste0(trues)) %>%
    group_by(trues, sfit) %>%
    dplyr::summarize(y = 1.0*mean(dnds) -0.1 ) %>%
    mutate(x = 0.25) %>%
    mutate(label = paste0(round(sfit,3), " (", trues,")"))

gsim <- dfdndsIntSim %>%
    mutate(trues = paste0(trues)) %>%
    ggplot2::ggplot(ggplot2::aes(x = CCF, y = dnds, fill = trues, group = trues)) +
    ggplot2::geom_point(aes(col = trues)) +
    ggplot2::geom_line(ggplot2::aes(y = dndsth, col = trues),  size = 1.0, alpha = 0.8) +
    ggplot2::ylab("Interval \n dN/dS") + ggplot2::xlab(expression(f[max])) +
    geom_hline(yintercept = 1.0, lty = 2) +
    xlab(expression(f[max])) +
    ylab("Interval dN/dS") +
    scale_fill_hue(l=40) +
    background_grid(major = "xy", minor = "none") +
    scale_fill_brewer(palette="Set1") +
    scale_colour_brewer(palette="Set1") +
    geom_text(data = labels, aes(x = x, y = y, label = label, col = trues), size = 5) +
    theme(legend.position="none") +
    ggtitle("Simulated cohort")


tp53 <- mutspergenepercancer %>%
    filter(gene == "TP53") %>%
    group_by(cancertype) %>%
    summarise(n = sum(n)) %>%
    summarise(x = mean(n)) %>%
    pull(x)

VHL <- mutspergenepercancer %>%
    filter(gene == "VHL", cancertype == "KIRC") %>%
    pull(n) %>%
    sum()

nonsense <- mutspergenepercancer %>%
    filter(impact == "Nonsense") %>%
    group_by(cancertype) %>%
    summarise(x = n()) %>%
    summarise(x = mean(x)) %>%
    pull()

missense <- mutspergenepercancer %>%
    filter(impact == "Missense") %>%
    group_by(cancertype) %>%
    summarise(x = n()) %>%
    summarise(x = mean(x)) %>%
    pull(x)


gsummary <- dfplot %>%
    filter(npassengers != 75) %>%
    #mutate(sample = fct_reorder(as.factor(as.character(npassengers)), npassengers)) %>%
    ggplot(aes(y = sfit, x = npassengers)) +
    geom_boxplot(aes(group = npassengers), alpha = 0.5) +
    #geom_sina(col = "plum4", alpha = 0.8, size = 0.7) +
    geom_hline(yintercept = 0.25, lty = 2, col = "dodgerblue4", size = 1.0) +
    annotate("text", x = 900, y = 0.45, label = "True s", col = "dodgerblue4", size = 6) +
    #geom_vline(xintercept = tp53, lty = 2, col = "firebrick4") +
    annotate("text", x = tp53 + 10, y = 1.5, label = "Mean # subclonal TP53 \n mutations per cancertype",
              col = "firebrick4") +
    annotate("segment", x=tp53,xend=tp53,y=1.5,yend=0.25,arrow=arrow(), color="firebrick4") +
    #geom_vline(xintercept = VHL, lty = 2, col = "darkslategray") +
    annotate("segment", x=VHL + 90,xend=VHL,y=1.0,yend=0.25,arrow=arrow(), color="darkslategray") +
    annotate("text", x = VHL + 90, y = 1.2, label = "# subclonal VHL \n mutations in KIRC",
              col = "darkslategray") +
    annotate("segment", x=missense + 90,xend=missense,y= -0.25,yend=0.25,arrow=arrow(), color="dodgerblue4") +
    annotate("text", x = missense + 90, y = -0.35, label = "Average # subclonal \n driver mutations per cancertype",
              col = "dodgerblue4") +
    xlab("Number of mutations") +
    ylab("Estimated s") +
    ylim(c(-0.5, 1.5)) +
    scale_x_log10(breaks = c(10, 20, 30, 50, 100, 200, 300, 500, 750, 1000)) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))


labels <- dfdndsIntmuSim %>%
    filter(trues == 0.5) %>%
    mutate(trues = paste0("s = ", trues), mud = paste0("ratio = ", mud/0.5)) %>%
    dplyr::group_by(mud, trues) %>%
    dplyr::summarize(s = paste0("s = ", round(first(sfit),2)), A = first(A))

ghitchhike <- dfdndsIntmuSim %>%
    filter(trues == 0.5) %>%
    mutate(trues = paste0("s = ", trues), mud = paste0("ratio = ", mud/0.5)) %>%
    ggplot2::ggplot(ggplot2::aes(x = CCF, y = dnds)) +
    ggplot2::geom_point(col = "grey4", alpha = 0.5, size = 0.5) +
    ggplot2::geom_line(ggplot2::aes(y = dndsth), col = "firebrick",  size = 0.5, alpha = 0.8) +
    ggplot2::ylab("Interval \n dN/dS") + ggplot2::xlab(expression(f[max])) +
    geom_hline(yintercept = 1.0, lty = 2) +
    xlab(expression(f[max])) +
    ylab("Interval dN/dS") +
    background_grid(major = "xy", minor = "none") +
    #scale_color_viridis(option="plasma") + scale_fill_viridis(option="plasma") +
    geom_text(data = labels, aes(x = 0.5, y = Inf, label = s), vjust  = 5) +
    facet_grid(~mud) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))

message("Generate final figure")
gtemp <- cowplot::plot_grid(gsim, gsummary, ncol = 2, align = T, labels = c("a", "b"), rel_widths = c(1, 1.3))
gInt <- plot_grid(gtemp, ghitchhike, ncol = 1, labels = c("", "c"))
save_plot(args$figure, gInt, base_height = 10, base_width = 10)
