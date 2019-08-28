library(argparse)
library(tidyverse)
library(readxl)

parser <- ArgumentParser(description = "Format SSB file to tidy format")
parser$add_argument('--inputfile', type='character',
                    help="Input file")
parser$add_argument('--outputfile', type='character',
                    help="Output file")
parser$add_argument('--step', type = 'double',
                    help="stepsize for interval dN/dS")
parser$add_argument('--minarea', type = 'double',
                    help="Min area for interval dN/dS")
parser$add_argument('--maxarea', type = 'double',
                    help="Min area for interval dN/dS")
args <- parser$parse_args()

dfall <- read_xlsx(args$inputfile, sheet = 1) %>%
    mutate(gene = "All")
dftp53 <- read_xlsx(args$inputfile, sheet = 2) %>%
    mutate(gene = "TP53")
dfnotch <- read_xlsx(args$inputfile, sheet = 3) %>%
    mutate(gene = "NOTCH1")

df <- bind_rows(dfall, dftp53) %>%
    bind_rows(., dfnotch)

df <- df %>%
    separate(ID, c("ID", "idx"), "_") %>%
    mutate(idx = as.numeric(idx)) %>%
    arrange(gene, idx)

message("Create vector of intervals for i-dN/dS")
minarea <- args$minarea
maxarea <- args$maxarea
step <- args$step
areacutoff <- seq(minarea, maxarea, step)
dfintervals <- data.frame(cutoff = areacutoff) %>%
    mutate(idx = 1:n())

df <- left_join(df, dfintervals)
write_csv(df, args$outputfile)
