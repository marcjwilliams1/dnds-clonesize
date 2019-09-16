library(tidyverse)
library(readxl)
library(cowplot)
library(dndscv)
library(GenomicRanges)
library(readr)
library(argparse)

parser <- ArgumentParser(description = "Calculate dN/dS normal")
parser$add_argument('--patientinfo', type='character',
                    help="Patient info xlsx file")
parser$add_argument('--oesophagusdata', type='character',
                    help="Oesophagus clone size data")
parser$add_argument('--hotspots', type='character',
                    help="Hotspot data")
parser$add_argument('--oesophagusdnds', type='character',
                    help="Oesophagus dnds")
parser$add_argument('--oesophagusdndsgenes', type='character',
                    help="Oesophagus dnds for neutral genes")
parser$add_argument('--step', type = 'double',
                    help="stepsize for interval dN/dS")
parser$add_argument('--minarea', type = 'double',
                    help="Min area for interval dN/dS")
parser$add_argument('--maxarea', type = 'double',
                    help="Min area for interval dN/dS")
args <- parser$parse_args()

args$step <- 1 * args$step

message("Read in meta data for the oesophagus")
dfdonor <- read_xlsx(args$patientinfo, skip = 1) %>%
  dplyr::rename(patient = PD)

message("Read in mutation data for the oesophagus")
df <- read_csv(args$oesophagusdata) %>%
  mutate(sumvaf = sumvaf * 2)

message("Create vector of intervals for i-dN/dS")
minarea <- args$minarea
maxarea <- args$maxarea
step <- args$step
areacutoff <- seq(minarea, maxarea, step)

message("Find unique genes")
target_genes <- unique(df$gene)

dfdnds.patient <- data.frame()
dfdnds.genes.patient <- data.frame()
df.hotspots <- data.frame()
i <- 1

message("Calculate global and per gene i-dN/dS per patient")
for (p in unique(df$donor)){
  message(paste0("Analysing patient ", p, "(patient ", i, "/", length(unique(df$donor)), " )"))
  for (cutoff in areacutoff){
    message(paste0("Max area is ", cutoff))
    x1 <- df %>%
      filter(donor == p) %>%
      filter(sumvaf < cutoff) #filter for mutations with sumvaf < cutoff
    x <- dndscv(x1, gene_list = target_genes,
                outp = 3, max_muts_per_gene_per_sample = Inf,
                max_coding_muts_per_sample = Inf,
                outmats = T)
    out <- x$globaldnds %>%
        mutate(areacutoff = cutoff, nmutations = length(x1$donor), patient = p)
    dfdnds.patient <- rbind(dfdnds.patient, out)
    combined <- left_join(x$sel_cv, x$sel_loc, by = c("gene_name", "n_syn", "n_mis", "n_non", "n_spl"))
    out2 <- combined %>%
      mutate(areacutoff = cutoff, patient = p)
    dfdnds.genes.patient <- rbind(dfdnds.genes.patient, out2)
    hotspots_res <- sitednds(x)
    if (is.null(hotspots_res$recursites)){
        message("No recurrent mutations...")
        next
    } else {
        message(paste0("Number of reccurent mutations: ", length(hotspots_res$recursites$chr)))
        hotspots <- hotspots_res$recursites %>%
          mutate(areacutoff = cutoff, patient = p)
        df.hotspots <- bind_rows(df.hotspots, hotspots)
    }
  }
  message("\n")
  i <- i + 1
}

message("Write output to file")
write_csv(dfdnds.patient, args$oesophagusdnds)
write_csv(dfdnds.genes.patient, args$oesophagusdndsgenes)
write_csv(dfhotspots, args$hotspots)
