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
parser$add_argument('--skindata', type='character',
                    help="Skin clone size data")
parser$add_argument('--oesophagusdnds', type='character',
                    help="Oesophagus dnds")
parser$add_argument('--oesophagusdndsgenes', type='character',
                    help="Oesophagus dnds for neutral genes")
parser$add_argument('--skindnds', type='character',
                    help="Skin dnds")
parser$add_argument('--skindndsgenes', type='character',
                    help="Skin dnds for neutral genes")
parser$add_argument('--oesophagusdndsneutral', type='character',
                    help="Oesophagus dnds")
parser$add_argument('--oesophagusdndsgenesneutral', type='character',
                    help="Oesophagus dnds per gene for neutral genes")
parser$add_argument('--singlepatient', type = 'character',
                    help="Single patient ID")
parser$add_argument('--singlepatientdnds', type = 'character',
                    help="Output file for dnds in a single patient")
parser$add_argument('--singlepatientdndsgenes', type = 'character',
                    help="Output file for dnds per gene in a single patient")
parser$add_argument('--step', type = 'double',
                    help="stepsize for interval dN/dS")
parser$add_argument('--minarea', type = 'double',
                    help="Min area for interval dN/dS")
parser$add_argument('--maxarea', type = 'double',
                    help="Min area for interval dN/dS")
args <- parser$parse_args()

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
                max_coding_muts_per_sample = Inf)
    out <- x$globaldnds %>%
        mutate(areacutoff = cutoff, nmutations = length(x1$donor), patient = p)
    dfdnds.patient <- rbind(dfdnds.patient, out)
    combined <- left_join(x$sel_cv, x$sel_loc, by = c("gene_name", "n_syn", "n_mis", "n_non", "n_spl"))
    out2 <- combined %>%
      mutate(areacutoff = cutoff, patient = p)
    dfdnds.genes.patient <- rbind(dfdnds.genes.patient, out2)
  }
  message("\n")
  i <- i + 1
}

message("Write output to file")
write_csv(dfdnds.patient, args$oesophagusdnds)
write_csv(dfdnds.genes.patient, args$oesophagusdndsgenes)

message("Make vector of genes that show no evidence of selection")
target_genesneutral <- target_genes[!target_genes %in%
                    c("NOTCH1", "NOTCH2", "NOTCH3", "TP53","NFE2L2", "ARID1A",
                    "CDH1",  "CREBBP", "FAT1", "SALL1", "KMT2D",
                    "ARID2","CCND2", "CCND1", "KMT2A", "SCN1A", "FBXW7",
                      "KMT2C", "PIK3CA", "SPHKAP", "CUL3", "TP63", "AJUBA",
                      "EPHA2", "IRF6", "FGFR2")]


minarea <- args$minarea
maxarea <- args$maxarea
step <- args$step
areacutoff <- seq(minarea, maxarea, step)
target_genes <- unique(df$gene)

dfdnds <- data.frame()
dfdnds.genes <- data.frame()

message("Calculate global i-dN/dS for all putatively neutral genes")
for (cutoff in areacutoff){
    message(paste0("Max area is ", cutoff))
    x1 <- df %>%
      filter(sumvaf < cutoff) #filter for mutations with ccf < cutoff
    x <- dndscv(x1, gene_list = target_genesneutral,
                outp = 3, max_muts_per_gene_per_sample = Inf,
                max_coding_muts_per_sample = Inf)
    out <- x$globaldnds %>%
      mutate(areacutoff = cutoff, nmutations = length(x1$donor))
    dfdnds <- rbind(dfdnds, out)
    combined <- left_join(x$sel_cv, x$sel_loc, by = c("gene_name", "n_syn", "n_mis", "n_non", "n_spl"))
    out2 <- combined %>%
      mutate(areacutoff = cutoff)
    dfdnds.genes <- rbind(dfdnds.genes, out2)
}

message("Write i-dN/dS for neutral genes to file.")
write_csv(dfdnds, args$oesophagusdndsneutral)
write_csv(dfdnds.genes, args$oesophagusdndsgenesneutral)

message("Filter for single patient")
df1 <- df %>%
   filter(donor == args$singlepatient) %>%
   mutate(sumvaf = sumvaf * 2) %>% #times by 2 to convert to area
   mutate(x = cut(sumvaf, breaks = c(0.004 * 2, 0.01 * 2, 0.03 * 2, 0.07 * 2, 0.2 * 2, 6 * 2))) %>%
   group_by(x) %>%
   mutate(medianvaf = median(sumvaf)) %>%
   ungroup()

target_genes <- unique(df1$gene)

dfdnds <- data.frame()
dfdnds.genes <- data.frame()
bins <- sort(unique(df1$x))
j <- 1

message("Calculate dN/dS in bins for single patient")
for (f in sort(unique(df1$medianvaf))){
  x1 <- df1 %>%
    filter(medianvaf ==f) #filter for mutations with ccf < cutoff
  x <- dndscv(x1, gene_list = target_genes,
              outp = 3, max_muts_per_gene_per_sample = Inf,
              max_coding_muts_per_sample = Inf)
  out <- x$globaldnds %>%
    mutate(medianbin = f, vafbin = bins[j], nmutations = length(x1$donor))
  dfdnds <- rbind(dfdnds, out)
  combined <- left_join(x$sel_cv, x$sel_loc, by = c("gene_name", "n_syn", "n_mis", "n_non", "n_spl"))
  out2 <- combined %>%
    mutate(medianbin = f, vafbin = bins[j])
  dfdnds.genes <- rbind(dfdnds.genes, out2)
  j <- j + 1
}

message("Write single patient to file")
write_csv(dfdnds, args$singlepatientdnds)
write_csv(dfdnds.genes, args$singlepatientdndsgenes)


################################################
# Analyse skin data
################################################

message("Read in skin data")
df <- read_csv(args$skindata)
data("dataset_normalskin_genes", package="dndscv")
target_genes

minarea <- args$minarea + 0.01
maxarea <- args$maxarea
step <- args$step
areacutoff <- seq(minarea, maxarea, step)

dfdnds.patient <- data.frame()
dfdnds.genes.patient <- data.frame()
i <- 1
for (p in unique(df$patient)){
  message(paste0("Analysing patient ", p, "(patient ", i, "/", length(unique(df$patient)), " )"))
    for (cutoff in areacutoff){
        message(paste0("Max area is ", cutoff))
        x1 <- df %>%
          filter(patient == p) %>%
          filter(clone.area < cutoff) #filter for mutations with ccf < cutoff
        x <- dndscv(x1, gene_list = target_genes,
                    outp = 3, max_muts_per_gene_per_sample = Inf,
                    max_coding_muts_per_sample = Inf)
        out <- x$globaldnds %>%
            mutate(areacutoff = cutoff, nmutations = length(x1$donor), patient = p)
        dfdnds.patient <- rbind(dfdnds.patient, out)
        combined <- left_join(x$sel_cv, x$sel_loc, by = c("gene_name", "n_syn", "n_mis", "n_non", "n_spl"))
        out2 <- combined %>%
          mutate(areacutoff = cutoff, patient = p)
        dfdnds.genes.patient <- bind_rows(dfdnds.genes.patient, out2)
  }
  message("\n")
  i <- i + 1
}

message("Write skin i-dN/dS values to file")
write_csv(dfdnds.patient, args$skindnds)
write_csv(dfdnds.genes.patient, args$skindndsgenes)
