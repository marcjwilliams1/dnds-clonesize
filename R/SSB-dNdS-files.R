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
parser$add_argument('--outputdir', type = 'character',
                    help="output directory")
parser$add_argument('--sample', type = 'character',
                    help="patient sample")
parser$add_argument('--step', type = 'double',
                    help="stepsize for interval dN/dS")
parser$add_argument('--minarea', type = 'double',
                    help="Min area for interval dN/dS")
parser$add_argument('--maxarea', type = 'double',
                    help="Min area for interval dN/dS")
args <- parser$parse_args()

message("Read in meta data for the oesophagus")
dfdonor <- read_xlsx(args$patientinfo, skip = 1) %>%
  dplyr::rename(patient = PD) %>%
  filter(patient == args$sample)

message("Read in mutation data for the oesophagus")
df <- read_csv(args$oesophagusdata) %>%
  mutate(sumvaf = sumvaf * 2) %>%
  dplyr::rename(patient = donor) %>%
  filter(patient == args$sample)

message("Create vector of intervals for i-dN/dS")
minarea <- args$minarea
maxarea <- args$maxarea
step <- args$step
areacutoff <- seq(minarea, maxarea, step)

for (i in 1:length(areacutoff)){
    message(paste0("Completed ", i , "/", length(areacutoff)))
    df %>%
      filter(sumvaf < areacutoff[i]) %>%
      write_delim(paste0(args$outputdir, args$sample, "_", i, ".txt"), delim = "\t")
}
