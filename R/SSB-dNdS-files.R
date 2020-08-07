library(tidyverse)
library(readxl)
library(cowplot)
library(dndscv)
library(GenomicRanges)
library(readr)

message("Read in meta data for the oesophagus")
dfdonor <- read_xlsx(snakemake@input$oesophaguspatientinfo, skip = 1) %>%
  dplyr::rename(patient = PD) %>%
  filter(patient == snakemake@wildcards$oes_sample)

message("Read in mutation data for the oesophagus")
df <- read_csv(snakemake@input$oesophagusdata) %>%
  mutate(sumvaf = sumvaf * 2) %>%
  dplyr::rename(patient = donor) %>%
  filter(patient == snakemake@wildcards$oes_sample)

message("Create vector of intervals for i-dN/dS")
minarea <- snakemake@params$minarea
maxarea <- snakemake@params$maxarea
step <- snakemake@params$step
areacutoff <- seq(minarea, maxarea, step)

for (i in 1:length(areacutoff)){
    message(paste0("Completed ", i , "/", length(areacutoff)))
    df %>%
      filter(sumvaf < areacutoff[i]) %>%
      write_delim(paste0(snakemake@output[[1]], snakemake@wildcards$oes_sample, "_", i, ".txt"), delim = "\t")
}
