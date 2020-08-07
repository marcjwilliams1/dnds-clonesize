library(tidyverse)
library(readxl)
library(cowplot)
library(dndscv)
library(GenomicRanges)
library(readr)

snakemake@params$step <- 1 * snakemake@params$step
snakemake@params$minarea <- 3 * snakemake@params$minarea

message("Read in meta data for the oesophagus")
dfdonor <- read_xlsx(snakemake@input$oesophaguspatientinfo, skip = 1) %>%
  dplyr::rename(patient = PD)

message("Read in mutation data for the oesophagus")
df <- read_csv(snakemake@input$oesophagusdata) %>%
  mutate(sumvaf = sumvaf * 2)

message("Create vector of intervals for i-dN/dS")
minarea <- snakemake@params$minarea
maxarea <- snakemake@params$maxarea
step <- snakemake@params$step
areacutoff <- seq(minarea, maxarea, step)

message("Find unique genes")
target_genes <- unique(df$gene)

dfdnds.patient <- data.frame()
dfdnds.genes.patient <- data.frame()
df.hotspots <- data.frame()
i <- 1

message("Calculate global and per gene i-dN/dS per patient")
for (p in rev(unique(df$donor))){
  message(paste0("Analysing patient ", p, "(patient ", i, "/", length(unique(df$donor)), " )"))
  for (cutoff in areacutoff){
    message(paste0("Max area is ", cutoff))
    x1 <- df %>%
      filter(donor == p) %>%
      filter(sumvaf < cutoff) #filter for mutations with sumvaf < cutoff
    print(head(x1))
    x <- dndscv(x1, gene_list = target_genes,
                outp = 3, max_muts_per_gene_per_sample = Inf,
                max_coding_muts_per_sample = Inf,
                outmats = T,
                uniquemuts = F)
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
          mutate(areacutoff = cutoff, patient = p,
                 nmuts = length(hotspots_res$recursites$chr))
        hotspots$chr <- as.character(hotspots$chr)
        df.hotspots <- bind_rows(df.hotspots, hotspots)
    }
  }
  message("\n")
  i <- i + 1
}

message("Write output to file")
write_csv(dfdnds.patient, snakemake@output$oesophagusdnds)
write_csv(dfdnds.genes.patient, snakemake@output$oesophagusdndsgenes)
write_csv(df.hotspots, snakemake@output$hotspots)
