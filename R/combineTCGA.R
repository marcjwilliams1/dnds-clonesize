library(TCGAbiolinks)
library(stringr)
library(tidyverse)
library(readr)
library(argparse)

parser <- ArgumentParser(description = "Generate Final Figures")
parser$add_argument('--MAFfile', type='character',
                    help="TCGA MAF file hg38")
parser$add_argument('--CNVfile', type='character',
                    help="TCGA CNV file hg38")
parser$add_argument('--ascatcellularity', type='character',
                    help="Ascat cellularity file")
parser$add_argument('--outputfile', type='character',
                    help="Ascat cellularity file")
args <- parser$parse_args()

message("Read in MAF to file...")
dfhg38 <- read_delim(args$MAFfile, delim = ",")

message("Select and rename some columns...")
dfhg38 <- dfhg38 %>%
    dplyr::mutate(VAF = t_alt_count/t_depth) %>%
    dplyr::rename(chr = Chromosome, start = Start_Position, end = End_Position) %>%
    dplyr::select(sampleid, chr, start, end, Reference_Allele, Tumor_Seq_Allele2, VAF,
                  t_depth, t_ref_count, t_alt_count, n_depth,
                  n_ref_count, n_alt_count, cancertype) %>%
    dplyr::mutate(nref = str_length(Reference_Allele), nalt = str_length(Tumor_Seq_Allele2)) %>%
    dplyr::mutate(mutation_type = ifelse((nref - nalt) == 0, "SNV", "INS/DEL")) %>%
    dplyr::select(-nref, -nalt)

message("Summarise the number of sample per tumour type")
dfhg38 %>%
    distinct(sampleid, cancertype) %>%
    group_by(cancertype) %>%
    summarise(n = n())

message("Read in CNV file")
dfhg38cnv <- read_delim(args$CNVfile, delim = ",")

message("Remove some columns and filter sampleid")
dfhg38cnv <- dfhg38cnv %>%
    dplyr::mutate(sampleid = str_sub(Sample, 1, 16)) %>%
    dplyr::select(-Sample) %>%
    filter(sampleid %in% dfhg38$sampleid)

message("Read in cellularity")
cellularity <- read.delim(args$ascatcellularity, header = T) %>%
  mutate(sampleid = str_sub(gsub("[.]", "-", Sample), 1, 16)) %>%
  select(-Sample) %>%
  dplyr::rename(ploidy = Ploidy, cellularity = Aberrant_Cell_Fraction.Purity.)

message("Combine all file types")
dfsnv <- dfhg38

#dfsnv <- dfhg38
df1temp <- left_join(dfsnv, cellularity, by = c("sampleid")) %>%
    filter(cellularity > 0.2)

dfhg38cnvt <- dfhg38cnv %>%
    dplyr::rename(chr = Chromosome) %>%#, start = Start, end = End) %>%
    dplyr::mutate(chr = paste0("chr", chr)) %>%
    dplyr::filter(sampleid %in% unique(df1temp$sampleid))

#join snv, cnv and cellularity
df2temp <- inner_join(df1temp, dfhg38cnvt, by = c("sampleid", "chr")) %>%
    filter(start >= Start & end <= End) %>%
    select(-Start, -End)

dfcombinedhg38 <- df2temp %>%
    select(-n_ref_count, -n_alt_count, -Num_Probes)

dfout <- dfcombinedhg38 %>%
    mutate(cellularity = ifelse(is.na(cellularity), 1, as.numeric(cellularity))) %>%
    #calculate CN by correcting for cellularit
    mutate(CN = 2^Segment_Mean * 2, CNcorrected = (2^Segment_Mean + cellularity - 1) * (2 / cellularity),
          absCN = round(CN), absCNcorrected = round(CNcorrected)) %>%
    #don't allow CN == 0
    mutate(absCN = ifelse(absCN == 0, 1, absCN), absCNcorrected = ifelse(absCNcorrected == 0, 1, absCNcorrected)) %>%
    mutate(MCN = ((CNcorrected - 2) * 1 + 2) * VAF/cellularity)

write_delim(dfout, args$outputfile, delim = ",")
