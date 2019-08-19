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

#get all TCGA cancer subtype codes
tcgacodes <- TCGAbiolinks:::getGDCprojects()$project_id
tcgacodes <- unlist(lapply(tcgacodes[str_detect(tcgacodes, "TCGA")],
                    function(x){strsplit(x, "-")[[1]][2]}))

message("Download mutect MAF files for all samples and save to csv file")
dfhg38 <- data.frame()
for (t in tcgacodes){
  maf <- GDCquery_Maf(t, pipelines = "mutect") %>%
    mutate(cancertype = t)
  dfhg38 <- rbind(dfhg38, maf)
}
dfhg38 <- dfhg38 %>%
    mutate(sampleid = str_sub(Tumor_Sample_Barcode, 1, 16)) #annotate so barcode is consistent with TCGA CNV id
write_delim(dfhg38, args$MAFfile, delim = ",")

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

message("Download copy number data...")
tcgacodes <- TCGAbiolinks:::getGDCprojects()$project_id
query <- GDCquery(project = tcgacodes[str_detect(tcgacodes, "TCGA")],
                data.category = "Copy Number Variation",
                data.type = "Copy Number Segment")

cnvsamples <- getResults(query)
GDCdownload(query)
dfhg38cnv <- GDCprepare(query)
write_delim(dfhg38cnv, args$CNVfile, delim = ",")

message("Remove some columns and filter sampleid")
dfhg38cnv <- dfhg38cnv %>%
    dplyr::mutate(sampleid = str_sub(Sample, 1, 16)) %>%
    dplyr::select(-Sample, -X1) %>%
    filter(sampleid %in% dfhg38$sampleid)

message("Read in cellularity")
cellularity <- read.delim(args$ascatcellularity, header = T) %>%
  mutate(sampleid = str_sub(gsub("[.]", "-", Sample), 1, 16)) %>%
  select(-Sample) %>%
  dplyr::rename(ploidy = Ploidy, cellularity = Aberrant_Cell_Fraction.Purity.)

 message("Combine all file types")
 dfsnv <- dfhg38 %>%
    dplyr::mutate(VAF = t_alt_count/t_depth) %>%
    dplyr::rename(chr = Chromosome, start = Start_Position, end = End_Position) %>%
    dplyr::select(sampleid, chr, start, end, Reference_Allele, Tumor_Seq_Allele2, VAF,
                  t_depth, t_ref_count, t_alt_count, n_depth,
                  n_ref_count, n_alt_count, cancertype) %>%
    dplyr::mutate(nref = str_length(Reference_Allele), nalt = str_length(Tumor_Seq_Allele2)) %>%
    dplyr::mutate(mutation_type = ifelse((nref - nalt) == 0, "SNV", "INS/DEL")) %>%
    dplyr::select(-nref, -nalt)

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
