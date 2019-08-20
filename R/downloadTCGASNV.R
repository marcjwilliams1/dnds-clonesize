library(TCGAbiolinks)
library(stringr)
library(tidyverse)
library(readr)
library(argparse)

parser <- ArgumentParser(description = "Generate Final Figures")
parser$add_argument('--MAFfile', type='character',
                    help="TCGA MAF file hg38")
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
message("Writing MAF to file...")
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
