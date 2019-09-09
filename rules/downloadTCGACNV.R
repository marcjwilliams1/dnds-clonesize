library(TCGAbiolinks)
library(stringr)
library(tidyverse)
library(readr)
library(argparse)

parser <- ArgumentParser(description = "Generate Final Figures")
parser$add_argument('--CNVfile', type='character',
                    help="TCGA CNV file hg38")
args <- parser$parse_args()

#get all TCGA cancer subtype codes
tcgacodes <- TCGAbiolinks:::getGDCprojects()$project_id
tcgacodes <- unlist(lapply(tcgacodes[str_detect(tcgacodes, "TCGA")],
                    function(x){strsplit(x, "-")[[1]][2]}))

message("Download copy number data...")
tcgacodes <- TCGAbiolinks:::getGDCprojects()$project_id
query <- GDCquery(project = tcgacodes[str_detect(tcgacodes, "TCGA")],
                data.category = "Copy Number Variation",
                data.type = "Copy Number Segment")

cnvsamples <- getResults(query)
GDCdownload(query)
dfhg38cnv <- GDCprepare(query)
message("Writing CNA to file...")
write_delim(dfhg38cnv, args$CNVfile, delim = ",")
