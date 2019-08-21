#load required libraries
library(MASS)
library(tidyverse)
library(dndscv)
library(GenomicRanges)
library(cowplot)
library(ggthemes)
library(broom)
library(viridis)
library(Hmisc)
library(jcolors)
library(argparse)

parser <- ArgumentParser(description = "Calculate dN/dS normal")
parser$add_argument('--data', type='character',
                    help="TCGA data file")
parser$add_argument('--drivergenelist', type='character',
                    help="driver gene list")
parser$add_argument('--essentialgenelist', type='character',
                    help="essential gene list")
parser$add_argument('--allgenes', type='character',
                    help="All genes")
parser$add_argument('--dndscvref', type='character',
                    help="dndscv Reference file")
parser$add_argument('--baseline', type='character',
                    help="baseline dn/ds values")
parser$add_argument('--baseline_validation', type='character',
                    help="validation for baseline dn/ds values")
parser$add_argument('--vafclonality', type='character',
                    help="Mutation VAFs with clonality assignment")
parser$add_argument('--dndsclonality', type='character',
                    help="Mutation VAFs with clonality assignment")
parser$add_argument('--dndsclonality_percancertype', type='character',
                    help="dn/ds as a function of clonality per cancertype")
parser$add_argument('--nmutations_gene', type='character',
                    help="Number of mutations per gene")
parser$add_argument('--nmutations_gene_percancertype', type='character',
                    help="Number of mutations per gene per cancertype")
parser$add_argument('--intervaldnds', type='character',
                    help="Interval dN/dS output")
args <- parser$parse_args()

###############################
# Read in data
###############################
message("Read in data")

df <- read_csv(args$data);

#read in list of drivers
drivers <- read.table(args$drivergenelist, header = FALSE) %>%
    filter(V1 != "FAM46C") #this gene is not compatible with dndscv

cellessential <- read.csv(args$essentialgenelist, header = T) %>%  #cell essential genes from Blomen et al
    filter(!(V1 %in% drivers$V1)) %>%
    filter(!(V1 %in% Cs(PET112, IKBKAP, UQCC, SLMO2, SLMO2, ATP5B, PPP2R4,
                      TOMM70A, ATP5C1, WIBG, C9orf114, C9orf114, ATP5A1,
                      WBSCR22, ATP5O, SHFM1, GTPBP5, GNB2L1, TOMM70A, SRPR,
                      CIRH1A, NHP2L1, NUPL1, UFD1L, ZNF259, FAM96B, DIEXF,
                      SKIV2L2, CCDC94, UTP11L, MKI67IP, NARFL, GLTSCR2, FDX1L))) #genes not compatible with dndscv

randomgenes <- read.table(args$allgenes, header = TRUE) %>%
    filter(!V1 %in% cellessential$V1, !V1 %in% drivers$V1) %>%
    sample_n(1000)

dfunfilt <- df

df <- df %>%
    mutate(VAFcorr = VAF/cellularity) %>%
    mutate(absCNcorrected = ifelse(absCNcorrected > 5, 5, absCNcorrected)) %>%
    filter(ploidy < 2.5) %>% #remove WGD samples where MCN is likely to be incorrect
    filter(!is.na(cellularity)) %>% #remove samples with no cellularity estimate
    #filter(cellularity < 1.0) %>% #remove samples cellularity estimate = 1 as this is likely incorrect
    filter(mutation_type == "SNV") %>% #only include SNVs
    group_by(sampleid) %>%
    mutate(sampledepth = mean(t_depth)) %>%
    ungroup() %>%
    mutate(cancertype = ifelse(cancertype == "READ", "COAD", cancertype)) %>% #combine COAD and READ
    group_by(sampleid) %>%
    mutate(effectivedepth = cellularity * sampledepth) %>%
    mutate(meaneffdepth = mean(effectivedepth),
          subclonalcutoff = 1.0 - (3 * (sqrt(1/meaneffdepth)))) %>%
    ungroup() %>%
    mutate(clonality = case_when(
                    MCN < 0.5 ~ "Subclonal",
                    MCN > 0.5 & MCN < 1.5 ~ "Clonal",
                    MCN > 1.5 ~ "Amplified"))

# print out number of samples per cancer type
df %>%
    filter(effectivedepth > 50.0) %>%
    distinct(sampleid, cancertype) %>%
    group_by(cancertype) %>%
    summarise(n = n()) %>% print()

###############################
# Calculate background dN/dS
###############################

dfback <- df %>%
    select(sampleid, chr, start, Reference_Allele, Tumor_Seq_Allele2, everything())

library(dndscv)
library(Hmisc)
dndsrandom <- data.frame()

for (i in 1:50){
    print(i)
    randomgenestemp <- read.table(args$allgenes, header = TRUE) %>%
        filter(!V1 %in% cellessential$V1, !V1 %in% drivers$V1) %>%
        sample_n(1000)
    x <- dndscv(dfback, gene_list = randomgenestemp$V1, outp = 1,
                refdb = args$dndscvref)
    out <- x$globaldnds %>% mutate(mutationtype = "randomgenes")
    dndsrandom <- rbind(dndsrandom, out)
}

baseline <- dndsrandom %>%
    group_by(name) %>%
    summarise(dndslq_bl = quantile(mle, 0.025), dnds_bl = mean(mle), dndsuq_bl = quantile(mle, 0.975))
write_csv(baseline, args$baseline)
baseline

# generate an additional 100 random samples and ensure mean is 1
dndsrandom2 <- data.frame()

for (i in 1:50){
    print(i)
    randomgenestemp <- read.table(args$allgenes, header = TRUE) %>%
        filter(!V1 %in% cellessential$V1, !V1 %in% drivers$V1) %>%
        sample_n(1000)
    x <- dndscv(dfback, gene_list = randomgenestemp$V1, outp = 1,
                refdb = args$dndscvref)
    out <- x$globaldnds %>% mutate(mutationtype = "randomgenes")
    dndsrandom2 <- rbind(dndsrandom2, out)
}

dndsrandom2 %>%
        left_join(., baseline, by = "name") %>%
        mutate(mle = mle - (dnds_bl - 1),
              cilow = cilow - (dnds_bl - 1),
              cihigh = cihigh - (dnds_bl - 1)) %>%
        write_csv(., args$baseline_validation)

message("Calculate clonal vs subclonal dN/dS")
dfsubc <- df %>%
    filter(meaneffdepth > 50.0) %>% #remove samples with effective depth < 50X
    select(sampleid, chr, start, Reference_Allele, Tumor_Seq_Allele2, everything())

print(length(unique(dfsubc$sampleid)))

write_csv(dfsubc, args$vafclonality)

dfdnds <- data.frame()

for (clon in c("Subclonal", "Clonal", "Amplified")){
    print(clon)
    x1 <- dfsubc %>% filter(clonality == clon)

    x <- dndscv(x1 , gene_list = as.character(drivers$V1), outp = 1,
               refdb = args$dndscvref)
    out <- x$globaldnds %>% mutate(clonality = clon, mutationtype = "drivers")
    dfdnds <- rbind(dfdnds, out)
    print(out)

    x <- dndscv(x1, gene_list = as.character(cellessential$V1), outp = 1,
               refdb = args$dndscvref)
    out <- x$globaldnds %>% mutate(clonality = clon, mutationtype = "cellessential")
    dfdnds <- rbind(dfdnds, out)

    x <- dndscv(x1, outp = 1,
               refdb = args$dndscvref)
    out <- x$globaldnds %>% mutate(clonality = clon, mutationtype = "all")
    dfdnds <- rbind(dfdnds, out)
}

dfdnds <- dfdnds %>%
    mutate(clonality = factor(clonality,
                              levels = c("Subclonal", "Clonal", "Amplified")))
write_csv(dfdnds, args$dndsclonality)


message('Calculate i-dN/dS')
dfsubc <- df %>%
    filter(meaneffdepth > 50.0) %>%
    select(sampleid, chr, start, Reference_Allele, Tumor_Seq_Allele2, everything())

dfdnds.interval <- data.frame()
cutoff <- seq(0.2, 0.5, 0.025)

for (i in cutoff){
    x1 <- dfsubc %>% filter(MCN < i)

    x <- dndscv(x1 , gene_list = as.character(drivers$V1), outp = 1,
               refdb = args$dndscvref)
    out <- x$globaldnds %>% mutate(MCN = i, mutationtype = "drivers")
    dfdnds.interval <- rbind(dfdnds.interval, out)
    print(out)

    x <- dndscv(x1, outp = 1,
               refdb = args$dndscvref)
    out <- x$globaldnds %>% mutate(MCN = i, mutationtype = "all")
    dfdnds.interval <- rbind(dfdnds.interval, out)
}
write.csv(dfdnds.interval, args$intervaldnds)

message("dN/dS clonality per cancer type")

dfsubc <- df %>%
    filter(meaneffdepth > 50.0) %>% #remove samples with effective depth < 50X
    select(sampleid, chr, start, Reference_Allele, Tumor_Seq_Allele2, everything())

cancertypes <- dfsubc %>%
    distinct(sampleid, cancertype) %>%
    group_by(cancertype) %>%
    summarise(n = n()) %>%
    filter(n >= 100) %>%
    pull(cancertype)

dfdnds.cancertype <- data.frame()

for (ct in cancertypes){
    print(ct)
    for (clon in c("Subclonal", "Clonal", "Amplified")){
        x1 <- dfsubc %>% filter(clonality == clon, cancertype == ct)

        x <- dndscv(x1 , gene_list = as.character(drivers$V1), outp = 1,
                   refdb = args$dndscvref)
        out <- x$globaldnds %>% mutate(clonality = clon, mutationtype = "drivers", cancertype = ct)
        dfdnds.cancertype <- rbind(dfdnds.cancertype, out)

        x <- dndscv(x1, outp = 1,
                   refdb = args$dndscvref)
        out <- x$globaldnds %>% mutate(clonality = clon, mutationtype = "all", cancertype = ct)
        dfdnds.cancertype <- rbind(dfdnds.cancertype, out)
    }
}

dfdnds.cancertype <- dfdnds.cancertype %>%
    mutate(clonality = factor(clonality,
                              levels = c("Subclonal", "Clonal", "Amplified")))
write_csv(dfdnds.cancertype, args$dndsclonality_percancertype)


message("Get number of mutations per gene")

dftemp <- df %>%
    select(sampleid, chr, start, Reference_Allele, Tumor_Seq_Allele2, everything())

#we'll use dndscv to annotate the mutations to be consistent
x <- dndscv(dftemp, outp = 1,
                refdb = args$dndscvref)

dfannot <- x$annotmuts
df2 <- df %>%
    dplyr::rename(pos = start, ref = Reference_Allele, sampleID = sampleid, mut = Tumor_Seq_Allele2)

df2 <- left_join(df2, dfannot, by = c("sampleID", "chr", "pos", "ref", "mut"))

df2 %>%
    filter(gene %in% drivers$V1, clonality == "Subclonal", MCN > 0.2)  %>%
    group_by(gene, impact) %>%
    summarise(n = n()) %>%
    write_csv(., args$nmutations_gene)

df2 %>%
    filter(gene %in% drivers$V1, clonality == "Subclonal", MCN > 0.2)  %>%
    group_by(gene, impact, cancertype) %>%
    summarise(n = n()) %>%
    write_csv(., args$nmutations_gene_percancertype)
