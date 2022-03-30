### This script will provide a full reverse
### pipeline from a A. thaliana gene ID to
### plotting the methylation state of its
### overlapping DMRs

### The pipeline will have 3 steps:
### 1) Get gene  number from ATID
### 2) Get DMRs from gene number
### 3) Plot DMRs using necessary data

### The inputs for this script are many:
### 1) output from ARPEGGIO about all 
###    DMRs overlapping with genes
### 2) reciprocal best hit file with matches
###    between our species and thaliana
### 3) Rdata file with candidate regions
###    showing differential methylation
### 4-9) Bismark cov files of the samples of 
###    interest (first diploid then polyploid)
### 10) ATID of interest
### 10) output folder

library(dmrseq)
library(tidyverse)
library(GenomicRanges)
library(annotatr)

### Set up command line arguments

comm_args <- commandArgs(trailingOnly = TRUE)

### Fist argument: ARPEGGIO's DM genes
dm_genes <- normalizePath(comm_args[1])

### Second argument: BLAST reciprocal best hit
blast_path <- normalizePath(comm_args[2])

### Third argument: R data from dmrseq
Rdata <- normalizePath(comm_args[3])

### Fourth to sixth argument: bismark cov set 1 (diploid)

cov_1_diploid <- comm_args[4]
cov_2_diploid <- comm_args[5]
cov_3_diploid <- comm_args[6]

### Seventh to ninth argument: bismark cov set 2 (polyploid)

cov_1_polyploid <- comm_args[7]
cov_2_polyploid <- comm_args[8]
cov_3_polyploid <- comm_args[9]

### Tenth argument A. thaliana geneID of interest
ATID <- comm_args[10]

### Eleventh argument: output folder
output_path <- comm_args[11]

### Twelveth argument: BED annotation
bed_anno <- comm_args[12]

### Thirteenth argument: context
context <- comm_args[13]

### set working directory

#setwd(normalizePath(output_path))

### Read files

dmr <- read.delim(dm_genes)
blast <- read.delim(blast_path)

# Load Rdata file

load(Rdata)

print("All files read and loaded!")

### reshape blast hit file to remove "scaffold" prefix in the first column

gene_list <- str_split_fixed(blast$Ahal, "\\.", 2)[,2]
if (length(gene_list)==0){
  gene_list <- str_split_fixed(blast$Alyr, "\\.", 2)[,2]
}

TAIR10_cds <- str_split_fixed(blast$TAIR10_cds, "\\.", 2)[,1]

blast_new <- data.frame(geneID = gene_list,
                        thalianaID = TAIR10_cds)

### match and pick corresponding thaliana ID

matches <- match(dmr$geneID, blast_new$geneID)

dmr_new <- dplyr::mutate(dmr, thalID=blast_new[matches, 2])

gene_of_interest <- dplyr::filter(dmr_new, thalID==ATID)

### If gene is not present, stop the script

if(nrow(gene_of_interest) == 0){
  stop("The gene of interest was not found")
}

print("The gene was found!")

### We filter the regions file to get the 
### corresponding DMR

if(nrow(gene_of_interest) > 1){
  thalIDregion <- regions
  for(i in 1:nrow(gene_of_interest)){
    regions_scaffold <- regions[regions@seqnames == gene_of_interest$seqname[i]]
    thalIDregion[i] <- regions_scaffold[regions_scaffold@ranges@start == gene_of_interest$region_start[i]]
  }
  thalIDregion <- thalIDregion[1:nrow(gene_of_interest)]
  
} else {
regions_scaffold <- regions[regions@seqnames == gene_of_interest$seqname]
thalIDregion <- regions_scaffold[regions_scaffold@ranges@start == gene_of_interest$region_start]
}


### Plot gene of interest

# combine cov files

cov_files <- c(cov_1_diploid,
               cov_2_diploid,
               cov_3_diploid,
               cov_1_polyploid,
               cov_2_polyploid,
               cov_3_polyploid)

# set up dmrseq input processing

print("Loading cov files...")

bismarkBSseq <- read.bismark(files = c(cov_files),
                             rmZeroCov = TRUE, 
                             strandCollapse = FALSE,
                             verbose = TRUE)

print("cov file loaded!")

sampleNames <- c(rep("par", 3), 
                 rep("kam", 3))

pData(bismarkBSseq)$Species <- sampleNames

# Filtering step

loci.idx <- which(DelayedMatrixStats::rowSums2(getCoverage(bismarkBSseq, type="Cov")==0) == 0)
sample.idx <- which(pData(bismarkBSseq)$Species %in% c("par", "kam"))

bs.filtered <- bismarkBSseq[loci.idx, sample.idx]

print("Filtering completed!")

# Plot DMR

if(nrow(gene_of_interest) > 1){
  annoTrack <- read_annotations(bed_anno)
  for(i in 1:nrow(gene_of_interest)){
  png(filename = paste0(context, "_", ATID, "_", i, ".png"), width = 600, height = 600)
  plotDMRs(bs.filtered, regions=thalIDregion[i], testCovariate="Species")
  dev.off()
  }
} else {
annoTrack <- read_annotations(bed_anno)
png(filename = paste0(context, "_", ATID, ".png"), width = 600, height = 600)
plotDMRs(bs.filtered, regions=thalIDregion, testCovariate="Species") #annoTrack = annoTrack)
dev.off()
}

print("Script completed :)")