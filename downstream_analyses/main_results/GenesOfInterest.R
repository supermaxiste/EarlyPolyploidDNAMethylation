# Extracting genes of interest showing changes in
# both expression and methylation 

# Packages

library(data.table)
library(ggplot2)
library(patchwork)
library(tidyverse)

# Methylation files

setwd("~/OneDrive/PhD/Project/Chapter_3/DMR_results/")

# Cold conditions

# halleri-side

HM_hal_v_synG1_DMG_CG <- read.delim("ARPEGGIO_results_syn1Vpro1_HM/dmrseq/CG_context/parent1_v_allo_sig_diff.txt")

HM_hal_v_synG1_DMG_CHG <- read.delim("ARPEGGIO_results_syn1Vpro1_HM/dmrseq/CHG_context/parent1_v_allo_sig_diff.txt")

HM_hal_v_synG1_DMG_CHH <- read.delim("ARPEGGIO_results_syn1Vpro1_HM/dmrseq/CHH_context/parent1_v_allo_sig_diff.txt")

HM_hal_v_synG4_DMG_CG <- read.delim("ARPEGGIO_results_syn4Vpro1_HM/dmrseq/CG_context/parent1_v_allo_sig_diff.txt")

HM_hal_v_synG4_DMG_CHG <- read.delim("ARPEGGIO_results_syn4Vpro1_HM/dmrseq/CHG_context/parent1_v_allo_sig_diff.txt")

HM_hal_v_synG4_DMG_CHH <- read.delim("ARPEGGIO_results_syn4Vpro1_HM/dmrseq/CHH_context/parent1_v_allo_sig_diff.txt")

# lyrata-side

HM_lyr_v_synG1_DMG_CG <- read.delim("ARPEGGIO_results_syn1Vpro1_HM/dmrseq/CG_context/parent2_v_allo_sig_diff.txt")

HM_lyr_v_synG1_DMG_CHG <- read.delim("ARPEGGIO_results_syn1Vpro1_HM/dmrseq/CHG_context/parent2_v_allo_sig_diff.txt")

HM_lyr_v_synG1_DMG_CHH <- read.delim("ARPEGGIO_results_syn1Vpro1_HM/dmrseq/CHH_context/parent2_v_allo_sig_diff.txt")

HM_lyr_v_synG4_DMG_CG <- read.delim("ARPEGGIO_results_syn4Vpro1_HM/dmrseq/CG_context/parent2_v_allo_sig_diff.txt")

HM_lyr_v_synG4_DMG_CHG <- read.delim("ARPEGGIO_results_syn4Vpro1_HM/dmrseq/CHG_context/parent2_v_allo_sig_diff.txt")

HM_lyr_v_synG4_DMG_CHH <- read.delim("ARPEGGIO_results_syn4Vpro1_HM/dmrseq/CHH_context/parent2_v_allo_sig_diff.txt")

# Hot conditions

# halleri-side

LL_hal_v_synG1_DMG_CG <- read.delim("ARPEGGIO_results_syn1Vpro1_LL/dmrseq/CG_context/parent1_v_allo_sig_diff.txt")

LL_hal_v_synG1_DMG_CHG <- read.delim("ARPEGGIO_results_syn1Vpro1_LL/dmrseq/CHG_context/parent1_v_allo_sig_diff.txt")

LL_hal_v_synG1_DMG_CHH <- read.delim("ARPEGGIO_results_syn1Vpro1_LL/dmrseq/CHH_context/parent1_v_allo_sig_diff.txt")

LL_hal_v_synG4_DMG_CG <- read.delim("ARPEGGIO_results_syn4Vpro1_LL/dmrseq/CG_context/parent1_v_allo_sig_diff.txt")

LL_hal_v_synG4_DMG_CHG <- read.delim("ARPEGGIO_results_syn4Vpro1_LL/dmrseq/CHG_context/parent1_v_allo_sig_diff.txt")

LL_hal_v_synG4_DMG_CHH <- read.delim("ARPEGGIO_results_syn4Vpro1_LL/dmrseq/CHH_context/parent1_v_allo_sig_diff.txt")

# lyrata-side

LL_lyr_v_synG1_DMG_CG <- read.delim("ARPEGGIO_results_syn1Vpro1_LL/dmrseq/CG_context/parent2_v_allo_sig_diff.txt")

LL_lyr_v_synG1_DMG_CHG <- read.delim("ARPEGGIO_results_syn1Vpro1_LL/dmrseq/CHG_context/parent2_v_allo_sig_diff.txt")

LL_lyr_v_synG1_DMG_CHH <- read.delim("ARPEGGIO_results_syn1Vpro1_LL/dmrseq/CHH_context/parent2_v_allo_sig_diff.txt")

LL_lyr_v_synG4_DMG_CG <- read.delim("ARPEGGIO_results_syn4Vpro1_LL/dmrseq/CG_context/parent2_v_allo_sig_diff.txt")

LL_lyr_v_synG4_DMG_CHG <- read.delim("ARPEGGIO_results_syn4Vpro1_LL/dmrseq/CHG_context/parent2_v_allo_sig_diff.txt")

LL_lyr_v_synG4_DMG_CHH <- read.delim("ARPEGGIO_results_syn4Vpro1_LL/dmrseq/CHH_context/parent2_v_allo_sig_diff.txt")

# Since methylation files do not show the geneID
# We will use the output DM files from ARPEGGIO
# to add the geneID

# Cold conditions

# halleri-side

HM_hal_v_synG1_DMG_CG_geneIDs <- read.delim("ARPEGGIO_results_syn1Vpro1_HM/dmrseq/CG_context/DM_genes_parent1_v_allo_CG_context.txt")

HM_hal_v_synG1_DMG_CHG_geneIDs <- read.delim("ARPEGGIO_results_syn1Vpro1_HM/dmrseq/CHG_context/DM_genes_parent1_v_allo_CHG_context.txt")

HM_hal_v_synG1_DMG_CHH_geneIDs <- read.delim("ARPEGGIO_results_syn1Vpro1_HM/dmrseq/CHH_context/DM_genes_parent1_v_allo_CHH_context.txt")

HM_hal_v_synG4_DMG_CG_geneIDs <- read.delim("ARPEGGIO_results_syn4Vpro1_HM/dmrseq/CG_context/DM_genes_parent1_v_allo_CG_context.txt")

HM_hal_v_synG4_DMG_CHG_geneIDs <- read.delim("ARPEGGIO_results_syn4Vpro1_HM/dmrseq/CHG_context/DM_genes_parent1_v_allo_CHG_context.txt")

HM_hal_v_synG4_DMG_CHH_geneIDs <- read.delim("ARPEGGIO_results_syn4Vpro1_HM/dmrseq/CHH_context/DM_genes_parent1_v_allo_CHH_context.txt")

# lyrata-side

HM_lyr_v_synG1_DMG_CG_geneIDs <- read.delim("ARPEGGIO_results_syn1Vpro1_HM/dmrseq/CG_context/DM_genes_parent2_v_allo_CG_context.txt")

HM_lyr_v_synG1_DMG_CHG_geneIDs <- read.delim("ARPEGGIO_results_syn1Vpro1_HM/dmrseq/CHG_context/DM_genes_parent2_v_allo_CHG_context.txt")

HM_lyr_v_synG1_DMG_CHH_geneIDs <- read.delim("ARPEGGIO_results_syn1Vpro1_HM/dmrseq/CHH_context/DM_genes_parent2_v_allo_CHH_context.txt")

HM_lyr_v_synG4_DMG_CG_geneIDs <- read.delim("ARPEGGIO_results_syn4Vpro1_HM/dmrseq/CG_context/DM_genes_parent2_v_allo_CG_context.txt")

HM_lyr_v_synG4_DMG_CHG_geneIDs <- read.delim("ARPEGGIO_results_syn4Vpro1_HM/dmrseq/CHG_context/DM_genes_parent2_v_allo_CHG_context.txt")

HM_lyr_v_synG4_DMG_CHH_geneIDs <- read.delim("ARPEGGIO_results_syn4Vpro1_HM/dmrseq/CHH_context/DM_genes_parent2_v_allo_CHH_context.txt")

# Hot conditions

# halleri-side

LL_hal_v_synG1_DMG_CG_geneIDs <- read.delim("ARPEGGIO_results_syn1Vpro1_LL/dmrseq/CG_context/DM_genes_parent1_v_allo_CG_context.txt")

LL_hal_v_synG1_DMG_CHG_geneIDs <- read.delim("ARPEGGIO_results_syn1Vpro1_LL/dmrseq/CHG_context/DM_genes_parent1_v_allo_CHG_context.txt")

LL_hal_v_synG1_DMG_CHH_geneIDs <- read.delim("ARPEGGIO_results_syn1Vpro1_LL/dmrseq/CHH_context/DM_genes_parent1_v_allo_CHH_context.txt")

LL_hal_v_synG4_DMG_CG_geneIDs <- read.delim("ARPEGGIO_results_syn4Vpro1_LL/dmrseq/CG_context/DM_genes_parent1_v_allo_CG_context.txt")

LL_hal_v_synG4_DMG_CHG_geneIDs <- read.delim("ARPEGGIO_results_syn4Vpro1_LL/dmrseq/CHG_context/DM_genes_parent1_v_allo_CHG_context.txt")

LL_hal_v_synG4_DMG_CHH_geneIDs <- read.delim("ARPEGGIO_results_syn4Vpro1_LL/dmrseq/CHH_context/DM_genes_parent1_v_allo_CHH_context.txt")

# lyrata-side

LL_lyr_v_synG1_DMG_CG_geneIDs <- read.delim("ARPEGGIO_results_syn1Vpro1_LL/dmrseq/CG_context/DM_genes_parent2_v_allo_CG_context.txt")

LL_lyr_v_synG1_DMG_CHG_geneIDs <- read.delim("ARPEGGIO_results_syn1Vpro1_LL/dmrseq/CHG_context/DM_genes_parent2_v_allo_CHG_context.txt")

LL_lyr_v_synG1_DMG_CHH_geneIDs <- read.delim("ARPEGGIO_results_syn1Vpro1_LL/dmrseq/CHH_context/DM_genes_parent2_v_allo_CHH_context.txt")

LL_lyr_v_synG4_DMG_CG_geneIDs <- read.delim("ARPEGGIO_results_syn4Vpro1_LL/dmrseq/CG_context/DM_genes_parent2_v_allo_CG_context.txt")

LL_lyr_v_synG4_DMG_CHG_geneIDs <- read.delim("ARPEGGIO_results_syn4Vpro1_LL/dmrseq/CHG_context/DM_genes_parent2_v_allo_CHG_context.txt")

LL_lyr_v_synG4_DMG_CHH_geneIDs <- read.delim("ARPEGGIO_results_syn4Vpro1_LL/dmrseq/CHH_context/DM_genes_parent2_v_allo_CHH_context.txt")

# We add the geneID to the dmrseq output files
# First we create a unique identifier for the DMR position

create_ID1 <- function(file){
  new_ID <- paste0(file$seqname,
                   "_",
                   file$start,
                   "_",
                   file$end)
  return(new_ID)
}

create_ID2 <- function(file){
  new_ID <- paste0(file$seqname,
                   "_",
                   file$region_start,
                   "_",
                   file$region_end)
  return(new_ID)
}

# halleri-side 

HM_hal_v_synG1_DMG_CG_geneIDs_pos <- create_ID2(HM_hal_v_synG1_DMG_CG_geneIDs)
HM_hal_v_synG1_DMG_CG_pos <- create_ID1(HM_hal_v_synG1_DMG_CG)
HM_hal_v_synG1_DMG_CHG_geneIDs_pos <- create_ID2(HM_hal_v_synG1_DMG_CHG_geneIDs)
HM_hal_v_synG1_DMG_CHG_pos <- create_ID1(HM_hal_v_synG1_DMG_CHG)
HM_hal_v_synG1_DMG_CHH_geneIDs_pos <- create_ID2(HM_hal_v_synG1_DMG_CHH_geneIDs)
HM_hal_v_synG1_DMG_CHH_pos <- create_ID1(HM_hal_v_synG1_DMG_CHH)

HM_hal_v_synG4_DMG_CG_geneIDs_pos <- create_ID2(HM_hal_v_synG4_DMG_CG_geneIDs)
HM_hal_v_synG4_DMG_CG_pos <- create_ID1(HM_hal_v_synG4_DMG_CG)
HM_hal_v_synG4_DMG_CHG_geneIDs_pos <- create_ID2(HM_hal_v_synG4_DMG_CHG_geneIDs)
HM_hal_v_synG4_DMG_CHG_pos <- create_ID1(HM_hal_v_synG4_DMG_CHG)
HM_hal_v_synG4_DMG_CHH_geneIDs_pos <- create_ID2(HM_hal_v_synG4_DMG_CHH_geneIDs)
HM_hal_v_synG4_DMG_CHH_pos <- create_ID1(HM_hal_v_synG4_DMG_CHH)

LL_hal_v_synG1_DMG_CG_geneIDs_pos <- create_ID2(LL_hal_v_synG1_DMG_CG_geneIDs)
LL_hal_v_synG1_DMG_CG_pos <- create_ID1(LL_hal_v_synG1_DMG_CG)
LL_hal_v_synG1_DMG_CHG_geneIDs_pos <- create_ID2(LL_hal_v_synG1_DMG_CHG_geneIDs)
LL_hal_v_synG1_DMG_CHG_pos <- create_ID1(LL_hal_v_synG1_DMG_CHG)
LL_hal_v_synG1_DMG_CHH_geneIDs_pos <- create_ID2(LL_hal_v_synG1_DMG_CHH_geneIDs)
LL_hal_v_synG1_DMG_CHH_pos <- create_ID1(LL_hal_v_synG1_DMG_CHH)

LL_hal_v_synG4_DMG_CG_geneIDs_pos <- create_ID2(LL_hal_v_synG4_DMG_CG_geneIDs)
LL_hal_v_synG4_DMG_CG_pos <- create_ID1(LL_hal_v_synG4_DMG_CG)
LL_hal_v_synG4_DMG_CHG_geneIDs_pos <- create_ID2(LL_hal_v_synG4_DMG_CHG_geneIDs)
LL_hal_v_synG4_DMG_CHG_pos <- create_ID1(LL_hal_v_synG4_DMG_CHG)
LL_hal_v_synG4_DMG_CHH_geneIDs_pos <- create_ID2(LL_hal_v_synG4_DMG_CHH_geneIDs)
LL_hal_v_synG4_DMG_CHH_pos <- create_ID1(LL_hal_v_synG4_DMG_CHH)

# lyrata side

HM_lyr_v_synG1_DMG_CG_geneIDs_pos <- create_ID2(HM_lyr_v_synG1_DMG_CG_geneIDs)
HM_lyr_v_synG1_DMG_CG_pos <- create_ID1(HM_lyr_v_synG1_DMG_CG)
HM_lyr_v_synG1_DMG_CHG_geneIDs_pos <- create_ID2(HM_lyr_v_synG1_DMG_CHG_geneIDs)
HM_lyr_v_synG1_DMG_CHG_pos <- create_ID1(HM_lyr_v_synG1_DMG_CHG)
HM_lyr_v_synG1_DMG_CHH_geneIDs_pos <- create_ID2(HM_lyr_v_synG1_DMG_CHH_geneIDs)
HM_lyr_v_synG1_DMG_CHH_pos <- create_ID1(HM_lyr_v_synG1_DMG_CHH)

HM_lyr_v_synG4_DMG_CG_geneIDs_pos <- create_ID2(HM_lyr_v_synG4_DMG_CG_geneIDs)
HM_lyr_v_synG4_DMG_CG_pos <- create_ID1(HM_lyr_v_synG4_DMG_CG)
HM_lyr_v_synG4_DMG_CHG_geneIDs_pos <- create_ID2(HM_lyr_v_synG4_DMG_CHG_geneIDs)
HM_lyr_v_synG4_DMG_CHG_pos <- create_ID1(HM_lyr_v_synG4_DMG_CHG)
HM_lyr_v_synG4_DMG_CHH_geneIDs_pos <- create_ID2(HM_lyr_v_synG4_DMG_CHH_geneIDs)
HM_lyr_v_synG4_DMG_CHH_pos <- create_ID1(HM_lyr_v_synG4_DMG_CHH)

LL_lyr_v_synG1_DMG_CG_geneIDs_pos <- create_ID2(LL_lyr_v_synG1_DMG_CG_geneIDs)
LL_lyr_v_synG1_DMG_CG_pos <- create_ID1(LL_lyr_v_synG1_DMG_CG)
LL_lyr_v_synG1_DMG_CHG_geneIDs_pos <- create_ID2(LL_lyr_v_synG1_DMG_CHG_geneIDs)
LL_lyr_v_synG1_DMG_CHG_pos <- create_ID1(LL_lyr_v_synG1_DMG_CHG)
LL_lyr_v_synG1_DMG_CHH_geneIDs_pos <- create_ID2(LL_lyr_v_synG1_DMG_CHH_geneIDs)
LL_lyr_v_synG1_DMG_CHH_pos <- create_ID1(LL_lyr_v_synG1_DMG_CHH)

LL_lyr_v_synG4_DMG_CG_geneIDs_pos <- create_ID2(LL_lyr_v_synG4_DMG_CG_geneIDs)
LL_lyr_v_synG4_DMG_CG_pos <- create_ID1(LL_lyr_v_synG4_DMG_CG)
LL_lyr_v_synG4_DMG_CHG_geneIDs_pos <- create_ID2(LL_lyr_v_synG4_DMG_CHG_geneIDs)
LL_lyr_v_synG4_DMG_CHG_pos <- create_ID1(LL_lyr_v_synG4_DMG_CHG)
LL_lyr_v_synG4_DMG_CHH_geneIDs_pos <- create_ID2(LL_lyr_v_synG4_DMG_CHH_geneIDs)
LL_lyr_v_synG4_DMG_CHH_pos <- create_ID1(LL_lyr_v_synG4_DMG_CHH)

# function to add raw methylation difference
# to DM files

add_rawMeth <- function(DM_pos_file, 
                        dmrseq_pos_file,
                        DM_file,
                        dmrseq_file){
  index <- match(DM_pos_file, dmrseq_pos_file)
  new_DM_file <- mutate(DM_file,
                        # raw methylation is based on synthetics as reference, sign needs to be inverted to have progenitors as reference
                        rawMeth = -dmrseq_file$rawDiff[index])
  return(new_DM_file)
}

#halleri side

HM_hal_v_synG1_DMG_CG_final <- add_rawMeth(HM_hal_v_synG1_DMG_CG_geneIDs_pos,
                                           HM_hal_v_synG1_DMG_CG_pos,
                                           HM_hal_v_synG1_DMG_CG_geneIDs,
                                           HM_hal_v_synG1_DMG_CG)

HM_hal_v_synG1_DMG_CHG_final <- add_rawMeth(HM_hal_v_synG1_DMG_CHG_geneIDs_pos,
                                            HM_hal_v_synG1_DMG_CHG_pos,
                                            HM_hal_v_synG1_DMG_CHG_geneIDs,
                                            HM_hal_v_synG1_DMG_CHG)

HM_hal_v_synG1_DMG_CHH_final <- add_rawMeth(HM_hal_v_synG1_DMG_CHH_geneIDs_pos,
                                            HM_hal_v_synG1_DMG_CHH_pos,
                                            HM_hal_v_synG1_DMG_CHH_geneIDs,
                                            HM_hal_v_synG1_DMG_CHH)

HM_hal_v_synG4_DMG_CG_final <- add_rawMeth(HM_hal_v_synG4_DMG_CG_geneIDs_pos,
                                           HM_hal_v_synG4_DMG_CG_pos,
                                           HM_hal_v_synG4_DMG_CG_geneIDs,
                                           HM_hal_v_synG4_DMG_CG)

HM_hal_v_synG4_DMG_CHG_final <- add_rawMeth(HM_hal_v_synG4_DMG_CHG_geneIDs_pos,
                                            HM_hal_v_synG4_DMG_CHG_pos,
                                            HM_hal_v_synG4_DMG_CHG_geneIDs,
                                            HM_hal_v_synG4_DMG_CHG)

HM_hal_v_synG4_DMG_CHH_final <- add_rawMeth(HM_hal_v_synG4_DMG_CHH_geneIDs_pos,
                                            HM_hal_v_synG4_DMG_CHH_pos,
                                            HM_hal_v_synG4_DMG_CHH_geneIDs,
                                            HM_hal_v_synG4_DMG_CHH)

LL_hal_v_synG1_DMG_CG_final <- add_rawMeth(LL_hal_v_synG1_DMG_CG_geneIDs_pos,
                                           LL_hal_v_synG1_DMG_CG_pos,
                                           LL_hal_v_synG1_DMG_CG_geneIDs,
                                           LL_hal_v_synG1_DMG_CG)

LL_hal_v_synG1_DMG_CHG_final <- add_rawMeth(LL_hal_v_synG1_DMG_CHG_geneIDs_pos,
                                            LL_hal_v_synG1_DMG_CHG_pos,
                                            LL_hal_v_synG1_DMG_CHG_geneIDs,
                                            LL_hal_v_synG1_DMG_CHG)

LL_hal_v_synG1_DMG_CHH_final <- add_rawMeth(LL_hal_v_synG1_DMG_CHH_geneIDs_pos,
                                            LL_hal_v_synG1_DMG_CHH_pos,
                                            LL_hal_v_synG1_DMG_CHH_geneIDs,
                                            LL_hal_v_synG1_DMG_CHH)

LL_hal_v_synG4_DMG_CG_final <- add_rawMeth(LL_hal_v_synG4_DMG_CG_geneIDs_pos,
                                           LL_hal_v_synG4_DMG_CG_pos,
                                           LL_hal_v_synG4_DMG_CG_geneIDs,
                                           LL_hal_v_synG4_DMG_CG)

LL_hal_v_synG4_DMG_CHG_final <- add_rawMeth(LL_hal_v_synG4_DMG_CHG_geneIDs_pos,
                                            LL_hal_v_synG4_DMG_CHG_pos,
                                            LL_hal_v_synG4_DMG_CHG_geneIDs,
                                            LL_hal_v_synG4_DMG_CHG)

LL_hal_v_synG4_DMG_CHH_final <- add_rawMeth(LL_hal_v_synG4_DMG_CHH_geneIDs_pos,
                                            LL_hal_v_synG4_DMG_CHH_pos,
                                            LL_hal_v_synG4_DMG_CHH_geneIDs,
                                            LL_hal_v_synG4_DMG_CHH)

#lyrata side

HM_lyr_v_synG1_DMG_CG_final <- add_rawMeth(HM_lyr_v_synG1_DMG_CG_geneIDs_pos,
                                           HM_lyr_v_synG1_DMG_CG_pos,
                                           HM_lyr_v_synG1_DMG_CG_geneIDs,
                                           HM_lyr_v_synG1_DMG_CG)

HM_lyr_v_synG1_DMG_CHG_final <- add_rawMeth(HM_lyr_v_synG1_DMG_CHG_geneIDs_pos,
                                            HM_lyr_v_synG1_DMG_CHG_pos,
                                            HM_lyr_v_synG1_DMG_CHG_geneIDs,
                                            HM_lyr_v_synG1_DMG_CHG)

HM_lyr_v_synG1_DMG_CHH_final <- add_rawMeth(HM_lyr_v_synG1_DMG_CHH_geneIDs_pos,
                                            HM_lyr_v_synG1_DMG_CHH_pos,
                                            HM_lyr_v_synG1_DMG_CHH_geneIDs,
                                            HM_lyr_v_synG1_DMG_CHH)

HM_lyr_v_synG4_DMG_CG_final <- add_rawMeth(HM_lyr_v_synG4_DMG_CG_geneIDs_pos,
                                           HM_lyr_v_synG4_DMG_CG_pos,
                                           HM_lyr_v_synG4_DMG_CG_geneIDs,
                                           HM_lyr_v_synG4_DMG_CG)

HM_lyr_v_synG4_DMG_CHG_final <- add_rawMeth(HM_lyr_v_synG4_DMG_CHG_geneIDs_pos,
                                            HM_lyr_v_synG4_DMG_CHG_pos,
                                            HM_lyr_v_synG4_DMG_CHG_geneIDs,
                                            HM_lyr_v_synG4_DMG_CHG)

HM_lyr_v_synG4_DMG_CHH_final <- add_rawMeth(HM_lyr_v_synG4_DMG_CHH_geneIDs_pos,
                                            HM_lyr_v_synG4_DMG_CHH_pos,
                                            HM_lyr_v_synG4_DMG_CHH_geneIDs,
                                            HM_lyr_v_synG4_DMG_CHH)

LL_lyr_v_synG1_DMG_CG_final <- add_rawMeth(LL_lyr_v_synG1_DMG_CG_geneIDs_pos,
                                           LL_lyr_v_synG1_DMG_CG_pos,
                                           LL_lyr_v_synG1_DMG_CG_geneIDs,
                                           LL_lyr_v_synG1_DMG_CG)

LL_lyr_v_synG1_DMG_CHG_final <- add_rawMeth(LL_lyr_v_synG1_DMG_CHG_geneIDs_pos,
                                            LL_lyr_v_synG1_DMG_CHG_pos,
                                            LL_lyr_v_synG1_DMG_CHG_geneIDs,
                                            LL_lyr_v_synG1_DMG_CHG)

LL_lyr_v_synG1_DMG_CHH_final <- add_rawMeth(LL_lyr_v_synG1_DMG_CHH_geneIDs_pos,
                                            LL_lyr_v_synG1_DMG_CHH_pos,
                                            LL_lyr_v_synG1_DMG_CHH_geneIDs,
                                            LL_lyr_v_synG1_DMG_CHH)

LL_lyr_v_synG4_DMG_CG_final <- add_rawMeth(LL_lyr_v_synG4_DMG_CG_geneIDs_pos,
                                           LL_lyr_v_synG4_DMG_CG_pos,
                                           LL_lyr_v_synG4_DMG_CG_geneIDs,
                                           LL_lyr_v_synG4_DMG_CG)

LL_lyr_v_synG4_DMG_CHG_final <- add_rawMeth(LL_lyr_v_synG4_DMG_CHG_geneIDs_pos,
                                            LL_lyr_v_synG4_DMG_CHG_pos,
                                            LL_lyr_v_synG4_DMG_CHG_geneIDs,
                                            LL_lyr_v_synG4_DMG_CHG)

LL_lyr_v_synG4_DMG_CHH_final <- add_rawMeth(LL_lyr_v_synG4_DMG_CHH_geneIDs_pos,
                                            LL_lyr_v_synG4_DMG_CHH_pos,
                                            LL_lyr_v_synG4_DMG_CHH_geneIDs,
                                            LL_lyr_v_synG4_DMG_CHH)

# Expression files 

setwd("~/OneDrive/PhD/Project/Chapter_3/Reports/RNAseqDEG")

# Cold conditions

HM_hal_v_synG1_DEG <- read.delim("HM_halG1_v_RS7G1_DEG.txt")
HM_hal_v_synG4_DEG <- read.delim("HM_halG1_v_RS7G4_DEG.txt")
HM_lyr_v_synG1_DEG <- read.delim("HM_lyrG1_v_RS7G1_DEG.txt")
HM_lyr_v_synG4_DEG <- read.delim("HM_lyrG1_v_RS7G4_DEG.txt")

# Hot conditions

LL_hal_v_synG1_DEG <- read.delim("LL_halG1_v_RS7G1_DEG.txt")
LL_hal_v_synG4_DEG <- read.delim("LL_halG1_v_RS7G4_DEG.txt")
LL_lyr_v_synG1_DEG <- read.delim("LL_lyrG1_v_RS7G1_DEG.txt")
LL_lyr_v_synG4_DEG <- read.delim("LL_lyrG1_v_RS7G4_DEG.txt")


# Find genes showing both methylation and expression changes
# and add logFC to DM file

add_logFC <- function(DM_final, DEG_file){
  index <- match(DM_final$geneID, DEG_file$geneID)
  new_DM_file <- mutate(DM_final,
                        logFC = DEG_file$logFC[index])
  return(new_DM_file)
}

# cold conditions 

HM_hal_v_synG1_DMG_CG_all <- add_logFC(HM_hal_v_synG1_DMG_CG_final,
                                       HM_hal_v_synG1_DEG)
HM_hal_v_synG1_DMG_CHG_all <- add_logFC(HM_hal_v_synG1_DMG_CHG_final,
                                        HM_hal_v_synG1_DEG)
HM_hal_v_synG1_DMG_CHH_all <- add_logFC(HM_hal_v_synG1_DMG_CHH_final,
                                        HM_hal_v_synG1_DEG)

HM_hal_v_synG4_DMG_CG_all <- add_logFC(HM_hal_v_synG4_DMG_CG_final,
                                       HM_hal_v_synG4_DEG)
HM_hal_v_synG4_DMG_CHG_all <- add_logFC(HM_hal_v_synG4_DMG_CHG_final,
                                        HM_hal_v_synG4_DEG)
HM_hal_v_synG4_DMG_CHH_all <- add_logFC(HM_hal_v_synG4_DMG_CHH_final,
                                        HM_hal_v_synG4_DEG)

HM_lyr_v_synG1_DMG_CG_all <- add_logFC(HM_lyr_v_synG1_DMG_CG_final,
                                       HM_lyr_v_synG1_DEG)
HM_lyr_v_synG1_DMG_CHG_all <- add_logFC(HM_lyr_v_synG1_DMG_CHG_final,
                                        HM_lyr_v_synG1_DEG)
HM_lyr_v_synG1_DMG_CHH_all <- add_logFC(HM_lyr_v_synG1_DMG_CHH_final,
                                        HM_lyr_v_synG1_DEG)

HM_lyr_v_synG4_DMG_CG_all <- add_logFC(HM_lyr_v_synG4_DMG_CG_final,
                                       HM_lyr_v_synG4_DEG)
HM_lyr_v_synG4_DMG_CHG_all <- add_logFC(HM_lyr_v_synG4_DMG_CHG_final,
                                        HM_lyr_v_synG4_DEG)
HM_lyr_v_synG4_DMG_CHH_all <- add_logFC(HM_lyr_v_synG4_DMG_CHH_final,
                                        HM_lyr_v_synG4_DEG)

# hot conditions 

LL_hal_v_synG1_DMG_CG_all <- add_logFC(LL_hal_v_synG1_DMG_CG_final,
                                       LL_hal_v_synG1_DEG)
LL_hal_v_synG1_DMG_CHG_all <- add_logFC(LL_hal_v_synG1_DMG_CHG_final,
                                        LL_hal_v_synG1_DEG)
LL_hal_v_synG1_DMG_CHH_all <- add_logFC(LL_hal_v_synG1_DMG_CHH_final,
                                        LL_hal_v_synG1_DEG)

LL_hal_v_synG4_DMG_CG_all <- add_logFC(LL_hal_v_synG4_DMG_CG_final,
                                       LL_hal_v_synG4_DEG)
LL_hal_v_synG4_DMG_CHG_all <- add_logFC(LL_hal_v_synG4_DMG_CHG_final,
                                        LL_hal_v_synG4_DEG)
LL_hal_v_synG4_DMG_CHH_all <- add_logFC(LL_hal_v_synG4_DMG_CHH_final,
                                        LL_hal_v_synG4_DEG)

LL_lyr_v_synG1_DMG_CG_all <- add_logFC(LL_lyr_v_synG1_DMG_CG_final,
                                       LL_lyr_v_synG1_DEG)
LL_lyr_v_synG1_DMG_CHG_all <- add_logFC(LL_lyr_v_synG1_DMG_CHG_final,
                                        LL_lyr_v_synG1_DEG)
LL_lyr_v_synG1_DMG_CHH_all <- add_logFC(LL_lyr_v_synG1_DMG_CHH_final,
                                        LL_lyr_v_synG1_DEG)

LL_lyr_v_synG4_DMG_CG_all <- add_logFC(LL_lyr_v_synG4_DMG_CG_final,
                                       LL_lyr_v_synG4_DEG)
LL_lyr_v_synG4_DMG_CHG_all <- add_logFC(LL_lyr_v_synG4_DMG_CHG_final,
                                        LL_lyr_v_synG4_DEG)
LL_lyr_v_synG4_DMG_CHH_all <- add_logFC(LL_lyr_v_synG4_DMG_CHH_final,
                                        LL_lyr_v_synG4_DEG)

# Merge DMGs and remove the ones with no expression change

# halleri

HM_hal_v_synG1_DMG_all <- na.omit(rbind(
  HM_hal_v_synG1_DMG_CG_all,
  HM_hal_v_synG1_DMG_CHG_all,
  HM_hal_v_synG1_DMG_CHH_all))

HM_hal_v_synG4_DMG_all <- na.omit(rbind(
  HM_hal_v_synG4_DMG_CG_all,
  HM_hal_v_synG4_DMG_CHG_all,
  HM_hal_v_synG4_DMG_CHH_all))

LL_hal_v_synG1_DMG_all <- na.omit(rbind(
  LL_hal_v_synG1_DMG_CG_all,
  LL_hal_v_synG1_DMG_CHG_all,
  LL_hal_v_synG1_DMG_CHH_all))

LL_hal_v_synG4_DMG_all <- na.omit(rbind(
  LL_hal_v_synG4_DMG_CG_all,
  LL_hal_v_synG4_DMG_CHG_all,
  LL_hal_v_synG4_DMG_CHH_all))

# lyrata

HM_lyr_v_synG1_DMG_all <- na.omit(rbind(
  HM_lyr_v_synG1_DMG_CG_all,
  HM_lyr_v_synG1_DMG_CHG_all,
  HM_lyr_v_synG1_DMG_CHH_all))

HM_lyr_v_synG4_DMG_all <- na.omit(rbind(
  HM_lyr_v_synG4_DMG_CG_all,
  HM_lyr_v_synG4_DMG_CHG_all,
  HM_lyr_v_synG4_DMG_CHH_all))

LL_lyr_v_synG1_DMG_all <- na.omit(rbind(
  LL_lyr_v_synG1_DMG_CG_all,
  LL_lyr_v_synG1_DMG_CHG_all,
  LL_lyr_v_synG1_DMG_CHH_all))

LL_lyr_v_synG4_DMG_all <- na.omit(rbind(
  LL_lyr_v_synG4_DMG_CG_all,
  LL_lyr_v_synG4_DMG_CHG_all,
  LL_lyr_v_synG4_DMG_CHH_all))

# Find overlapping genes

overlap_hal_1 <- base::intersect(unique(HM_hal_v_synG1_DMG_all$geneID),
                            unique(HM_hal_v_synG4_DMG_all$geneID))

overlap_hal_2 <- base::intersect(unique(LL_hal_v_synG1_DMG_all$geneID),
                            unique(LL_hal_v_synG4_DMG_all$geneID))

overlap_hal <- intersect(overlap_hal_1,
                         overlap_hal_2)

overlap_lyr_1 <- base::intersect(unique(HM_lyr_v_synG1_DMG_all$geneID),
                                 unique(HM_lyr_v_synG4_DMG_all$geneID))

overlap_lyr_2 <- base::intersect(unique(LL_lyr_v_synG1_DMG_all$geneID),
                                 unique(LL_lyr_v_synG4_DMG_all$geneID))

overlap_lyr <- intersect(overlap_lyr_1,
                         overlap_lyr_2)

# Import reciprocal best hit files
# to thaliana genes

rbh_Ahal <- read.delim("~/OneDrive/PhD/Project/Chapter2/rbh_hit_list_Ahal_v2_2_to_TAIR10_cds")

rbh_Ahal$Ahal <- unlist(strsplit(as.character(rbh_Ahal$Ahal), split = "[.]"))[rep(c(F,T), 21152)]

rbh_Alyr <- read.delim("~/OneDrive/PhD/Project/Chapter2/rbh_hit_list_Alyr_v2_2_to_TAIR10_cds")

rbh_Alyr$Alyr <- unlist(strsplit(as.character(rbh_Alyr$Alyr), split = "[.]"))[rep(c(F,T), 20962)]

# Match overlaps to A. thaliana gene IDs

hal_overlap_IDs <- na.omit(rbh_Ahal[match(overlap_hal, rbh_Ahal$Ahal),2])

lyr_overlap_IDs <- na.omit(rbh_Alyr[match(overlap_lyr, rbh_Alyr$Alyr),2])

all_overlap_IDs <- intersect(hal_overlap_IDs,
                             lyr_overlap_IDs)

## Create data frames with all overlaps and their expression
## or methylation change

hal_overlap_IDs_all <- na.omit(rbh_Ahal[match(overlap_hal, rbh_Ahal$Ahal),])

lyr_overlap_IDs_all <- na.omit(rbh_Alyr[match(overlap_lyr, rbh_Alyr$Alyr),])

get_rawMeth_hal <- function(file_all, overlap){
  index <- match(overlap$Ahal,
                 file_all$geneID)
  rawMeth <- file_all$rawMeth[index]
  return(rawMeth)
}

get_expr_hal <- function(DEG_file, overlap){
  index <- match(overlap$Ahal,
                 DEG_file$geneID)
  expr <- DEG_file$logFC[index]
  return(expr)
}

get_rawMeth_lyr <- function(file_all, overlap){
  index <- match(overlap$Alyr,
                 file_all$geneID)
  rawMeth <- file_all$rawMeth[index]
  return(rawMeth)
}

get_expr_lyr <- function(DEG_file, overlap){
  index <- match(overlap$Alyr,
                 DEG_file$geneID)
  expr <- DEG_file$logFC[index]
  return(expr)
}


hal_overlaps_DEMG <- data.frame(
  geneID = hal_overlap_IDs,
  HMrawMethCG_G1 = get_rawMeth_hal(HM_hal_v_synG1_DMG_CG_all,
                                   hal_overlap_IDs_all),
  HMrawMethCHG_G1 = get_rawMeth_hal(HM_hal_v_synG1_DMG_CHG_all,
                                    hal_overlap_IDs_all),
  HMrawMethCHH_G1 = get_rawMeth_hal(HM_hal_v_synG1_DMG_CHH_all,
                                    hal_overlap_IDs_all),
  HMlogFC_G1 = get_expr_hal(HM_hal_v_synG1_DEG,
                            hal_overlap_IDs_all),
  
  HMrawMethCG_G4 = get_rawMeth_hal(HM_hal_v_synG4_DMG_CG_all,
                                   hal_overlap_IDs_all),
  HMrawMethCHG_G4 = get_rawMeth_hal(HM_hal_v_synG4_DMG_CHG_all,
                                    hal_overlap_IDs_all),
  HMrawMethCHH_G4 = get_rawMeth_hal(HM_hal_v_synG4_DMG_CHH_all,
                                    hal_overlap_IDs_all),
  HMlogFC_G4 = get_expr_hal(HM_hal_v_synG4_DEG,
                            hal_overlap_IDs_all),
  LLrawMethCG_G1 = get_rawMeth_hal(LL_hal_v_synG1_DMG_CG_all,
                                   hal_overlap_IDs_all),
  LLrawMethCHG_G1 = get_rawMeth_hal(LL_hal_v_synG1_DMG_CHG_all,
                                    hal_overlap_IDs_all),
  LLrawMethCHH_G1 = get_rawMeth_hal(LL_hal_v_synG1_DMG_CHH_all,
                                    hal_overlap_IDs_all),
  LLlogFC_G1 = get_expr_hal(LL_hal_v_synG1_DEG,
                          hal_overlap_IDs_all),
  LLrawMethCG_G4 = get_rawMeth_hal(LL_hal_v_synG4_DMG_CG_all,
                                   hal_overlap_IDs_all),
  LLrawMethCHG_G4 = get_rawMeth_hal(LL_hal_v_synG4_DMG_CHG_all,
                                    hal_overlap_IDs_all),
  LLrawMethCHH_G4 = get_rawMeth_hal(LL_hal_v_synG4_DMG_CHH_all,
                                    hal_overlap_IDs_all),
  LLlogFC_G4 = get_expr_hal(LL_hal_v_synG4_DEG,
                          hal_overlap_IDs_all)
)

lyr_overlaps_DEMG <- data.frame(
  geneID = lyr_overlap_IDs,
  HMrawMethCG_G1 = get_rawMeth_lyr(HM_lyr_v_synG1_DMG_CG_all,
                                   lyr_overlap_IDs_all),
  HMrawMethCHG_G1 = get_rawMeth_lyr(HM_lyr_v_synG1_DMG_CHG_all,
                                    lyr_overlap_IDs_all),
  HMrawMethCHH_G1 = get_rawMeth_lyr(HM_lyr_v_synG1_DMG_CHH_all,
                                    lyr_overlap_IDs_all),
  HMlogFC_G1 = get_expr_lyr(HM_lyr_v_synG1_DEG,
                            lyr_overlap_IDs_all),
  
  HMrawMethCG_G4 = get_rawMeth_lyr(HM_lyr_v_synG4_DMG_CG_all,
                                   lyr_overlap_IDs_all),
  HMrawMethCHG_G4 = get_rawMeth_lyr(HM_lyr_v_synG4_DMG_CHG_all,
                                    lyr_overlap_IDs_all),
  HMrawMethCHH_G4 = get_rawMeth_lyr(HM_lyr_v_synG4_DMG_CHH_all,
                                    lyr_overlap_IDs_all),
  HMlogFC_G4 = get_expr_lyr(HM_lyr_v_synG4_DEG,
                            lyr_overlap_IDs_all),
  LLrawMethCG_G1 = get_rawMeth_lyr(LL_lyr_v_synG1_DMG_CG_all,
                                   lyr_overlap_IDs_all),
  LLrawMethCHG_G1 = get_rawMeth_lyr(LL_lyr_v_synG1_DMG_CHG_all,
                                    lyr_overlap_IDs_all),
  LLrawMethCHH_G1 = get_rawMeth_lyr(LL_lyr_v_synG1_DMG_CHH_all,
                                    lyr_overlap_IDs_all),
  LLlogFC_G1 = get_expr_lyr(LL_lyr_v_synG1_DEG,
                            lyr_overlap_IDs_all),
  LLrawMethCG_G4 = get_rawMeth_lyr(LL_lyr_v_synG4_DMG_CG_all,
                                   lyr_overlap_IDs_all),
  LLrawMethCHG_G4 = get_rawMeth_lyr(LL_lyr_v_synG4_DMG_CHG_all,
                                    lyr_overlap_IDs_all),
  LLrawMethCHH_G4 = get_rawMeth_lyr(LL_lyr_v_synG4_DMG_CHH_all,
                                    lyr_overlap_IDs_all),
  LLlogFC_G4 = get_expr_lyr(LL_lyr_v_synG4_DEG,
                            lyr_overlap_IDs_all)
)

# Export

# write.table(hal_overlap_IDs,
#             file = "Ahal_genesOfInterest.txt",
#             quote = F,
#             row.names = F,
#             col.names = F)
# 
# write.table(lyr_overlap_IDs,
#             file = "Alyr_genesOfInterest.txt",
#             quote = F,
#             row.names = F,
#             col.names = F)
