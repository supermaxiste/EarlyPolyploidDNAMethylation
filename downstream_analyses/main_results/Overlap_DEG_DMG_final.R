### This script will compute and plot the overlap
### between differentially methylated and expressed
### genes. Inputs are from dmrseq outputs for
### differential methylation analysis and edgeR for
### differential expression analysis


### Import libraries

library(tidyverse)
library(patchwork)
library(VennDiagram)
library(ggVennDiagram)
library(RColorBrewer)

### Import methylation data

setwd("/Users/ste/OneDrive/PhD/Project/Chapter_3/DMR_results/")

# Cold conditions (HM)

### halleri side

HM_halG1_v_RS7G1_DMG_CG <- read.delim("ARPEGGIO_results_syn1Vpro1_HM/dmrseq/CG_context/DM_genes_parent1_v_allo_CG_context.txt")
HM_halG1_v_RS7G1_DMG_CHG <- read.delim("ARPEGGIO_results_syn1Vpro1_HM/dmrseq/CHG_context/DM_genes_parent1_v_allo_CHG_context.txt")
HM_halG1_v_RS7G1_DMG_CHH <- read.delim("ARPEGGIO_results_syn1Vpro1_HM/dmrseq/CHH_context/DM_genes_parent1_v_allo_CHH_context.txt")

HM_halG1_v_RS7G4_DMG_CG <- read.delim("ARPEGGIO_results_syn4Vpro1_HM/dmrseq/CG_context/DM_genes_parent1_v_allo_CG_context.txt")
HM_halG1_v_RS7G4_DMG_CHG <- read.delim("ARPEGGIO_results_syn4Vpro1_HM/dmrseq/CHG_context/DM_genes_parent1_v_allo_CHG_context.txt")
HM_halG1_v_RS7G4_DMG_CHH <- read.delim("ARPEGGIO_results_syn4Vpro1_HM/dmrseq/CHH_context/DM_genes_parent1_v_allo_CHH_context.txt")

HM_RS7G1_v_RS7G4_DMG_CG_hal <- read.delim("ARPEGGIO_results_syn4v1_HM/dmrseq/CG_context/DM_genes_A_v_B_polyploid_CG_context_1.txt")
HM_RS7G1_v_RS7G4_DMG_CHG_hal <- read.delim("ARPEGGIO_results_syn4v1_HM/dmrseq/CHG_context/DM_genes_A_v_B_polyploid_CHG_context_1.txt")
HM_RS7G1_v_RS7G4_DMG_CHH_hal <- read.delim("ARPEGGIO_results_syn4v1_HM/dmrseq/CHH_context/DM_genes_A_v_B_polyploid_CHH_context_1.txt")

HM_ALKG1_v_halG1_DMG_CG <- read.delim("ARPEGGIO_results_alk1Vpro1_HM/dmrseq/CG_context/DM_genes_parent1_v_allo_CG_context.txt")
HM_ALKG1_v_halG1_DMG_CHG <- read.delim("ARPEGGIO_results_alk1Vpro1_HM/dmrseq/CHG_context/DM_genes_parent1_v_allo_CHG_context.txt")
HM_ALKG1_v_halG1_DMG_CHH <- read.delim("ARPEGGIO_results_alk1Vpro1_HM/dmrseq/CHH_context/DM_genes_parent1_v_allo_CHH_context.txt")

HM_ALKG1_v_RS7G1_DMG_CG_hal <- read.delim("ARPEGGIO_results_alk1Vsyn1_HM/dmrseq/CG_context/DM_genes_A_v_B_polyploid_CG_context_1.txt")
HM_ALKG1_v_RS7G1_DMG_CHG_hal <- read.delim("ARPEGGIO_results_alk1Vsyn1_HM/dmrseq/CHG_context/DM_genes_A_v_B_polyploid_CHG_context_1.txt")
HM_ALKG1_v_RS7G1_DMG_CHH_hal <- read.delim("ARPEGGIO_results_alk1Vsyn1_HM/dmrseq/CHH_context/DM_genes_A_v_B_polyploid_CHH_context_1.txt")

HM_ALKG1_v_RS7G4_DMG_CG_hal <- read.delim("ARPEGGIO_results_alk1Vsyn4_HM/dmrseq/CG_context/DM_genes_A_v_B_polyploid_CG_context_1.txt")
HM_ALKG1_v_RS7G4_DMG_CHG_hal <- read.delim("ARPEGGIO_results_alk1Vsyn4_HM/dmrseq/CHG_context/DM_genes_A_v_B_polyploid_CHG_context_1.txt")
HM_ALKG1_v_RS7G4_DMG_CHH_hal <- read.delim("ARPEGGIO_results_alk1Vsyn4_HM/dmrseq/CHH_context/DM_genes_A_v_B_polyploid_CHH_context_1.txt")

HM_TKSG1_v_halG1_DMG_CG <- read.delim("ARPEGGIO_results_tks1Vpro1_HM/dmrseq/CG_context/DM_genes_parent1_v_allo_CG_context.txt")
HM_TKSG1_v_halG1_DMG_CHG <- read.delim("ARPEGGIO_results_tks1Vpro1_HM/dmrseq/CHG_context/DM_genes_parent1_v_allo_CHG_context.txt")
HM_TKSG1_v_halG1_DMG_CHH <- read.delim("ARPEGGIO_results_tks1Vpro1_HM/dmrseq/CHH_context/DM_genes_parent1_v_allo_CHH_context.txt")

HM_TKSG1_v_RS7G1_DMG_CG_hal <- read.delim("ARPEGGIO_results_tks1Vsyn1_HM/dmrseq/CG_context/DM_genes_A_v_B_polyploid_CG_context_1.txt")
HM_TKSG1_v_RS7G1_DMG_CHG_hal <- read.delim("ARPEGGIO_results_tks1Vsyn1_HM/dmrseq/CHG_context/DM_genes_A_v_B_polyploid_CHG_context_1.txt")
HM_TKSG1_v_RS7G1_DMG_CHH_hal <- read.delim("ARPEGGIO_results_tks1Vsyn1_HM/dmrseq/CHH_context/DM_genes_A_v_B_polyploid_CHH_context_1.txt")

HM_TKSG1_v_RS7G4_DMG_CG_hal <- read.delim("ARPEGGIO_results_tks1Vsyn4_HM/dmrseq/CG_context/DM_genes_A_v_B_polyploid_CG_context_1.txt")
HM_TKSG1_v_RS7G4_DMG_CHG_hal <- read.delim("ARPEGGIO_results_tks1Vsyn4_HM/dmrseq/CHG_context/DM_genes_A_v_B_polyploid_CHG_context_1.txt")
HM_TKSG1_v_RS7G4_DMG_CHH_hal <- read.delim("ARPEGGIO_results_tks1Vsyn4_HM/dmrseq/CHH_context/DM_genes_A_v_B_polyploid_CHH_context_1.txt")


### lyrata side

HM_lyrG1_v_RS7G1_DMG_CG <- read.delim("ARPEGGIO_results_syn1Vpro1_HM/dmrseq/CG_context/DM_genes_parent2_v_allo_CG_context.txt")
HM_lyrG1_v_RS7G1_DMG_CHG <- read.delim("ARPEGGIO_results_syn1Vpro1_HM/dmrseq/CHG_context/DM_genes_parent2_v_allo_CHG_context.txt")
HM_lyrG1_v_RS7G1_DMG_CHH <- read.delim("ARPEGGIO_results_syn1Vpro1_HM/dmrseq/CHH_context/DM_genes_parent2_v_allo_CHH_context.txt")

HM_lyrG1_v_RS7G4_DMG_CG <- read.delim("ARPEGGIO_results_syn4Vpro1_HM/dmrseq/CG_context/DM_genes_parent2_v_allo_CG_context.txt")
HM_lyrG1_v_RS7G4_DMG_CHG <- read.delim("ARPEGGIO_results_syn4Vpro1_HM/dmrseq/CHG_context/DM_genes_parent2_v_allo_CHG_context.txt")
HM_lyrG1_v_RS7G4_DMG_CHH <- read.delim("ARPEGGIO_results_syn4Vpro1_HM/dmrseq/CHH_context/DM_genes_parent2_v_allo_CHH_context.txt")

HM_RS7G1_v_RS7G4_DMG_CG_lyr <- read.delim("ARPEGGIO_results_syn4v1_HM/dmrseq/CG_context/DM_genes_A_v_B_polyploid_CG_context_2.txt")
HM_RS7G1_v_RS7G4_DMG_CHG_lyr <- read.delim("ARPEGGIO_results_syn4v1_HM/dmrseq/CHG_context/DM_genes_A_v_B_polyploid_CHG_context_2.txt")
HM_RS7G1_v_RS7G4_DMG_CHH_lyr <- read.delim("ARPEGGIO_results_syn4v1_HM/dmrseq/CHH_context/DM_genes_A_v_B_polyploid_CHH_context_2.txt")

HM_ALKG1_v_lyrG1_DMG_CG <- read.delim("ARPEGGIO_results_alk1Vpro1_HM/dmrseq/CG_context/DM_genes_parent2_v_allo_CG_context.txt")
HM_ALKG1_v_lyrG1_DMG_CHG <- read.delim("ARPEGGIO_results_alk1Vpro1_HM/dmrseq/CHG_context/DM_genes_parent2_v_allo_CHG_context.txt")
HM_ALKG1_v_lyrG1_DMG_CHH <- read.delim("ARPEGGIO_results_alk1Vpro1_HM/dmrseq/CHH_context/DM_genes_parent2_v_allo_CHH_context.txt")

HM_ALKG1_v_RS7G1_DMG_CG_lyr <- read.delim("ARPEGGIO_results_alk1Vsyn1_HM/dmrseq/CG_context/DM_genes_A_v_B_polyploid_CG_context_2.txt")
HM_ALKG1_v_RS7G1_DMG_CHG_lyr <- read.delim("ARPEGGIO_results_alk1Vsyn1_HM/dmrseq/CHG_context/DM_genes_A_v_B_polyploid_CHG_context_2.txt")
HM_ALKG1_v_RS7G1_DMG_CHH_lyr <- read.delim("ARPEGGIO_results_alk1Vsyn1_HM/dmrseq/CHH_context/DM_genes_A_v_B_polyploid_CHH_context_2.txt")

HM_ALKG1_v_RS7G4_DMG_CG_lyr <- read.delim("ARPEGGIO_results_alk1Vsyn4_HM/dmrseq/CG_context/DM_genes_A_v_B_polyploid_CG_context_2.txt")
HM_ALKG1_v_RS7G4_DMG_CHG_lyr <- read.delim("ARPEGGIO_results_alk1Vsyn4_HM/dmrseq/CHG_context/DM_genes_A_v_B_polyploid_CHG_context_2.txt")
HM_ALKG1_v_RS7G4_DMG_CHH_lyr <- read.delim("ARPEGGIO_results_alk1Vsyn4_HM/dmrseq/CHH_context/DM_genes_A_v_B_polyploid_CHH_context_2.txt")

HM_TKSG1_v_lyrG1_DMG_CG <- read.delim("ARPEGGIO_results_tks1Vpro1_HM/dmrseq/CG_context/DM_genes_parent2_v_allo_CG_context.txt")
HM_TKSG1_v_lyrG1_DMG_CHG <- read.delim("ARPEGGIO_results_tks1Vpro1_HM/dmrseq/CHG_context/DM_genes_parent2_v_allo_CHG_context.txt")
HM_TKSG1_v_lyrG1_DMG_CHH <- read.delim("ARPEGGIO_results_tks1Vpro1_HM/dmrseq/CHH_context/DM_genes_parent2_v_allo_CHH_context.txt")

HM_TKSG1_v_RS7G1_DMG_CG_lyr <- read.delim("ARPEGGIO_results_tks1Vsyn1_HM/dmrseq/CG_context/DM_genes_A_v_B_polyploid_CG_context_2.txt")
HM_TKSG1_v_RS7G1_DMG_CHG_lyr <- read.delim("ARPEGGIO_results_tks1Vsyn1_HM/dmrseq/CHG_context/DM_genes_A_v_B_polyploid_CHG_context_2.txt")
HM_TKSG1_v_RS7G1_DMG_CHH_lyr <- read.delim("ARPEGGIO_results_tks1Vsyn1_HM/dmrseq/CHH_context/DM_genes_A_v_B_polyploid_CHH_context_2.txt")

HM_TKSG1_v_RS7G4_DMG_CG_lyr <- read.delim("ARPEGGIO_results_tks1Vsyn4_HM/dmrseq/CG_context/DM_genes_A_v_B_polyploid_CG_context_2.txt")
HM_TKSG1_v_RS7G4_DMG_CHG_lyr <- read.delim("ARPEGGIO_results_tks1Vsyn4_HM/dmrseq/CHG_context/DM_genes_A_v_B_polyploid_CHG_context_2.txt")
HM_TKSG1_v_RS7G4_DMG_CHH_lyr <- read.delim("ARPEGGIO_results_tks1Vsyn4_HM/dmrseq/CHH_context/DM_genes_A_v_B_polyploid_CHH_context_2.txt")


### Hot conditions (LL)

### halleri side

LL_halG1_v_RS7G1_DMG_CG <- read.delim("ARPEGGIO_results_syn1Vpro1_LL/dmrseq/CG_context/DM_genes_parent1_v_allo_CG_context.txt")
LL_halG1_v_RS7G1_DMG_CHG <- read.delim("ARPEGGIO_results_syn1Vpro1_LL/dmrseq/CHG_context/DM_genes_parent1_v_allo_CHG_context.txt")
LL_halG1_v_RS7G1_DMG_CHH <- read.delim("ARPEGGIO_results_syn1Vpro1_LL/dmrseq/CHH_context/DM_genes_parent1_v_allo_CHH_context.txt")

LL_halG1_v_RS7G4_DMG_CG <- read.delim("ARPEGGIO_results_syn4Vpro1_LL/dmrseq/CG_context/DM_genes_parent1_v_allo_CG_context.txt")
LL_halG1_v_RS7G4_DMG_CHG <- read.delim("ARPEGGIO_results_syn4Vpro1_LL/dmrseq/CHG_context/DM_genes_parent1_v_allo_CHG_context.txt")
LL_halG1_v_RS7G4_DMG_CHH <- read.delim("ARPEGGIO_results_syn4Vpro1_LL/dmrseq/CHH_context/DM_genes_parent1_v_allo_CHH_context.txt")

LL_RS7G1_v_RS7G4_DMG_CG_hal <- read.delim("ARPEGGIO_results_syn4v1_LL/dmrseq/CG_context/DM_genes_A_v_B_polyploid_CG_context_1.txt")
LL_RS7G1_v_RS7G4_DMG_CHG_hal <- read.delim("ARPEGGIO_results_syn4v1_LL/dmrseq/CHG_context/DM_genes_A_v_B_polyploid_CHG_context_1.txt")
LL_RS7G1_v_RS7G4_DMG_CHH_hal <- read.delim("ARPEGGIO_results_syn4v1_LL/dmrseq/CHH_context/DM_genes_A_v_B_polyploid_CHH_context_1.txt")

LL_ALKG1_v_halG1_DMG_CG <- read.delim("ARPEGGIO_results_alk1Vpro1_LL/dmrseq/CG_context/DM_genes_parent1_v_allo_CG_context.txt")
LL_ALKG1_v_halG1_DMG_CHG <- read.delim("ARPEGGIO_results_alk1Vpro1_LL/dmrseq/CHG_context/DM_genes_parent1_v_allo_CHG_context.txt")
LL_ALKG1_v_halG1_DMG_CHH <- read.delim("ARPEGGIO_results_alk1Vpro1_LL/dmrseq/CHH_context/DM_genes_parent1_v_allo_CHH_context.txt")

LL_ALKG1_v_RS7G1_DMG_CG_hal <- read.delim("ARPEGGIO_results_alk1Vsyn1_LL/dmrseq/CG_context/DM_genes_A_v_B_polyploid_CG_context_1.txt")
LL_ALKG1_v_RS7G1_DMG_CHG_hal <- read.delim("ARPEGGIO_results_alk1Vsyn1_LL/dmrseq/CHG_context/DM_genes_A_v_B_polyploid_CHG_context_1.txt")
LL_ALKG1_v_RS7G1_DMG_CHH_hal <- read.delim("ARPEGGIO_results_alk1Vsyn1_LL/dmrseq/CHH_context/DM_genes_A_v_B_polyploid_CHH_context_1.txt")

LL_ALKG1_v_RS7G4_DMG_CG_hal <- read.delim("ARPEGGIO_results_alk1Vsyn4_LL/dmrseq/CG_context/DM_genes_A_v_B_polyploid_CG_context_1.txt")
LL_ALKG1_v_RS7G4_DMG_CHG_hal <- read.delim("ARPEGGIO_results_alk1Vsyn4_LL/dmrseq/CHG_context/DM_genes_A_v_B_polyploid_CHG_context_1.txt")
LL_ALKG1_v_RS7G4_DMG_CHH_hal <- read.delim("ARPEGGIO_results_alk1Vsyn4_LL/dmrseq/CHH_context/DM_genes_A_v_B_polyploid_CHH_context_1.txt")

LL_TKSG1_v_halG1_DMG_CG <- read.delim("ARPEGGIO_results_tks1Vpro1_LL/dmrseq/CG_context/DM_genes_parent1_v_allo_CG_context.txt")
LL_TKSG1_v_halG1_DMG_CHG <- read.delim("ARPEGGIO_results_tks1Vpro1_LL/dmrseq/CHG_context/DM_genes_parent1_v_allo_CHG_context.txt")
LL_TKSG1_v_halG1_DMG_CHH <- read.delim("ARPEGGIO_results_tks1Vpro1_LL/dmrseq/CHH_context/DM_genes_parent1_v_allo_CHH_context.txt")

LL_TKSG1_v_RS7G1_DMG_CG_hal <- read.delim("ARPEGGIO_results_tks1Vsyn1_LL/dmrseq/CG_context/DM_genes_A_v_B_polyploid_CG_context_1.txt")
LL_TKSG1_v_RS7G1_DMG_CHG_hal <- read.delim("ARPEGGIO_results_tks1Vsyn1_LL/dmrseq/CHG_context/DM_genes_A_v_B_polyploid_CHG_context_1.txt")
LL_TKSG1_v_RS7G1_DMG_CHH_hal <- read.delim("ARPEGGIO_results_tks1Vsyn1_LL/dmrseq/CHH_context/DM_genes_A_v_B_polyploid_CHH_context_1.txt")

LL_TKSG1_v_RS7G4_DMG_CG_hal <- read.delim("ARPEGGIO_results_tks1Vsyn4_LL/dmrseq/CG_context/DM_genes_A_v_B_polyploid_CG_context_1.txt")
LL_TKSG1_v_RS7G4_DMG_CHG_hal <- read.delim("ARPEGGIO_results_tks1Vsyn4_LL/dmrseq/CHG_context/DM_genes_A_v_B_polyploid_CHG_context_1.txt")
LL_TKSG1_v_RS7G4_DMG_CHH_hal <- read.delim("ARPEGGIO_results_tks1Vsyn4_LL/dmrseq/CHH_context/DM_genes_A_v_B_polyploid_CHH_context_1.txt")


### lyrata side

LL_lyrG1_v_RS7G1_DMG_CG <- read.delim("ARPEGGIO_results_syn1Vpro1_LL/dmrseq/CG_context/DM_genes_parent2_v_allo_CG_context.txt")
LL_lyrG1_v_RS7G1_DMG_CHG <- read.delim("ARPEGGIO_results_syn1Vpro1_LL/dmrseq/CHG_context/DM_genes_parent2_v_allo_CHG_context.txt")
LL_lyrG1_v_RS7G1_DMG_CHH <- read.delim("ARPEGGIO_results_syn1Vpro1_LL/dmrseq/CHH_context/DM_genes_parent2_v_allo_CHH_context.txt")

LL_lyrG1_v_RS7G4_DMG_CG <- read.delim("ARPEGGIO_results_syn4Vpro1_LL/dmrseq/CG_context/DM_genes_parent2_v_allo_CG_context.txt")
LL_lyrG1_v_RS7G4_DMG_CHG <- read.delim("ARPEGGIO_results_syn4Vpro1_LL/dmrseq/CHG_context/DM_genes_parent2_v_allo_CHG_context.txt")
LL_lyrG1_v_RS7G4_DMG_CHH <- read.delim("ARPEGGIO_results_syn4Vpro1_LL/dmrseq/CHH_context/DM_genes_parent2_v_allo_CHH_context.txt")

LL_RS7G1_v_RS7G4_DMG_CG_lyr <- read.delim("ARPEGGIO_results_syn4v1_LL/dmrseq/CG_context/DM_genes_A_v_B_polyploid_CG_context_2.txt")
LL_RS7G1_v_RS7G4_DMG_CHG_lyr <- read.delim("ARPEGGIO_results_syn4v1_LL/dmrseq/CHG_context/DM_genes_A_v_B_polyploid_CHG_context_2.txt")
LL_RS7G1_v_RS7G4_DMG_CHH_lyr <- read.delim("ARPEGGIO_results_syn4v1_LL/dmrseq/CHH_context/DM_genes_A_v_B_polyploid_CHH_context_2.txt")

LL_ALKG1_v_lyrG1_DMG_CG <- read.delim("ARPEGGIO_results_alk1Vpro1_LL/dmrseq/CG_context/DM_genes_parent2_v_allo_CG_context.txt")
LL_ALKG1_v_lyrG1_DMG_CHG <- read.delim("ARPEGGIO_results_alk1Vpro1_LL/dmrseq/CHG_context/DM_genes_parent2_v_allo_CHG_context.txt")
LL_ALKG1_v_lyrG1_DMG_CHH <- read.delim("ARPEGGIO_results_alk1Vpro1_LL/dmrseq/CHH_context/DM_genes_parent2_v_allo_CHH_context.txt")

LL_ALKG1_v_RS7G1_DMG_CG_lyr <- read.delim("ARPEGGIO_results_alk1Vsyn1_LL/dmrseq/CG_context/DM_genes_A_v_B_polyploid_CG_context_2.txt")
LL_ALKG1_v_RS7G1_DMG_CHG_lyr <- read.delim("ARPEGGIO_results_alk1Vsyn1_LL/dmrseq/CHG_context/DM_genes_A_v_B_polyploid_CHG_context_2.txt")
LL_ALKG1_v_RS7G1_DMG_CHH_lyr <- read.delim("ARPEGGIO_results_alk1Vsyn1_LL/dmrseq/CHH_context/DM_genes_A_v_B_polyploid_CHH_context_2.txt")

LL_ALKG1_v_RS7G4_DMG_CG_lyr <- read.delim("ARPEGGIO_results_alk1Vsyn4_LL/dmrseq/CG_context/DM_genes_A_v_B_polyploid_CG_context_2.txt")
LL_ALKG1_v_RS7G4_DMG_CHG_lyr <- read.delim("ARPEGGIO_results_alk1Vsyn4_LL/dmrseq/CHG_context/DM_genes_A_v_B_polyploid_CHG_context_2.txt")
LL_ALKG1_v_RS7G4_DMG_CHH_lyr <- read.delim("ARPEGGIO_results_alk1Vsyn4_LL/dmrseq/CHH_context/DM_genes_A_v_B_polyploid_CHH_context_2.txt")

LL_TKSG1_v_lyrG1_DMG_CG <- read.delim("ARPEGGIO_results_tks1Vpro1_LL/dmrseq/CG_context/DM_genes_parent2_v_allo_CG_context.txt")
LL_TKSG1_v_lyrG1_DMG_CHG <- read.delim("ARPEGGIO_results_tks1Vpro1_LL/dmrseq/CHG_context/DM_genes_parent2_v_allo_CHG_context.txt")
LL_TKSG1_v_lyrG1_DMG_CHH <- read.delim("ARPEGGIO_results_tks1Vpro1_LL/dmrseq/CHH_context/DM_genes_parent2_v_allo_CHH_context.txt")

LL_TKSG1_v_RS7G1_DMG_CG_lyr <- read.delim("ARPEGGIO_results_tks1Vsyn1_LL/dmrseq/CG_context/DM_genes_A_v_B_polyploid_CG_context_2.txt")
LL_TKSG1_v_RS7G1_DMG_CHG_lyr <- read.delim("ARPEGGIO_results_tks1Vsyn1_LL/dmrseq/CHG_context/DM_genes_A_v_B_polyploid_CHG_context_2.txt")
LL_TKSG1_v_RS7G1_DMG_CHH_lyr <- read.delim("ARPEGGIO_results_tks1Vsyn1_LL/dmrseq/CHH_context/DM_genes_A_v_B_polyploid_CHH_context_2.txt")

LL_TKSG1_v_RS7G4_DMG_CG_lyr <- read.delim("ARPEGGIO_results_tks1Vsyn4_LL/dmrseq/CG_context/DM_genes_A_v_B_polyploid_CG_context_2.txt")
LL_TKSG1_v_RS7G4_DMG_CHG_lyr <- read.delim("ARPEGGIO_results_tks1Vsyn4_LL/dmrseq/CHG_context/DM_genes_A_v_B_polyploid_CHG_context_2.txt")
LL_TKSG1_v_RS7G4_DMG_CHH_lyr <- read.delim("ARPEGGIO_results_tks1Vsyn4_LL/dmrseq/CHH_context/DM_genes_A_v_B_polyploid_CHH_context_2.txt")

# We filter all DMGs that are part of low coverage regions

# Import file with low coverage regions

HM_hal_lowC_scaffolds  <- read.table("~/OneDrive/PhD/Project/Chapter_3/coverage_analysis/HM_RS7_G4_hal_LowCovScaffolds.txt", quote="\"", comment.char="")
HM_lyr_lowC_scaffolds  <- read.table("~/OneDrive/PhD/Project/Chapter_3/coverage_analysis/HM_RS7_G4_lyr_LowCovScaffolds.txt", quote="\"", comment.char="")

LL_hal_lowC_scaffolds  <- read.table("~/OneDrive/PhD/Project/Chapter_3/coverage_analysis/LL_RS7_G4_hal_LowCovScaffolds.txt", quote="\"", comment.char="")
LL_lyr_lowC_scaffolds  <- read.table("~/OneDrive/PhD/Project/Chapter_3/coverage_analysis/LL_RS7_G4_lyr_LowCovScaffolds.txt", quote="\"", comment.char="")

HM_halG1_v_RS7G1_DMG_CG <- filter(HM_halG1_v_RS7G1_DMG_CG, 
                                  !(seqname %in% HM_hal_lowC_scaffolds$V1))
HM_halG1_v_RS7G1_DMG_CHG <- filter(HM_halG1_v_RS7G1_DMG_CHG, 
                                   !(seqname %in% HM_hal_lowC_scaffolds$V1))
HM_halG1_v_RS7G1_DMG_CHH <- filter(HM_halG1_v_RS7G1_DMG_CHH, 
                                   !(seqname %in% HM_hal_lowC_scaffolds$V1))

HM_halG1_v_RS7G4_DMG_CG <- filter(HM_halG1_v_RS7G4_DMG_CG, 
                                  !(seqname %in% HM_hal_lowC_scaffolds$V1))
HM_halG1_v_RS7G4_DMG_CHG <- filter(HM_halG1_v_RS7G4_DMG_CHG, 
                                   !(seqname %in% HM_hal_lowC_scaffolds$V1))
HM_halG1_v_RS7G4_DMG_CHH <- filter(HM_halG1_v_RS7G4_DMG_CHH, 
                                   !(seqname %in% HM_hal_lowC_scaffolds$V1))

HM_RS7G1_v_RS7G4_DMG_CG_hal <- filter(HM_RS7G1_v_RS7G4_DMG_CG_hal, 
                                      !(seqname %in% HM_hal_lowC_scaffolds$V1))
HM_RS7G1_v_RS7G4_DMG_CHG_hal <- filter(HM_RS7G1_v_RS7G4_DMG_CHG_hal, 
                                       !(seqname %in% HM_hal_lowC_scaffolds$V1))
HM_RS7G1_v_RS7G4_DMG_CHH_hal <- filter(HM_RS7G1_v_RS7G4_DMG_CHH_hal, 
                                       !(seqname %in% HM_hal_lowC_scaffolds$V1))

HM_ALKG1_v_halG1_DMG_CG <- filter(HM_ALKG1_v_halG1_DMG_CG, 
                                  !(seqname %in% HM_hal_lowC_scaffolds$V1))
HM_ALKG1_v_halG1_DMG_CHG <- filter(HM_ALKG1_v_halG1_DMG_CHG, 
                                   !(seqname %in% HM_hal_lowC_scaffolds$V1))
HM_ALKG1_v_halG1_DMG_CHH <- filter(HM_ALKG1_v_halG1_DMG_CHH, 
                                   !(seqname %in% HM_hal_lowC_scaffolds$V1))

HM_ALKG1_v_RS7G1_DMG_CG_hal <- filter(HM_ALKG1_v_RS7G1_DMG_CG_hal, 
                                      !(seqname %in% HM_hal_lowC_scaffolds$V1))
HM_ALKG1_v_RS7G1_DMG_CHG_hal <- filter(HM_ALKG1_v_RS7G1_DMG_CHG_hal, 
                                       !(seqname %in% HM_hal_lowC_scaffolds$V1))
HM_ALKG1_v_RS7G1_DMG_CHH_hal <- filter(HM_ALKG1_v_RS7G1_DMG_CHH_hal, 
                                       !(seqname %in% HM_hal_lowC_scaffolds$V1))

HM_ALKG1_v_RS7G4_DMG_CG_hal <- filter(HM_ALKG1_v_RS7G4_DMG_CG_hal, 
                                      !(seqname %in% HM_hal_lowC_scaffolds$V1))
HM_ALKG1_v_RS7G4_DMG_CHG_hal <- filter(HM_ALKG1_v_RS7G4_DMG_CHG_hal, 
                                       !(seqname %in% HM_hal_lowC_scaffolds$V1))
HM_ALKG1_v_RS7G4_DMG_CHH_hal <- filter(HM_ALKG1_v_RS7G4_DMG_CHH_hal, 
                                       !(seqname %in% HM_hal_lowC_scaffolds$V1))

HM_TKSG1_v_halG1_DMG_CG <- filter(HM_TKSG1_v_halG1_DMG_CG, 
                                  !(seqname %in% HM_hal_lowC_scaffolds$V1))
HM_TKSG1_v_halG1_DMG_CHG <- filter(HM_TKSG1_v_halG1_DMG_CHG, 
                                   !(seqname %in% HM_hal_lowC_scaffolds$V1))
HM_TKSG1_v_halG1_DMG_CHH <- filter(HM_TKSG1_v_halG1_DMG_CHH, 
                                   !(seqname %in% HM_hal_lowC_scaffolds$V1))

HM_TKSG1_v_RS7G1_DMG_CG_hal <- filter(HM_TKSG1_v_RS7G1_DMG_CG_hal, 
                                      !(seqname %in% HM_hal_lowC_scaffolds$V1))
HM_TKSG1_v_RS7G1_DMG_CHG_hal <- filter(HM_TKSG1_v_RS7G1_DMG_CHG_hal, 
                                       !(seqname %in% HM_hal_lowC_scaffolds$V1))
HM_TKSG1_v_RS7G1_DMG_CHH_hal <- filter(HM_TKSG1_v_RS7G1_DMG_CHH_hal, 
                                       !(seqname %in% HM_hal_lowC_scaffolds$V1))

HM_TKSG1_v_RS7G4_DMG_CG_hal <- filter(HM_TKSG1_v_RS7G4_DMG_CG_hal, 
                                      !(seqname %in% HM_hal_lowC_scaffolds$V1))
HM_TKSG1_v_RS7G4_DMG_CHG_hal <- filter(HM_TKSG1_v_RS7G4_DMG_CHG_hal, 
                                       !(seqname %in% HM_hal_lowC_scaffolds$V1))
HM_TKSG1_v_RS7G4_DMG_CHH_hal <- filter(HM_TKSG1_v_RS7G4_DMG_CHH_hal, 
                                       !(seqname %in% HM_hal_lowC_scaffolds$V1))


### lyrata side

HM_lyrG1_v_RS7G1_DMG_CG <- filter(HM_lyrG1_v_RS7G1_DMG_CG, 
                                  !(seqname %in% HM_lyr_lowC_scaffolds$V1))
HM_lyrG1_v_RS7G1_DMG_CHG <- filter(HM_lyrG1_v_RS7G1_DMG_CHG, 
                                   !(seqname %in% HM_lyr_lowC_scaffolds$V1))
HM_lyrG1_v_RS7G1_DMG_CHH <- filter(HM_lyrG1_v_RS7G1_DMG_CHH, 
                                   !(seqname %in% HM_lyr_lowC_scaffolds$V1))

HM_lyrG1_v_RS7G4_DMG_CG <- filter(HM_lyrG1_v_RS7G4_DMG_CG, 
                                  !(seqname %in% HM_lyr_lowC_scaffolds$V1))
HM_lyrG1_v_RS7G4_DMG_CHG <- filter(HM_lyrG1_v_RS7G4_DMG_CHG, 
                                   !(seqname %in% HM_lyr_lowC_scaffolds$V1))
HM_lyrG1_v_RS7G4_DMG_CHH <- filter(HM_lyrG1_v_RS7G4_DMG_CHH, 
                                   !(seqname %in% HM_lyr_lowC_scaffolds$V1))

HM_RS7G1_v_RS7G4_DMG_CG_lyr <- filter(HM_RS7G1_v_RS7G4_DMG_CG_lyr, 
                                      !(seqname %in% HM_lyr_lowC_scaffolds$V1))
HM_RS7G1_v_RS7G4_DMG_CHG_lyr <- filter(HM_RS7G1_v_RS7G4_DMG_CHG_lyr, 
                                       !(seqname %in% HM_lyr_lowC_scaffolds$V1))
HM_RS7G1_v_RS7G4_DMG_CHH_lyr <- filter(HM_RS7G1_v_RS7G4_DMG_CHH_lyr, 
                                       !(seqname %in% HM_lyr_lowC_scaffolds$V1))

HM_ALKG1_v_lyrG1_DMG_CG <- filter(HM_ALKG1_v_lyrG1_DMG_CG, 
                                  !(seqname %in% HM_lyr_lowC_scaffolds$V1))
HM_ALKG1_v_lyrG1_DMG_CHG <- filter(HM_ALKG1_v_lyrG1_DMG_CHG, 
                                   !(seqname %in% HM_lyr_lowC_scaffolds$V1))
HM_ALKG1_v_lyrG1_DMG_CHH <- filter(HM_ALKG1_v_lyrG1_DMG_CHH, 
                                   !(seqname %in% HM_lyr_lowC_scaffolds$V1))

HM_ALKG1_v_RS7G1_DMG_CG_lyr <- filter(HM_ALKG1_v_RS7G1_DMG_CG_lyr, 
                                      !(seqname %in% HM_lyr_lowC_scaffolds$V1))
HM_ALKG1_v_RS7G1_DMG_CHG_lyr <- filter(HM_ALKG1_v_RS7G1_DMG_CHG_lyr, 
                                       !(seqname %in% HM_lyr_lowC_scaffolds$V1))
HM_ALKG1_v_RS7G1_DMG_CHH_lyr <- filter(HM_ALKG1_v_RS7G1_DMG_CHH_lyr, 
                                       !(seqname %in% HM_lyr_lowC_scaffolds$V1))

HM_ALKG1_v_RS7G4_DMG_CG_lyr <- filter(HM_ALKG1_v_RS7G4_DMG_CG_lyr, 
                                      !(seqname %in% HM_lyr_lowC_scaffolds$V1))
HM_ALKG1_v_RS7G4_DMG_CHG_lyr <- filter(HM_ALKG1_v_RS7G4_DMG_CHG_lyr, 
                                       !(seqname %in% HM_lyr_lowC_scaffolds$V1))
HM_ALKG1_v_RS7G4_DMG_CHH_lyr <- filter(HM_ALKG1_v_RS7G4_DMG_CHH_lyr, 
                                       !(seqname %in% HM_lyr_lowC_scaffolds$V1))

HM_TKSG1_v_lyrG1_DMG_CG <- filter(HM_TKSG1_v_lyrG1_DMG_CG, 
                                  !(seqname %in% HM_lyr_lowC_scaffolds$V1))
HM_TKSG1_v_lyrG1_DMG_CHG <- filter(HM_TKSG1_v_lyrG1_DMG_CHG, 
                                   !(seqname %in% HM_lyr_lowC_scaffolds$V1))
HM_TKSG1_v_lyrG1_DMG_CHH <- filter(HM_TKSG1_v_lyrG1_DMG_CHH, 
                                   !(seqname %in% HM_lyr_lowC_scaffolds$V1))

HM_TKSG1_v_RS7G1_DMG_CG_lyr <- filter(HM_TKSG1_v_RS7G1_DMG_CG_lyr, 
                                      !(seqname %in% HM_lyr_lowC_scaffolds$V1))
HM_TKSG1_v_RS7G1_DMG_CHG_lyr <- filter(HM_TKSG1_v_RS7G1_DMG_CHG_lyr, 
                                       !(seqname %in% HM_lyr_lowC_scaffolds$V1))
HM_TKSG1_v_RS7G1_DMG_CHH_lyr <- filter(HM_TKSG1_v_RS7G1_DMG_CHH_lyr, 
                                       !(seqname %in% HM_lyr_lowC_scaffolds$V1))

HM_TKSG1_v_RS7G4_DMG_CG_lyr <- filter(HM_TKSG1_v_RS7G4_DMG_CG_lyr, 
                                      !(seqname %in% HM_lyr_lowC_scaffolds$V1))
HM_TKSG1_v_RS7G4_DMG_CHG_lyr <- filter(HM_TKSG1_v_RS7G4_DMG_CHG_lyr, 
                                       !(seqname %in% HM_lyr_lowC_scaffolds$V1))
HM_TKSG1_v_RS7G4_DMG_CHH_lyr <- filter(HM_TKSG1_v_RS7G4_DMG_CHH_lyr, 
                                       !(seqname %in% HM_lyr_lowC_scaffolds$V1))


### Hot conditions (LL)

### halleri side

LL_halG1_v_RS7G1_DMG_CG <- filter(LL_halG1_v_RS7G1_DMG_CG, 
                                  !(seqname %in% LL_hal_lowC_scaffolds$V1))
LL_halG1_v_RS7G1_DMG_CHG <- filter(LL_halG1_v_RS7G1_DMG_CHG, 
                                   !(seqname %in% LL_hal_lowC_scaffolds$V1))
LL_halG1_v_RS7G1_DMG_CHH <- filter(LL_halG1_v_RS7G1_DMG_CHH, 
                                   !(seqname %in% LL_hal_lowC_scaffolds$V1))

LL_halG1_v_RS7G4_DMG_CG <- filter(LL_halG1_v_RS7G4_DMG_CG, 
                                  !(seqname %in% LL_hal_lowC_scaffolds$V1))
LL_halG1_v_RS7G4_DMG_CHG <- filter(LL_halG1_v_RS7G4_DMG_CHG, 
                                   !(seqname %in% LL_hal_lowC_scaffolds$V1))
LL_halG1_v_RS7G4_DMG_CHH <- filter(LL_halG1_v_RS7G4_DMG_CHH, 
                                   !(seqname %in% LL_hal_lowC_scaffolds$V1))

LL_RS7G1_v_RS7G4_DMG_CG_hal <- filter(LL_RS7G1_v_RS7G4_DMG_CG_hal, 
                                      !(seqname %in% LL_hal_lowC_scaffolds$V1))
LL_RS7G1_v_RS7G4_DMG_CHG_hal <- filter(LL_RS7G1_v_RS7G4_DMG_CHG_hal, 
                                       !(seqname %in% LL_hal_lowC_scaffolds$V1))
LL_RS7G1_v_RS7G4_DMG_CHH_hal <- filter(LL_RS7G1_v_RS7G4_DMG_CHH_hal, 
                                       !(seqname %in% LL_hal_lowC_scaffolds$V1))

LL_ALKG1_v_halG1_DMG_CG <- filter(LL_ALKG1_v_halG1_DMG_CG, 
                                  !(seqname %in% LL_hal_lowC_scaffolds$V1))
LL_ALKG1_v_halG1_DMG_CHG <- filter(LL_ALKG1_v_halG1_DMG_CHG, 
                                   !(seqname %in% LL_hal_lowC_scaffolds$V1))
LL_ALKG1_v_halG1_DMG_CHH <- filter(LL_ALKG1_v_halG1_DMG_CHH, 
                                   !(seqname %in% LL_hal_lowC_scaffolds$V1))

LL_ALKG1_v_RS7G1_DMG_CG_hal <- filter(LL_ALKG1_v_RS7G1_DMG_CG_hal, 
                                      !(seqname %in% LL_hal_lowC_scaffolds$V1))
LL_ALKG1_v_RS7G1_DMG_CHG_hal <- filter(LL_ALKG1_v_RS7G1_DMG_CHG_hal, 
                                       !(seqname %in% LL_hal_lowC_scaffolds$V1))
LL_ALKG1_v_RS7G1_DMG_CHH_hal <- filter(LL_ALKG1_v_RS7G1_DMG_CHH_hal, 
                                       !(seqname %in% LL_hal_lowC_scaffolds$V1))

LL_ALKG1_v_RS7G4_DMG_CG_hal <- filter(LL_ALKG1_v_RS7G4_DMG_CG_hal, 
                                      !(seqname %in% LL_hal_lowC_scaffolds$V1))
LL_ALKG1_v_RS7G4_DMG_CHG_hal <- filter(LL_ALKG1_v_RS7G4_DMG_CHG_hal, 
                                       !(seqname %in% LL_hal_lowC_scaffolds$V1))
LL_ALKG1_v_RS7G4_DMG_CHH_hal <- filter(LL_ALKG1_v_RS7G4_DMG_CHH_hal, 
                                       !(seqname %in% LL_hal_lowC_scaffolds$V1))

LL_TKSG1_v_halG1_DMG_CG <- filter(LL_TKSG1_v_halG1_DMG_CG, 
                                  !(seqname %in% LL_hal_lowC_scaffolds$V1))
LL_TKSG1_v_halG1_DMG_CHG <- filter(LL_TKSG1_v_halG1_DMG_CHG, 
                                   !(seqname %in% LL_hal_lowC_scaffolds$V1))
LL_TKSG1_v_halG1_DMG_CHH <- filter(LL_TKSG1_v_halG1_DMG_CHH, 
                                   !(seqname %in% LL_hal_lowC_scaffolds$V1))

LL_TKSG1_v_RS7G1_DMG_CG_hal <- filter(LL_TKSG1_v_RS7G1_DMG_CG_hal, 
                                      !(seqname %in% LL_hal_lowC_scaffolds$V1))
LL_TKSG1_v_RS7G1_DMG_CHG_hal <- filter(LL_TKSG1_v_RS7G1_DMG_CHG_hal, 
                                       !(seqname %in% LL_hal_lowC_scaffolds$V1))
LL_TKSG1_v_RS7G1_DMG_CHH_hal <- filter(LL_TKSG1_v_RS7G1_DMG_CHH_hal, 
                                       !(seqname %in% LL_hal_lowC_scaffolds$V1))

LL_TKSG1_v_RS7G4_DMG_CG_hal <- filter(LL_TKSG1_v_RS7G4_DMG_CG_hal, 
                                      !(seqname %in% LL_hal_lowC_scaffolds$V1))
LL_TKSG1_v_RS7G4_DMG_CHG_hal <- filter(LL_TKSG1_v_RS7G4_DMG_CHG_hal, 
                                       !(seqname %in% LL_hal_lowC_scaffolds$V1))
LL_TKSG1_v_RS7G4_DMG_CHH_hal <- filter(LL_TKSG1_v_RS7G4_DMG_CHH_hal, 
                                       !(seqname %in% LL_hal_lowC_scaffolds$V1))


### lyrata side

LL_lyrG1_v_RS7G1_DMG_CG <- filter(LL_lyrG1_v_RS7G1_DMG_CG, 
                                  !(seqname %in% LL_lyr_lowC_scaffolds$V1))
LL_lyrG1_v_RS7G1_DMG_CHG <- filter(LL_lyrG1_v_RS7G1_DMG_CHG, 
                                   !(seqname %in% LL_lyr_lowC_scaffolds$V1))
LL_lyrG1_v_RS7G1_DMG_CHH <- filter(LL_lyrG1_v_RS7G1_DMG_CHH, 
                                   !(seqname %in% LL_lyr_lowC_scaffolds$V1))

LL_lyrG1_v_RS7G4_DMG_CG <- filter(LL_lyrG1_v_RS7G4_DMG_CG, 
                                  !(seqname %in% LL_lyr_lowC_scaffolds$V1))
LL_lyrG1_v_RS7G4_DMG_CHG <- filter(LL_lyrG1_v_RS7G4_DMG_CHG, 
                                   !(seqname %in% LL_lyr_lowC_scaffolds$V1))
LL_lyrG1_v_RS7G4_DMG_CHH <- filter(LL_lyrG1_v_RS7G4_DMG_CHH, 
                                   !(seqname %in% LL_lyr_lowC_scaffolds$V1))

LL_RS7G1_v_RS7G4_DMG_CG_lyr <- filter(LL_RS7G1_v_RS7G4_DMG_CG_lyr, 
                                      !(seqname %in% LL_lyr_lowC_scaffolds$V1))
LL_RS7G1_v_RS7G4_DMG_CHG_lyr <- filter(LL_RS7G1_v_RS7G4_DMG_CHG_lyr, 
                                       !(seqname %in% LL_lyr_lowC_scaffolds$V1))
LL_RS7G1_v_RS7G4_DMG_CHH_lyr <- filter(LL_RS7G1_v_RS7G4_DMG_CHH_lyr, 
                                       !(seqname %in% LL_lyr_lowC_scaffolds$V1))

LL_ALKG1_v_lyrG1_DMG_CG <- filter(LL_ALKG1_v_lyrG1_DMG_CG, 
                                  !(seqname %in% LL_lyr_lowC_scaffolds$V1))
LL_ALKG1_v_lyrG1_DMG_CHG <- filter(LL_ALKG1_v_lyrG1_DMG_CHG, 
                                   !(seqname %in% LL_lyr_lowC_scaffolds$V1))
LL_ALKG1_v_lyrG1_DMG_CHH <- filter(LL_ALKG1_v_lyrG1_DMG_CHH, 
                                   !(seqname %in% LL_lyr_lowC_scaffolds$V1))

LL_ALKG1_v_RS7G1_DMG_CG_lyr <- filter(LL_ALKG1_v_RS7G1_DMG_CG_lyr, 
                                      !(seqname %in% LL_lyr_lowC_scaffolds$V1))
LL_ALKG1_v_RS7G1_DMG_CHG_lyr <- filter(LL_ALKG1_v_RS7G1_DMG_CHG_lyr, 
                                       !(seqname %in% LL_lyr_lowC_scaffolds$V1))
LL_ALKG1_v_RS7G1_DMG_CHH_lyr <- filter(LL_ALKG1_v_RS7G1_DMG_CHH_lyr, 
                                       !(seqname %in% LL_lyr_lowC_scaffolds$V1))

LL_ALKG1_v_RS7G4_DMG_CG_lyr <- filter(LL_ALKG1_v_RS7G4_DMG_CG_lyr, 
                                      !(seqname %in% LL_lyr_lowC_scaffolds$V1))
LL_ALKG1_v_RS7G4_DMG_CHG_lyr <- filter(LL_ALKG1_v_RS7G4_DMG_CHG_lyr, 
                                       !(seqname %in% LL_lyr_lowC_scaffolds$V1))
LL_ALKG1_v_RS7G4_DMG_CHH_lyr <- filter(LL_ALKG1_v_RS7G4_DMG_CHH_lyr, 
                                       !(seqname %in% LL_lyr_lowC_scaffolds$V1))

LL_TKSG1_v_lyrG1_DMG_CG <- filter(LL_TKSG1_v_lyrG1_DMG_CG, 
                                  !(seqname %in% LL_lyr_lowC_scaffolds$V1))
LL_TKSG1_v_lyrG1_DMG_CHG <- filter(LL_TKSG1_v_lyrG1_DMG_CHG, 
                                   !(seqname %in% LL_lyr_lowC_scaffolds$V1))
LL_TKSG1_v_lyrG1_DMG_CHH <- filter(LL_TKSG1_v_lyrG1_DMG_CHH, 
                                   !(seqname %in% LL_lyr_lowC_scaffolds$V1))

LL_TKSG1_v_RS7G1_DMG_CG_lyr <- filter(LL_TKSG1_v_RS7G1_DMG_CG_lyr, 
                                      !(seqname %in% LL_lyr_lowC_scaffolds$V1))
LL_TKSG1_v_RS7G1_DMG_CHG_lyr <- filter(LL_TKSG1_v_RS7G1_DMG_CHG_lyr, 
                                       !(seqname %in% LL_lyr_lowC_scaffolds$V1))
LL_TKSG1_v_RS7G1_DMG_CHH_lyr <- filter(LL_TKSG1_v_RS7G1_DMG_CHH_lyr, 
                                       !(seqname %in% LL_lyr_lowC_scaffolds$V1))

LL_TKSG1_v_RS7G4_DMG_CG_lyr <- filter(LL_TKSG1_v_RS7G4_DMG_CG_lyr, 
                                      !(seqname %in% LL_lyr_lowC_scaffolds$V1))
LL_TKSG1_v_RS7G4_DMG_CHG_lyr <- filter(LL_TKSG1_v_RS7G4_DMG_CHG_lyr, 
                                       !(seqname %in% LL_lyr_lowC_scaffolds$V1))
LL_TKSG1_v_RS7G4_DMG_CHH_lyr <- filter(LL_TKSG1_v_RS7G4_DMG_CHH_lyr, 
                                       !(seqname %in% LL_lyr_lowC_scaffolds$V1))

### Import expression data

setwd("/Users/ste/OneDrive/PhD/Project/Chapter_3/Reports/RNAseqDEG/")

# Cold conditions (HM)

### halleri side

HM_halG1_v_RS7G1_DEG <- read.delim("HM_halG1_v_RS7G1_DEG.txt")
HM_halG1_v_RS7G4_DEG <- read.delim("HM_halG1_v_RS7G4_DEG.txt")

HM_RS7G1_v_RS7G4_hal_DEG <- read.delim("HM_RS7G1_v_G4_hal_DEG.txt")

HM_ALKG1_v_halG1_DEG <- read.delim("HM_ALKG1_v_halG1_DEG.txt")
HM_ALKG1_v_RS7G1_hal_DEG <- read.delim("HM_ALKG1_v_RS7G1_hal_DEG.txt")
HM_ALKG1_v_RS7G4_hal_DEG <- read.delim("HM_ALKG1_v_RS7G4_hal_DEG.txt")

HM_TKSG1_v_halG1_DEG <- read.delim("HM_TKSG1_v_halG1_DEG.txt")
HM_TKSG1_v_RS7G1_hal_DEG <- read.delim("HM_TKSG1_v_RS7G1_hal_DEG.txt")
HM_TKSG1_v_RS7G4_hal_DEG <- read.delim("HM_TKSG1_v_RS7G4_hal_DEG.txt")


### lyrata side

HM_lyrG1_v_RS7G1_DEG <- read.delim("HM_lyrG1_v_RS7G1_DEG.txt")
HM_lyrG1_v_RS7G4_DEG <- read.delim("HM_lyrG1_v_RS7G4_DEG.txt")

HM_RS7G1_v_RS7G4_lyr_DEG <- read.delim("HM_RS7G1_v_G4_lyr_DEG.txt")

HM_ALKG1_v_lyrG1_DEG <- read.delim("HM_ALKG1_v_lyrG1_DEG.txt")
HM_ALKG1_v_RS7G1_lyr_DEG <- read.delim("HM_ALKG1_v_RS7G1_lyr_DEG.txt")
HM_ALKG1_v_RS7G4_lyr_DEG <- read.delim("HM_ALKG1_v_RS7G4_lyr_DEG.txt")

HM_TKSG1_v_lyrG1_DEG <- read.delim("HM_TKSG1_v_lyrG1_DEG.txt")
HM_TKSG1_v_RS7G1_lyr_DEG <- read.delim("HM_TKSG1_v_RS7G1_lyr_DEG.txt")
HM_TKSG1_v_RS7G4_lyr_DEG <- read.delim("HM_TKSG1_v_RS7G4_lyr_DEG.txt")


# Hot conditions (LL)

### halleri side

LL_halG1_v_RS7G1_DEG <- read.delim("LL_halG1_v_RS7G1_DEG.txt")
LL_halG1_v_RS7G4_DEG <- read.delim("LL_halG1_v_RS7G4_DEG.txt")

LL_RS7G1_v_RS7G4_hal_DEG <- read.delim("LL_RS7G1_v_G4_hal_DEG.txt")

LL_ALKG1_v_halG1_DEG <- read.delim("LL_ALKG1_v_halG1_DEG.txt")
LL_ALKG1_v_RS7G1_hal_DEG <- read.delim("LL_ALKG1_v_RS7G1_hal_DEG.txt")
LL_ALKG1_v_RS7G4_hal_DEG <- read.delim("LL_ALKG1_v_RS7G4_hal_DEG.txt")

LL_TKSG1_v_halG1_DEG <- read.delim("LL_TKSG1_v_halG1_DEG.txt")
LL_TKSG1_v_RS7G1_hal_DEG <- read.delim("LL_TKSG1_v_RS7G1_hal_DEG.txt")
LL_TKSG1_v_RS7G4_hal_DEG <- read.delim("LL_TKSG1_v_RS7G4_hal_DEG.txt")

### lyrata side

LL_lyrG1_v_RS7G1_DEG <- read.delim("LL_lyrG1_v_RS7G1_DEG.txt")
LL_lyrG1_v_RS7G4_DEG <- read.delim("LL_lyrG1_v_RS7G4_DEG.txt")

LL_RS7G1_v_RS7G4_lyr_DEG <- read.delim("LL_RS7G1_v_G4_lyr_DEG.txt")

LL_ALKG1_v_lyrG1_DEG <- read.delim("LL_ALKG1_v_lyrG1_DEG.txt")
LL_ALKG1_v_RS7G1_lyr_DEG <- read.delim("LL_ALKG1_v_RS7G1_lyr_DEG.txt")
LL_ALKG1_v_RS7G4_lyr_DEG <- read.delim("LL_ALKG1_v_RS7G4_lyr_DEG.txt")

LL_TKSG1_v_lyrG1_DEG <- read.delim("LL_TKSG1_v_lyrG1_DEG.txt")
LL_TKSG1_v_RS7G1_lyr_DEG <- read.delim("LL_TKSG1_v_RS7G1_lyr_DEG.txt")
LL_TKSG1_v_RS7G4_lyr_DEG <- read.delim("LL_TKSG1_v_RS7G4_lyr_DEG.txt")

# We filter DEG falling in low coverage regions as well

# We import annotations to find the corresponding scaffold of DEGs

Ahal_v2_2 <- read.delim("~/OneDrive/PhD/Project/bs_data/Ahal_genome/Ahal_v2_2.gff",
                        header=FALSE,
                        col.names = c("scaffold", "tool", "context", "start", "end", "number",
                                      "strand", "dot", "extra"))
Alyr_v2_2 <- read.delim("~/OneDrive/PhD/Project/bs_data/Alyr_genome/Alyr_v2_2_renamed.gff", 
                        header=FALSE,
                        col.names = c("scaffold", "tool", "context", "start", "end", "number",
                                      "strand", "dot", "extra"))

Ahal_anno <- select(Ahal_v2_2, scaffold, context, start, end, strand, extra)
rm(Ahal_v2_2)


Alyr_anno <- select(Alyr_v2_2, scaffold, context, start, end, strand, extra)
rm(Alyr_v2_2)

# We select only genes

Ahal_anno <- filter(Ahal_anno, context == "gene")
Alyr_anno <- filter(Alyr_anno, context == "gene")

# We pick the geneID from the extra column

split1 <- str_split_fixed(Ahal_anno$extra, "=", 2)[,2]
geneID_hal <- str_split_fixed(split1, ";", 2)[,1]

split1 <- str_split_fixed(Alyr_anno$extra, "=", 2)[,2]
geneID_lyr <- str_split_fixed(split1, ";", 2)[,1]

# We add the geneID to the annotation

Ahal_anno <- mutate(Ahal_anno, geneID = geneID_hal)
Alyr_anno <- mutate(Alyr_anno, geneID = geneID_lyr)

# We add scaffolds to all DEGs and then filter

add_scaffolds <- function(DEG, anno){
  
  scaffolds <- match(DEG$geneID, anno$geneID)
  DEG_new <- mutate(DEG, scaffold = anno$scaffold[scaffolds])
  return(DEG_new)
  
}

# Cold conditions (HM)

### halleri side

HM_halG1_v_RS7G1_DEG <- add_scaffolds(HM_halG1_v_RS7G1_DEG, Ahal_anno)
HM_halG1_v_RS7G4_DEG <- add_scaffolds(HM_halG1_v_RS7G4_DEG, Ahal_anno)

HM_RS7G1_v_RS7G4_hal_DEG <- add_scaffolds(HM_RS7G1_v_RS7G4_hal_DEG, Ahal_anno)

HM_ALKG1_v_halG1_DEG <- add_scaffolds(HM_ALKG1_v_halG1_DEG, Ahal_anno)
HM_ALKG1_v_RS7G1_hal_DEG <- add_scaffolds(HM_ALKG1_v_RS7G1_hal_DEG, Ahal_anno)
HM_ALKG1_v_RS7G4_hal_DEG <- add_scaffolds(HM_ALKG1_v_RS7G4_hal_DEG, Ahal_anno)

HM_TKSG1_v_halG1_DEG <- add_scaffolds(HM_TKSG1_v_halG1_DEG, Ahal_anno)
HM_TKSG1_v_RS7G1_hal_DEG <- add_scaffolds(HM_TKSG1_v_RS7G1_hal_DEG, Ahal_anno)
HM_TKSG1_v_RS7G4_hal_DEG <- add_scaffolds(HM_TKSG1_v_RS7G4_hal_DEG, Ahal_anno)


### lyrata side

HM_lyrG1_v_RS7G1_DEG <- add_scaffolds(HM_lyrG1_v_RS7G1_DEG, Alyr_anno)
HM_lyrG1_v_RS7G4_DEG <- add_scaffolds(HM_lyrG1_v_RS7G4_DEG, Alyr_anno)

HM_RS7G1_v_RS7G4_lyr_DEG <- add_scaffolds(HM_RS7G1_v_RS7G4_lyr_DEG, Alyr_anno)

HM_ALKG1_v_lyrG1_DEG <- add_scaffolds(HM_ALKG1_v_lyrG1_DEG, Alyr_anno)
HM_ALKG1_v_RS7G1_lyr_DEG <- add_scaffolds(HM_ALKG1_v_RS7G1_lyr_DEG, Alyr_anno)
HM_ALKG1_v_RS7G4_lyr_DEG <- add_scaffolds(HM_ALKG1_v_RS7G4_lyr_DEG, Alyr_anno)

HM_TKSG1_v_lyrG1_DEG <- add_scaffolds(HM_TKSG1_v_lyrG1_DEG, Alyr_anno)
HM_TKSG1_v_RS7G1_lyr_DEG <- add_scaffolds(HM_TKSG1_v_RS7G1_lyr_DEG, Alyr_anno)
HM_TKSG1_v_RS7G4_lyr_DEG <- add_scaffolds(HM_TKSG1_v_RS7G4_lyr_DEG, Alyr_anno)


# Hot conditions (LL)

### halleri side

LL_halG1_v_RS7G1_DEG <- add_scaffolds(LL_halG1_v_RS7G1_DEG, Ahal_anno)
LL_halG1_v_RS7G4_DEG <- add_scaffolds(LL_halG1_v_RS7G4_DEG, Ahal_anno)

LL_RS7G1_v_RS7G4_hal_DEG <- add_scaffolds(LL_RS7G1_v_RS7G4_hal_DEG, Ahal_anno)

LL_ALKG1_v_halG1_DEG <- add_scaffolds(LL_ALKG1_v_halG1_DEG, Ahal_anno)
LL_ALKG1_v_RS7G1_hal_DEG <- add_scaffolds(LL_ALKG1_v_RS7G1_hal_DEG, Ahal_anno)
LL_ALKG1_v_RS7G4_hal_DEG <- add_scaffolds(LL_ALKG1_v_RS7G4_hal_DEG, Ahal_anno)

LL_TKSG1_v_halG1_DEG <- add_scaffolds(LL_TKSG1_v_halG1_DEG, Ahal_anno)
LL_TKSG1_v_RS7G1_hal_DEG <- add_scaffolds(LL_TKSG1_v_RS7G1_hal_DEG, Ahal_anno)
LL_TKSG1_v_RS7G4_hal_DEG <- add_scaffolds(LL_TKSG1_v_RS7G4_hal_DEG, Ahal_anno)

### lyrata side

LL_lyrG1_v_RS7G1_DEG <- add_scaffolds(LL_lyrG1_v_RS7G1_DEG, Alyr_anno)
LL_lyrG1_v_RS7G4_DEG <- add_scaffolds(LL_lyrG1_v_RS7G4_DEG, Alyr_anno)

LL_RS7G1_v_RS7G4_lyr_DEG <- add_scaffolds(LL_RS7G1_v_RS7G4_lyr_DEG, Alyr_anno)

LL_ALKG1_v_lyrG1_DEG <- add_scaffolds(LL_ALKG1_v_lyrG1_DEG, Alyr_anno)
LL_ALKG1_v_RS7G1_lyr_DEG <- add_scaffolds(LL_ALKG1_v_RS7G1_lyr_DEG, Alyr_anno)
LL_ALKG1_v_RS7G4_lyr_DEG <- add_scaffolds(LL_ALKG1_v_RS7G4_lyr_DEG, Alyr_anno)

LL_TKSG1_v_lyrG1_DEG <- add_scaffolds(LL_TKSG1_v_lyrG1_DEG, Alyr_anno)
LL_TKSG1_v_RS7G1_lyr_DEG <- add_scaffolds(LL_TKSG1_v_RS7G1_lyr_DEG, Alyr_anno)
LL_TKSG1_v_RS7G4_lyr_DEG <- add_scaffolds(LL_TKSG1_v_RS7G4_lyr_DEG, Alyr_anno)

# Filtering

# Cold conditions (HM)

### halleri side

HM_halG1_v_RS7G1_DEG <- filter(HM_halG1_v_RS7G1_DEG, 
                               !(scaffold %in% HM_hal_lowC_scaffolds$V1))
HM_halG1_v_RS7G4_DEG <- filter(HM_halG1_v_RS7G4_DEG, 
                               !(scaffold %in% HM_hal_lowC_scaffolds$V1))

HM_RS7G1_v_RS7G4_hal_DEG <- filter(HM_RS7G1_v_RS7G4_hal_DEG, 
                                   !(scaffold %in% HM_hal_lowC_scaffolds$V1))

HM_ALKG1_v_halG1_DEG <- filter(HM_ALKG1_v_halG1_DEG, 
                               !(scaffold %in% HM_hal_lowC_scaffolds$V1))
HM_ALKG1_v_RS7G1_hal_DEG <- filter(HM_ALKG1_v_RS7G1_hal_DEG, 
                                   !(scaffold %in% HM_hal_lowC_scaffolds$V1))
HM_ALKG1_v_RS7G4_hal_DEG <- filter(HM_ALKG1_v_RS7G4_hal_DEG, 
                                   !(scaffold %in% HM_hal_lowC_scaffolds$V1))

HM_TKSG1_v_halG1_DEG <- filter(HM_TKSG1_v_halG1_DEG, 
                               !(scaffold %in% HM_hal_lowC_scaffolds$V1))
HM_TKSG1_v_RS7G1_hal_DEG <- filter(HM_TKSG1_v_RS7G1_hal_DEG, 
                                   !(scaffold %in% HM_hal_lowC_scaffolds$V1))
HM_TKSG1_v_RS7G4_hal_DEG <- filter(HM_TKSG1_v_RS7G4_hal_DEG, 
                                   !(scaffold %in% HM_hal_lowC_scaffolds$V1))


### lyrata side

HM_lyrG1_v_RS7G1_DEG <- filter(HM_lyrG1_v_RS7G1_DEG, 
                               !(scaffold %in% HM_lyr_lowC_scaffolds$V1))
HM_lyrG1_v_RS7G4_DEG <- filter(HM_lyrG1_v_RS7G4_DEG, 
                               !(scaffold %in% HM_lyr_lowC_scaffolds$V1))

HM_RS7G1_v_RS7G4_lyr_DEG <- filter(HM_RS7G1_v_RS7G4_lyr_DEG, 
                                   !(scaffold %in% HM_lyr_lowC_scaffolds$V1))

HM_ALKG1_v_lyrG1_DEG <- filter(HM_ALKG1_v_lyrG1_DEG, 
                               !(scaffold %in% HM_lyr_lowC_scaffolds$V1))
HM_ALKG1_v_RS7G1_lyr_DEG <- filter(HM_ALKG1_v_RS7G1_lyr_DEG, 
                                   !(scaffold %in% HM_lyr_lowC_scaffolds$V1))
HM_ALKG1_v_RS7G4_lyr_DEG <- filter(HM_ALKG1_v_RS7G4_lyr_DEG, 
                                   !(scaffold %in% HM_lyr_lowC_scaffolds$V1))

HM_TKSG1_v_lyrG1_DEG <- filter(HM_TKSG1_v_lyrG1_DEG, 
                               !(scaffold %in% HM_lyr_lowC_scaffolds$V1))
HM_TKSG1_v_RS7G1_lyr_DEG <- filter(HM_TKSG1_v_RS7G1_lyr_DEG, 
                                   !(scaffold %in% HM_lyr_lowC_scaffolds$V1))
HM_TKSG1_v_RS7G4_lyr_DEG <- filter(HM_TKSG1_v_RS7G4_lyr_DEG, 
                                   !(scaffold %in% HM_lyr_lowC_scaffolds$V1))


# Hot conditions (LL)

### halleri side

LL_halG1_v_RS7G1_DEG <- filter(LL_halG1_v_RS7G1_DEG, 
                               !(scaffold %in% LL_hal_lowC_scaffolds$V1))
LL_halG1_v_RS7G4_DEG <- filter(LL_halG1_v_RS7G4_DEG, 
                               !(scaffold %in% LL_hal_lowC_scaffolds$V1))

LL_RS7G1_v_RS7G4_hal_DEG <- filter(LL_RS7G1_v_RS7G4_hal_DEG, 
                                   !(scaffold %in% LL_hal_lowC_scaffolds$V1))

LL_ALKG1_v_halG1_DEG <- filter(LL_ALKG1_v_halG1_DEG, 
                               !(scaffold %in% LL_hal_lowC_scaffolds$V1))
LL_ALKG1_v_RS7G1_hal_DEG <- filter(LL_ALKG1_v_RS7G1_hal_DEG, 
                                   !(scaffold %in% LL_hal_lowC_scaffolds$V1))
LL_ALKG1_v_RS7G4_hal_DEG <- filter(LL_ALKG1_v_RS7G4_hal_DEG, 
                                   !(scaffold %in% LL_hal_lowC_scaffolds$V1))

LL_TKSG1_v_halG1_DEG <- filter(LL_TKSG1_v_halG1_DEG, 
                               !(scaffold %in% LL_hal_lowC_scaffolds$V1))
LL_TKSG1_v_RS7G1_hal_DEG <- filter(LL_TKSG1_v_RS7G1_hal_DEG, 
                                   !(scaffold %in% LL_hal_lowC_scaffolds$V1))
LL_TKSG1_v_RS7G4_hal_DEG <- filter(LL_TKSG1_v_RS7G4_hal_DEG, 
                                   !(scaffold %in% LL_hal_lowC_scaffolds$V1))

### lyrata side

LL_lyrG1_v_RS7G1_DEG <- filter(LL_lyrG1_v_RS7G1_DEG, 
                               !(scaffold %in% LL_lyr_lowC_scaffolds$V1))
LL_lyrG1_v_RS7G4_DEG <- filter(LL_lyrG1_v_RS7G4_DEG, 
                               !(scaffold %in% LL_lyr_lowC_scaffolds$V1))

LL_RS7G1_v_RS7G4_lyr_DEG <- filter(LL_RS7G1_v_RS7G4_lyr_DEG, 
                                   !(scaffold %in% LL_lyr_lowC_scaffolds$V1))

LL_ALKG1_v_lyrG1_DEG <- filter(LL_ALKG1_v_lyrG1_DEG, 
                               !(scaffold %in% LL_lyr_lowC_scaffolds$V1))
LL_ALKG1_v_RS7G1_lyr_DEG <- filter(LL_ALKG1_v_RS7G1_lyr_DEG, 
                                   !(scaffold %in% LL_lyr_lowC_scaffolds$V1))
LL_ALKG1_v_RS7G4_lyr_DEG <- filter(LL_ALKG1_v_RS7G4_lyr_DEG, 
                                   !(scaffold %in% LL_lyr_lowC_scaffolds$V1))

LL_TKSG1_v_lyrG1_DEG <- filter(LL_TKSG1_v_lyrG1_DEG, 
                               !(scaffold %in% LL_lyr_lowC_scaffolds$V1))
LL_TKSG1_v_RS7G1_lyr_DEG <- filter(LL_TKSG1_v_RS7G1_lyr_DEG, 
                                   !(scaffold %in% LL_lyr_lowC_scaffolds$V1))
LL_TKSG1_v_RS7G4_lyr_DEG <- filter(LL_TKSG1_v_RS7G4_lyr_DEG, 
                                   !(scaffold %in% LL_lyr_lowC_scaffolds$V1))


### Consolidate methylation data
### 1) merge all contexts
### 2) remove duplicated genes

consolidate <- function(DMG1, DMG2, DMG3){
  new_DMG <- rbind(DMG1, DMG2, DMG3)
  new_DMG <- arrange(new_DMG, geneID)
  new_DMG <- new_DMG[!duplicated(new_DMG$geneID),]
  return(new_DMG)
}

### Cold conditions

## halleri side

HM_halG1_v_RS7G1_DMG <- consolidate(HM_halG1_v_RS7G1_DMG_CG,
                                    HM_halG1_v_RS7G1_DMG_CHG,
                                    HM_halG1_v_RS7G1_DMG_CHH)

HM_halG1_v_RS7G4_DMG <- consolidate(HM_halG1_v_RS7G4_DMG_CG,
                                    HM_halG1_v_RS7G4_DMG_CHG,
                                    HM_halG1_v_RS7G4_DMG_CHH)

HM_RS7G1_v_RS7G4_DMG_hal <- consolidate(HM_RS7G1_v_RS7G4_DMG_CG_hal,
                                    HM_RS7G1_v_RS7G4_DMG_CHG_hal,
                                    HM_RS7G1_v_RS7G4_DMG_CHH_hal)

HM_ALKG1_v_halG1_DMG <- consolidate(HM_ALKG1_v_halG1_DMG_CG,
                                    HM_ALKG1_v_halG1_DMG_CHG,
                                    HM_ALKG1_v_halG1_DMG_CHH)

HM_ALKG1_v_RS7G1_DMG_hal <- consolidate(HM_ALKG1_v_RS7G1_DMG_CG_hal,
                                    HM_ALKG1_v_RS7G1_DMG_CHG_hal,
                                    HM_ALKG1_v_RS7G1_DMG_CHH_hal)

HM_ALKG1_v_RS7G4_DMG_hal <- consolidate(HM_ALKG1_v_RS7G4_DMG_CG_hal,
                                    HM_ALKG1_v_RS7G4_DMG_CHG_hal,
                                    HM_ALKG1_v_RS7G4_DMG_CHH_hal)

HM_TKSG1_v_halG1_DMG <- consolidate(HM_TKSG1_v_halG1_DMG_CG,
                                    HM_TKSG1_v_halG1_DMG_CHG,
                                    HM_TKSG1_v_halG1_DMG_CHH)

HM_TKSG1_v_RS7G1_DMG_hal <- consolidate(HM_TKSG1_v_RS7G1_DMG_CG_hal,
                                    HM_TKSG1_v_RS7G1_DMG_CHG_hal,
                                    HM_TKSG1_v_RS7G1_DMG_CHH_hal)

HM_TKSG1_v_RS7G4_DMG_hal <- consolidate(HM_TKSG1_v_RS7G4_DMG_CG_hal,
                                    HM_TKSG1_v_RS7G4_DMG_CHG_hal,
                                    HM_TKSG1_v_RS7G4_DMG_CHH_hal)

## lyrata side

HM_lyrG1_v_RS7G1_DMG <- consolidate(HM_lyrG1_v_RS7G1_DMG_CG,
                                    HM_lyrG1_v_RS7G1_DMG_CHG,
                                    HM_lyrG1_v_RS7G1_DMG_CHH)

HM_lyrG1_v_RS7G4_DMG <- consolidate(HM_lyrG1_v_RS7G4_DMG_CG,
                                    HM_lyrG1_v_RS7G4_DMG_CHG,
                                    HM_lyrG1_v_RS7G4_DMG_CHH)

HM_RS7G1_v_RS7G4_DMG_lyr <- consolidate(HM_RS7G1_v_RS7G4_DMG_CG_lyr,
                                        HM_RS7G1_v_RS7G4_DMG_CHG_lyr,
                                        HM_RS7G1_v_RS7G4_DMG_CHH_lyr)

HM_ALKG1_v_lyrG1_DMG <- consolidate(HM_ALKG1_v_lyrG1_DMG_CG,
                                    HM_ALKG1_v_lyrG1_DMG_CHG,
                                    HM_ALKG1_v_lyrG1_DMG_CHH)

HM_ALKG1_v_RS7G1_DMG_lyr <- consolidate(HM_ALKG1_v_RS7G1_DMG_CG_lyr,
                                        HM_ALKG1_v_RS7G1_DMG_CHG_lyr,
                                        HM_ALKG1_v_RS7G1_DMG_CHH_lyr)

HM_ALKG1_v_RS7G4_DMG_lyr <- consolidate(HM_ALKG1_v_RS7G4_DMG_CG_lyr,
                                        HM_ALKG1_v_RS7G4_DMG_CHG_lyr,
                                        HM_ALKG1_v_RS7G4_DMG_CHH_lyr)

HM_TKSG1_v_lyrG1_DMG <- consolidate(HM_TKSG1_v_lyrG1_DMG_CG,
                                    HM_TKSG1_v_lyrG1_DMG_CHG,
                                    HM_TKSG1_v_lyrG1_DMG_CHH)

HM_TKSG1_v_RS7G1_DMG_lyr <- consolidate(HM_TKSG1_v_RS7G1_DMG_CG_lyr,
                                        HM_TKSG1_v_RS7G1_DMG_CHG_lyr,
                                        HM_TKSG1_v_RS7G1_DMG_CHH_lyr)

HM_TKSG1_v_RS7G4_DMG_lyr <- consolidate(HM_TKSG1_v_RS7G4_DMG_CG_lyr,
                                        HM_TKSG1_v_RS7G4_DMG_CHG_lyr,
                                        HM_TKSG1_v_RS7G4_DMG_CHH_lyr)


### Hot conditions

## halleri side

LL_halG1_v_RS7G1_DMG <- consolidate(LL_halG1_v_RS7G1_DMG_CG,
                                    LL_halG1_v_RS7G1_DMG_CHG,
                                    LL_halG1_v_RS7G1_DMG_CHH)

LL_halG1_v_RS7G4_DMG <- consolidate(LL_halG1_v_RS7G4_DMG_CG,
                                    LL_halG1_v_RS7G4_DMG_CHG,
                                    LL_halG1_v_RS7G4_DMG_CHH)

LL_RS7G1_v_RS7G4_DMG_hal <- consolidate(LL_RS7G1_v_RS7G4_DMG_CG_hal,
                                        LL_RS7G1_v_RS7G4_DMG_CHG_hal,
                                        LL_RS7G1_v_RS7G4_DMG_CHH_hal)

LL_ALKG1_v_halG1_DMG <- consolidate(LL_ALKG1_v_halG1_DMG_CG,
                                    LL_ALKG1_v_halG1_DMG_CHG,
                                    LL_ALKG1_v_halG1_DMG_CHH)

LL_ALKG1_v_RS7G1_DMG_hal <- consolidate(LL_ALKG1_v_RS7G1_DMG_CG_hal,
                                        LL_ALKG1_v_RS7G1_DMG_CHG_hal,
                                        LL_ALKG1_v_RS7G1_DMG_CHH_hal)

LL_ALKG1_v_RS7G4_DMG_hal <- consolidate(LL_ALKG1_v_RS7G4_DMG_CG_hal,
                                        LL_ALKG1_v_RS7G4_DMG_CHG_hal,
                                        LL_ALKG1_v_RS7G4_DMG_CHH_hal)

LL_TKSG1_v_halG1_DMG <- consolidate(LL_TKSG1_v_halG1_DMG_CG,
                                    LL_TKSG1_v_halG1_DMG_CHG,
                                    LL_TKSG1_v_halG1_DMG_CHH)

LL_TKSG1_v_RS7G1_DMG_hal <- consolidate(LL_TKSG1_v_RS7G1_DMG_CG_hal,
                                        LL_TKSG1_v_RS7G1_DMG_CHG_hal,
                                        LL_TKSG1_v_RS7G1_DMG_CHH_hal)

LL_TKSG1_v_RS7G4_DMG_hal <- consolidate(LL_TKSG1_v_RS7G4_DMG_CG_hal,
                                        LL_TKSG1_v_RS7G4_DMG_CHG_hal,
                                        LL_TKSG1_v_RS7G4_DMG_CHH_hal)

## lyrata side

LL_lyrG1_v_RS7G1_DMG <- consolidate(LL_lyrG1_v_RS7G1_DMG_CG,
                                    LL_lyrG1_v_RS7G1_DMG_CHG,
                                    LL_lyrG1_v_RS7G1_DMG_CHH)

LL_lyrG1_v_RS7G4_DMG <- consolidate(LL_lyrG1_v_RS7G4_DMG_CG,
                                    LL_lyrG1_v_RS7G4_DMG_CHG,
                                    LL_lyrG1_v_RS7G4_DMG_CHH)

LL_RS7G1_v_RS7G4_DMG_lyr <- consolidate(LL_RS7G1_v_RS7G4_DMG_CG_lyr,
                                        LL_RS7G1_v_RS7G4_DMG_CHG_lyr,
                                        LL_RS7G1_v_RS7G4_DMG_CHH_lyr)

LL_ALKG1_v_lyrG1_DMG <- consolidate(LL_ALKG1_v_lyrG1_DMG_CG,
                                    LL_ALKG1_v_lyrG1_DMG_CHG,
                                    LL_ALKG1_v_lyrG1_DMG_CHH)

LL_ALKG1_v_RS7G1_DMG_lyr <- consolidate(LL_ALKG1_v_RS7G1_DMG_CG_lyr,
                                        LL_ALKG1_v_RS7G1_DMG_CHG_lyr,
                                        LL_ALKG1_v_RS7G1_DMG_CHH_lyr)

LL_ALKG1_v_RS7G4_DMG_lyr <- consolidate(LL_ALKG1_v_RS7G4_DMG_CG_lyr,
                                        LL_ALKG1_v_RS7G4_DMG_CHG_lyr,
                                        LL_ALKG1_v_RS7G4_DMG_CHH_lyr)

LL_TKSG1_v_lyrG1_DMG <- consolidate(LL_TKSG1_v_lyrG1_DMG_CG,
                                    LL_TKSG1_v_lyrG1_DMG_CHG,
                                    LL_TKSG1_v_lyrG1_DMG_CHH)

LL_TKSG1_v_RS7G1_DMG_lyr <- consolidate(LL_TKSG1_v_RS7G1_DMG_CG_lyr,
                                        LL_TKSG1_v_RS7G1_DMG_CHG_lyr,
                                        LL_TKSG1_v_RS7G1_DMG_CHH_lyr)

LL_TKSG1_v_RS7G4_DMG_lyr <- consolidate(LL_TKSG1_v_RS7G4_DMG_CG_lyr,
                                        LL_TKSG1_v_RS7G4_DMG_CHG_lyr,
                                        LL_TKSG1_v_RS7G4_DMG_CHH_lyr)

### Plot overlaps

setwd("~/OneDrive/PhD/Project/Chapter_3/Pictures/Paper_picsV3")

# Prepare a palette of 4 colors with R colorbrewer:
myCol <- brewer.pal(4, "Paired")

# Function to output venn diagram for two sets
venn <- function(set1, set2, title, filename){
  venn.diagram(
  x = list(set1, set2),
  category.names = c("DEG" , "DMG"),
  cat.cex = 0.7,
  cat.pos = 1,
  filename = filename,
  output=TRUE,
  main = title,
  main.cex = 0.5,
  print.mode = c("raw", "percent"),

  # Output features
  imagetype="png" ,
  height = 480 ,
  width = 480 ,
  resolution = 400,
  compression = "lzw",

  # Circles
  lwd = 2,
  lty = 'blank',
  fill = myCol[1:2],

  # Numbers
  cex = 0.5,
  fontface = "bold"
)
}

# Function to output venn diagram for four sets

venn4 <- function(set1, set2, set3, set4, title, filename){
  venn.diagram(
    x = list(set1, set2, set3, set4),
    category.names = c("DEG" , "DMG", "DEG" , "DMG"),
    cat.cex = 0.7,
    cat.pos = 1,
    filename = filename,
    output=TRUE,
    main = title,
    main.cex = 0.5,
    print.mode = c("raw", "percent"),
    
    # Output features
    imagetype="svg" ,
    height = 480 , 
    width = 480 , 
    resolution = 400,
    compression = "lzw",
    
    # Circles
    lwd = 2,
    lty = 'blank',
    fill = myCol[1:4],
    
    # Numbers
    cex = 0.5,
    fontface = "bold"
  )
}

# venn_DMG <- function(set1, set2, title, filename){
#   venn.diagram(
#     x = list(set1, set2),
#     category.names = c("DMG lyr" , "DMG hal"),
#     cat.cex = 0.7,
#     cat.pos = 1,
#     filename = filename,
#     output=TRUE,
#     main = title,
#     main.cex = 0.5,
#     
#     # Output features
#     imagetype="png" ,
#     height = 480 , 
#     width = 480 , 
#     resolution = 400,
#     compression = "lzw",
#     
#     # Circles
#     lwd = 2,
#     lty = 'blank',
#     fill = myCol[1:2],
#     
#     # Numbers
#     cex = 0.5,
#     fontface = "bold"
#   )
# }

## HM halleri G1 vs synthetic G1

venn_HM_hal1Vsyn1_DEG_v_DMG <- venn(set1 = HM_halG1_v_RS7G1_DEG$geneID,
                                    set2 = HM_halG1_v_RS7G1_DMG$geneID,
                                    title = "HM_halG1_v_RS7G1",
                                    filename = NULL)

ggsave(venn_HM_hal1Vsyn1_DEG_v_DMG, 
       filename = "HM_hal1Vsyn1_DEG_v_DMG.svg",
       device = "svg",
       width = 1.5,
       height = 1.5)

## HM halleri DMG vs lyrata DEG G1 vs synthetic G1

venn_HM_lyrDEGhalDMG1Vsyn1_DEG_v_DMG <- venn(set1 = HM_lyrG1_v_RS7G1_DEG$geneID,
                                             set2 = HM_halG1_v_RS7G1_DMG$geneID,
                                             title = "HM_lyrDEGhalDMGG1_v_RS7G1",
                                             filename = NULL)

ggsave(venn_HM_lyrDEGhalDMG1Vsyn1_DEG_v_DMG, 
       filename = "HM_lyrDEGhalDMG1Vsyn1_DEG_v_DMG.svg",
       device = "svg",
       width = 1.5,
       height = 1.5)

## HM lyrata DMG vs halleri DEG G1 vs synthetic G1

venn_HM_lyrDMGhalDEG1Vsyn1_DEG_v_DMG <- venn(set1 = HM_lyrG1_v_RS7G1_DMG$geneID,
     set2 = HM_halG1_v_RS7G1_DEG$geneID,
     title = "HM_lyrDMGhalDEGG1_v_RS7G1",
     filename = NULL)

ggsave(venn_HM_lyrDMGhalDEG1Vsyn1_DEG_v_DMG, 
       filename = "HM_lyrDEGhalDMG1Vsyn1_DEG_v_DMG.svg",
       device = "svg",
       width = 1.5,
       height = 1.5)

## HM halleri DMG vs lyrata DMG (G1)

venn_HM_lyrDMGhalDMG1Vsyn1 <- venn_DMG(set1 = HM_lyrG1_v_RS7G1_DMG$geneID,
     set2 = HM_halG1_v_RS7G1_DMG$geneID,
     title = "HM_lyrhalDMGG1_v_RS7G1",
     filename = NULL)

ggsave(venn_HM_lyrDMGhalDMG1Vsyn1, 
       filename = "HM_lyrDMGhalDMG1Vsyn1.svg",
       device = "svg",
       width = 1.5,
       height = 1.5)

## HM halleri DMG vs lyrata DMG (G4)

venn_HM_lyrDMGhalDMG1Vsyn4 <- venn_DMG(set1 = HM_lyrG1_v_RS7G4_DMG$geneID,
     set2 = HM_halG1_v_RS7G4_DMG$geneID,
     title = "HM_lyrhalDMGG4_v_RS7G1",
     filename = NULL)

ggsave(venn_HM_lyrDMGhalDMG1Vsyn4, 
       filename = "HM_lyrDMGhalDMG1Vsyn4.svg",
       device = "svg",
       width = 1.5,
       height = 1.5)

## LL halleri DMG vs lyrata DMG (G1)

venn_LL_lyrDMGhalDMG1Vsyn1 <- venn_DMG(set1 = LL_lyrG1_v_RS7G1_DMG$geneID,
         set2 = LL_halG1_v_RS7G1_DMG$geneID,
         title = "LL_lyrhalDMGG1_v_RS7G1",
         filename = NULL)

ggsave(venn_LL_lyrDMGhalDMG1Vsyn1, 
       filename = "LL_lyrDMGhalDMG1Vsyn1.svg",
       device = "svg",
       width = 1.5,
       height = 1.5)

## LL halleri DMG vs lyrata DMG (G4)

venn_LL_lyrDMGhalDMG1Vsyn4 <- venn_DMG(set1 = LL_lyrG1_v_RS7G4_DMG$geneID,
         set2 = LL_halG1_v_RS7G4_DMG$geneID,
         title = "LL_lyrhalDMGG4_v_RS7G1",
         filename = NULL)

ggsave(venn_LL_lyrDMGhalDMG1Vsyn4, 
       filename = "LL_lyrDMGhalDMG1Vsyn4.svg",
       device = "svg",
       width = 1.5,
       height = 1.5)

## HM halleri G1 vs synthetic G4

venn_HM_hal1Vsyn4_DEG_v_DMG <- venn(set1 = HM_halG1_v_RS7G4_DEG$geneID,
     set2 = HM_halG1_v_RS7G4_DMG$geneID,
     title = "HM_halG1_v_RS7G4",
     filename = NULL)

ggsave(venn_HM_hal1Vsyn4_DEG_v_DMG, 
       filename = "HM_hal1Vsyn4_DEG_v_DMG.svg",
       device = "svg",
       width = 1.5,
       height = 1.5)

## HM lyrata G1 vs synthetic G1

venn_HM_lyr1Vsyn1_DEG_v_DMG <- venn(set1 = HM_lyrG1_v_RS7G1_DEG$geneID,
     set2 = HM_lyrG1_v_RS7G1_DMG$geneID,
     title = "HM_lyrG1_v_RS7G1",
     filename = NULL)

ggsave(venn_HM_lyr1Vsyn1_DEG_v_DMG, 
       filename = "HM_lyr1Vsyn1_DEG_v_DMG.svg",
       device = "svg",
       width = 1.5,
       height = 1.5)

## HM lyrata G1 vs synthetic G4

venn_HM_lyr1Vsyn4_DEG_v_DMG <- venn(set1 = HM_lyrG1_v_RS7G4_DEG$geneID,
     set2 = HM_lyrG1_v_RS7G4_DMG$geneID,
     title = "HM_lyrG1_v_RS7G4",
     filename = NULL)

ggsave(venn_HM_lyr1Vsyn4_DEG_v_DMG, 
       filename = "HM_lyr1Vsyn4_DEG_v_DMG.svg",
       device = "svg",
       width = 1.5,
       height = 1.5)

## LL halleri G1 vs synthetic G1

venn_LL_hal1Vsyn1_DEG_v_DMG <- venn(set1 = LL_halG1_v_RS7G1_DEG$geneID,
     set2 = LL_halG1_v_RS7G1_DMG$geneID,
     title = "LL_halG1_v_RS7G1",
     filename = NULL)

ggsave(venn_LL_hal1Vsyn1_DEG_v_DMG, 
       filename = "LL_hal1Vsyn1_DEG_v_DMG.svg",
       device = "svg",
       width = 1.5,
       height = 1.5)

## LL halleri G1 vs synthetic G4

venn_LL_hal1Vsyn4_DEG_v_DMG <- venn(set1 = LL_halG1_v_RS7G4_DEG$geneID,
     set2 = LL_halG1_v_RS7G4_DMG$geneID,
     title = "LL_halG1_v_RS7G4",
     filename = NULL)

ggsave(venn_LL_hal1Vsyn4_DEG_v_DMG, 
       filename = "LL_hal1Vsyn4_DEG_v_DMG.svg",
       device = "svg",
       width = 1.5,
       height = 1.5)

## LL lyrata G1 vs synthetic G1

venn_LL_lyr1Vsyn1_DEG_v_DMG <- venn(set1 = LL_lyrG1_v_RS7G1_DEG$geneID,
     set2 = LL_lyrG1_v_RS7G1_DMG$geneID,
     title = "LL_lyrG1_v_RS7G1",
     filename = NULL)

ggsave(venn_LL_lyr1Vsyn1_DEG_v_DMG, 
       filename = "LL_lyr1Vsyn1_DEG_v_DMG.svg",
       device = "svg",
       width = 1.5,
       height = 1.5)

## LL lyrata G1 vs synthetic G4

venn_LL_lyr1Vsyn4_DEG_v_DMG <- venn(set1 = LL_lyrG1_v_RS7G4_DEG$geneID,
     set2 = LL_lyrG1_v_RS7G4_DMG$geneID,
     title = "LL_lyrG1_v_RS7G4",
     filename = NULL)

ggsave(venn_LL_lyr1Vsyn4_DEG_v_DMG, 
       filename = "LL_lyr1Vsyn4_DEG_v_DMG.svg",
       device = "svg",
       width = 1.5,
       height = 1.5)

## HM synthetic G1 vs synthetic G4

venn_HM_syn1Vsyn4_DEG_v_DMG <- venn(set1 = HM_RS7G1_v_RS7G4_hal_DEG$geneID,
     set2 = HM_RS7G1_v_RS7G4_DMG_hal$geneID,
     title = "HM_RS7G1_v_RS7G4 (hal)",
     filename = NULL)

ggsave(venn_HM_syn1Vsyn4_DEG_v_DMG, 
       filename = "HM_syn1Vsyn4_DEG_v_DMG.svg",
       device = "svg",
       width = 1.5,
       height = 1.5)

## HM ALK G1 vs pro G1 (hal)

venn_HM_alk1Vhal1_DEG_v_DMG <- venn(set1 = HM_ALKG1_v_halG1_DEG$geneID,
     set2 = HM_ALKG1_v_halG1_DMG$geneID,
     title = "HM_ALKG1_v_halG1",
     filename = NULL)

ggsave(venn_HM_alk1Vhal1_DEG_v_DMG, 
       filename = "HM_alk1Vhal1_DEG_v_DMG.svvg",
       device = "svg",
       width = 1.5,
       height = 1.5)

## HM ALK G1 vs synthetic G1 (hal)

venn_HM_alk1Vsyn1_DEG_v_DMG_hal <- venn(set1 = HM_ALKG1_v_RS7G1_hal_DEG$geneID,
     set2 = HM_ALKG1_v_RS7G1_DMG_hal$geneID,
     title = "HM_ALKG1_v_RS7G1 (hal)",
     filename = NULL)

ggsave(venn_HM_alk1Vsyn1_DEG_v_DMG_hal, 
       filename = "HM_alk1Vsyn1_DEG_v_DMG_hal.svg",
       device = "svg",
       width = 1.5,
       height = 1.5)

## HM ALK G1 vs synthetic G4 (hal)

venn_HM_alk1Vsyn4_DEG_v_DMG_hal <- venn(set1 = HM_ALKG1_v_RS7G4_hal_DEG$geneID,
     set2 = HM_ALKG1_v_RS7G4_DMG_hal$geneID,
     title = "HM_ALKG1_v_RS7G4 (hal)",
     filename = NULL)

ggsave(venn_HM_alk1Vsyn4_DEG_v_DMG_hal, 
       filename = "HM_alk1Vsyn4_DEG_v_DMG_hal.svg",
       device = "svg",
       width = 1.5,
       height = 1.5)

## HM TKS G1 vs pro G1 (hal)

venn_HM_tks1Vhal1_DEG_v_DMG <- venn(set1 = HM_TKSG1_v_halG1_DEG$geneID,
     set2 = HM_TKSG1_v_halG1_DMG$geneID,
     title = "HM_TKSG1_v_halG1",
     filename = NULL)

ggsave(venn_HM_tks1Vhal1_DEG_v_DMG, 
       filename = "HM_tks1Vhal1_DEG_v_DMG.svg",
       device = "svg",
       width = 1.5,
       height = 1.5)

## HM TKS G1 vs synthetic G1 (hal)

venn_HM_tks1Vsyn1_DEG_v_DMG_hal <- venn(set1 = HM_TKSG1_v_RS7G1_hal_DEG$geneID,
     set2 = HM_TKSG1_v_RS7G1_DMG_hal$geneID,
     title = "HM_TKSG1_v_RS7G1 (hal)",
     filename = NULL)

ggsave(venn_HM_tks1Vsyn1_DEG_v_DMG_hal, 
       filename = "HM_tks1Vsyn1_DEG_v_DMG_hal.svg",
       device = "svg",
       width = 1.5,
       height = 1.5)

## HM TKS G1 vs synthetic G4 (hal)

venn_HM_tks1Vsyn4_DEG_v_DMG_hal <- venn(set1 = HM_TKSG1_v_RS7G4_hal_DEG$geneID,
     set2 = HM_TKSG1_v_RS7G4_DMG_hal$geneID,
     title = "HM_TKSG1_v_RS7G4 (hal)",
     filename = NULL)

ggsave(venn_HM_tks1Vsyn4_DEG_v_DMG_hal, 
       filename = "HM_tks1Vsyn4_DEG_v_DMG_hal.svg",
       device = "svg",
       width = 1.5,
       height = 1.5)

## LL ALK G1 vs pro G1 (hal)

venn_LL_alk1Vhal1_DEG_v_DMG <- venn(set1 = LL_ALKG1_v_halG1_DEG$geneID,
     set2 = LL_ALKG1_v_halG1_DMG$geneID,
     title = "LL_ALKG1_v_halG1",
     filename = NULL)

ggsave(venn_LL_alk1Vhal1_DEG_v_DMG, 
       filename = "LL_alk1Vhal1_DEG_v_DMG.svg",
       device = "svg",
       width = 1.5,
       height = 1.5)

## LL ALK G1 vs synthetic G1 (hal)

venn_LL_alk1Vsyn1_DEG_v_DMG_hal <- venn(set1 = LL_ALKG1_v_RS7G1_hal_DEG$geneID,
     set2 = LL_ALKG1_v_RS7G1_DMG_hal$geneID,
     title = "LL_ALKG1_v_RS7G1 (hal)",
     filename = NULL)

ggsave(venn_LL_alk1Vsyn1_DEG_v_DMG_hal, 
       filename = "LL_alk1Vsyn1_DEG_v_DMG_hal.svg",
       device = "svg",
       width = 1.5,
       height = 1.5)

## LL ALK G1 vs synthetic G4 (hal)

venn_LL_alk1Vsyn4_DEG_v_DMG_hal <- venn(set1 = LL_ALKG1_v_RS7G4_hal_DEG$geneID,
     set2 = LL_ALKG1_v_RS7G4_DMG_hal$geneID,
     title = "LL_ALKG1_v_RS7G4 (hal)",
     filename = NULL)

ggsave(venn_LL_alk1Vsyn4_DEG_v_DMG_hal, 
       filename = "LL_alk1Vsyn4_DEG_v_DMG_hal.svg",
       device = "svg",
       width = 1.5,
       height = 1.5)

## LL TKS G1 vs pro G1 (hal)

venn_LL_tks1Vhal1_DEG_v_DMG <- venn(set1 = LL_TKSG1_v_halG1_DEG$geneID,
     set2 = LL_TKSG1_v_halG1_DMG$geneID,
     title = "LL_TKSG1_v_halG1",
     filename = NULL)

ggsave(venn_LL_tks1Vhal1_DEG_v_DMG, 
       filename = "LL_tks1Vhal1_DEG_v_DMG.svg",
       device = "svg",
       width = 1.5,
       height = 1.5)

## LL TKS G1 vs synthetic G1 (hal)

venn_LL_tks1Vsyn1_DEG_v_DMG_hal <- venn(set1 = LL_TKSG1_v_RS7G1_hal_DEG$geneID,
     set2 = LL_TKSG1_v_RS7G1_DMG_hal$geneID,
     title = "LL_TKSG1_v_RS7G1 (hal)",
     filename = NULL)

ggsave(venn_LL_tks1Vsyn1_DEG_v_DMG_hal, 
       filename = "LL_tks1Vsyn1_DEG_v_DMG_hal.svg",
       device = "svg",
       width = 1.5,
       height = 1.5)

## LL TKS G1 vs synthetic G4 (hal)

venn_LL_tks1Vsyn4_DEG_v_DMG_hal <- venn(set1 = LL_TKSG1_v_RS7G4_hal_DEG$geneID,
     set2 = LL_TKSG1_v_RS7G4_DMG_hal$geneID,
     title = "LL_TKSG1_v_RS7G4 (hal)",
     filename = NULL)

ggsave(venn_LL_tks1Vsyn4_DEG_v_DMG_hal, 
       filename = "LL_tks1Vsyn4_DEG_v_DMG_hal.svg",
       device = "svg",
       width = 1.5,
       height = 1.5)

## HM ALK G1 vs pro G1 (lyr)

venn_HM_alk1Vlyr1_DEG_v_DMG <- venn(set1 = HM_ALKG1_v_lyrG1_DEG$geneID,
     set2 = HM_ALKG1_v_lyrG1_DMG$geneID,
     title = "HM_ALKG1_v_lyrG1",
     filename = NULL)

ggsave(venn_HM_alk1Vlyr1_DEG_v_DMG, 
       filename = "HM_alk1Vlyr1_DEG_v_DMG.svg",
       device = "svg",
       width = 1.5,
       height = 1.5)

## HM ALK G1 vs synthetic G1 (lyr)

venn_HM_alk1Vsyn1_DEG_v_DMG_lyr <- venn(set1 = HM_ALKG1_v_RS7G1_lyr_DEG$geneID,
     set2 = HM_ALKG1_v_RS7G1_DMG_lyr$geneID,
     title = "HM_ALKG1_v_RS7G1 (lyr)",
     filename = NULL)

ggsave(venn_HM_alk1Vsyn1_DEG_v_DMG_lyr, 
       filename = "HM_alk1Vsyn1_DEG_v_DMG_lyr.svg",
       device = "svg",
       width = 1.5,
       height = 1.5)

## HM ALK G1 vs synthetic G4 (lyr)

venn_HM_alk1Vsyn4_DEG_v_DMG_lyr <- venn(set1 = HM_ALKG1_v_RS7G4_lyr_DEG$geneID,
     set2 = HM_ALKG1_v_RS7G4_DMG_lyr$geneID,
     title = "HM_ALKG1_v_RS7G4 (lyr)",
     filename = NULL)

ggsave(venn_HM_alk1Vsyn4_DEG_v_DMG_lyr, 
       filename = "HM_alk1Vsyn4_DEG_v_DMG_lyr.svg",
       device = "svg",
       width = 1.5,
       height = 1.5)

## HM TKS G1 vs pro G1 (lyr)

venn_HM_tks1Vlyr1_DEG_v_DMG <- venn(set1 = HM_TKSG1_v_lyrG1_DEG$geneID,
     set2 = HM_TKSG1_v_lyrG1_DMG$geneID,
     title = "HM_TKSG1_v_lyrG1",
     filename = NULL)

ggsave(venn_HM_tks1Vlyr1_DEG_v_DMG, 
       filename = "HM_tks1Vlyr1_DEG_v_DMG.svg",
       device = "svg",
       width = 1.5,
       height = 1.5)

## HM TKS G1 vs synthetic G1 (lyr)

venn_HM_tks1Vsyn1_DEG_v_DMG_lyr <- venn(set1 = HM_TKSG1_v_RS7G1_lyr_DEG$geneID,
     set2 = HM_TKSG1_v_RS7G1_DMG_lyr$geneID,
     title = "HM_TKSG1_v_RS7G1 (lyr)",
     filename = NULL)

ggsave(venn_HM_tks1Vsyn1_DEG_v_DMG_lyr, 
       filename = "HM_tks1Vsyn1_DEG_v_DMG_lyr.svg",
       device = "svg",
       width = 1.5,
       height = 1.5)

## HM TKS G1 vs synthetic G4 (lyr)

venn_HM_tks1Vsyn4_DEG_v_DMG_lyr <- venn(set1 = HM_TKSG1_v_RS7G4_lyr_DEG$geneID,
     set2 = HM_TKSG1_v_RS7G4_DMG_lyr$geneID,
     title = "HM_TKSG1_v_RS7G4 (lyr)",
     filename = NULL)

ggsave(venn_HM_tks1Vsyn4_DEG_v_DMG_lyr, 
       filename = "HM_tks1Vsyn4_DEG_v_DMG_lyr.svg",
       device = "svg",
       width = 1.5,
       height = 1.5)

## LL ALK G1 vs pro G1 (lyr)

venn_LL_alk1Vlyr1_DEG_v_DMG <- venn(set1 = LL_ALKG1_v_lyrG1_DEG$geneID,
     set2 = LL_ALKG1_v_lyrG1_DMG$geneID,
     title = "LL_ALKG1_v_lyrG1",
     filename = NULL)

ggsave(venn_LL_alk1Vlyr1_DEG_v_DMG, 
       filename = "LL_alk1Vlyr1_DEG_v_DMG.svg",
       device = "svg",
       width = 1.5,
       height = 1.5)

## LL ALK G1 vs synthetic G1 (lyr)

venn_LL_alk1Vsyn1_DEG_v_DMG_lyr <- venn(set1 = LL_ALKG1_v_RS7G1_lyr_DEG$geneID,
     set2 = LL_ALKG1_v_RS7G1_DMG_lyr$geneID,
     title = "LL_ALKG1_v_RS7G1 (lyr)",
     filename = NULL)

ggsave(venn_LL_alk1Vsyn1_DEG_v_DMG_lyr, 
       filename = "LL_alk1Vsyn1_DEG_v_DMG_lyr.svg",
       device = "svg",
       width = 1.5,
       height = 1.5)

## LL ALK G1 vs synthetic G4 (lyr)

venn_LL_alk1Vsyn4_DEG_v_DMG_lyr <- venn(set1 = LL_ALKG1_v_RS7G4_lyr_DEG$geneID,
     set2 = LL_ALKG1_v_RS7G4_DMG_lyr$geneID,
     title = "LL_ALKG1_v_RS7G4 (lyr)",
     filename = NULL)

ggsave(venn_LL_alk1Vsyn4_DEG_v_DMG_lyr, 
       filename = "LL_alk1Vsyn4_DEG_v_DMG_lyr.svg",
       device = "svg",
       width = 1.5,
       height = 1.5)

## LL TKS G1 vs pro G1 (lyr)

venn_LL_tks1Vlyr1_DEG_v_DMG <- venn(set1 = LL_TKSG1_v_lyrG1_DEG$geneID,
     set2 = LL_TKSG1_v_lyrG1_DMG$geneID,
     title = "LL_TKSG1_v_lyrG1",
     filename = NULL)

ggsave(venn_LL_tks1Vlyr1_DEG_v_DMG, 
       filename = "LL_tks1Vlyr1_DEG_v_DMG.svg",
       device = "svg",
       width = 1.5,
       height = 1.5)

## LL TKS G1 vs synthetic G1 (lyr)

venn_LL_tks1Vsyn1_DEG_v_DMG_lyr <- venn(set1 = LL_TKSG1_v_RS7G1_lyr_DEG$geneID,
     set2 = LL_TKSG1_v_RS7G1_DMG_lyr$geneID,
     title = "LL_TKSG1_v_RS7G1 (lyr)",
     filename = NULL)

ggsave(venn_LL_tks1Vsyn1_DEG_v_DMG_lyr, 
       filename = "LL_tks1Vsyn1_DEG_v_DMG_lyr.svg",
       device = "svg",
       width = 1.5,
       height = 1.5)

## LL TKS G1 vs synthetic G4 (lyr)

venn_LL_tks1Vsyn4_DEG_v_DMG_lyr <- venn(set1 = LL_TKSG1_v_RS7G4_lyr_DEG$geneID,
     set2 = LL_TKSG1_v_RS7G4_DMG_lyr$geneID,
     title = "LL_TKSG1_v_RS7G4 (lyr)",
     filename = NULL)

ggsave(venn_LL_tks1Vsyn4_DEG_v_DMG_lyr, 
       filename = "LL_tks1Vsyn4_DEG_v_DMG_lyr.svg",
       device = "svg",
       width = 1.5,
       height = 1.5)
