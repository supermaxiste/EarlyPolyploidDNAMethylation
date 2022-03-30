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

setwd("files")

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

HM_hal_lowC_scaffolds  <- read.table("files/coverage_data/HM_RS7_G4_hal_LowCovScaffolds.txt", quote="\"", comment.char="")
HM_lyr_lowC_scaffolds  <- read.table("files/coverage_data/HM_RS7_G4_lyr_LowCovScaffolds.txt", quote="\"", comment.char="")

LL_hal_lowC_scaffolds  <- read.table("files/coverage_data/LL_RS7_G4_hal_LowCovScaffolds.txt", quote="\"", comment.char="")
LL_lyr_lowC_scaffolds  <- read.table("files/coverage_data/LL_RS7_G4_lyr_LowCovScaffolds.txt", quote="\"", comment.char="")

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


### Import expression data

setwd("DEGs/")

# Cold conditions (HM)

### halleri side

HM_halG1_v_RS7G1_DEG <- read.delim("HM_halG1_v_RS7G1_DEG.txt")
HM_halG1_v_RS7G4_DEG <- read.delim("HM_halG1_v_RS7G4_DEG.txt")


### lyrata side

HM_lyrG1_v_RS7G1_DEG <- read.delim("HM_lyrG1_v_RS7G1_DEG.txt")
HM_lyrG1_v_RS7G4_DEG <- read.delim("HM_lyrG1_v_RS7G4_DEG.txt")


# Hot conditions (LL)

### halleri side

LL_halG1_v_RS7G1_DEG <- read.delim("LL_halG1_v_RS7G1_DEG.txt")
LL_halG1_v_RS7G4_DEG <- read.delim("LL_halG1_v_RS7G4_DEG.txt")

### lyrata side

LL_lyrG1_v_RS7G1_DEG <- read.delim("LL_lyrG1_v_RS7G1_DEG.txt")
LL_lyrG1_v_RS7G4_DEG <- read.delim("LL_lyrG1_v_RS7G4_DEG.txt")

# We filter DEG falling in low coverage regions as well

# We import annotations to find the corresponding scaffold of DEGs

Ahal_v2_2 <- read.delim("Ahal_v2_2.gff",
                        header=FALSE,
                        col.names = c("scaffold", "tool", "context", "start", "end", "number",
                                      "strand", "dot", "extra"))
Alyr_v2_2 <- read.delim("Alyr_v2_2_renamed.gff",
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

### lyrata side

HM_lyrG1_v_RS7G1_DEG <- add_scaffolds(HM_lyrG1_v_RS7G1_DEG, Alyr_anno)
HM_lyrG1_v_RS7G4_DEG <- add_scaffolds(HM_lyrG1_v_RS7G4_DEG, Alyr_anno)

# Hot conditions (LL)

### halleri side

LL_halG1_v_RS7G1_DEG <- add_scaffolds(LL_halG1_v_RS7G1_DEG, Ahal_anno)
LL_halG1_v_RS7G4_DEG <- add_scaffolds(LL_halG1_v_RS7G4_DEG, Ahal_anno)

### lyrata side

LL_lyrG1_v_RS7G1_DEG <- add_scaffolds(LL_lyrG1_v_RS7G1_DEG, Alyr_anno)
LL_lyrG1_v_RS7G4_DEG <- add_scaffolds(LL_lyrG1_v_RS7G4_DEG, Alyr_anno)

# Filtering

# Cold conditions (HM)

### halleri side

HM_halG1_v_RS7G1_DEG <- filter(HM_halG1_v_RS7G1_DEG,
                               !(scaffold %in% HM_hal_lowC_scaffolds$V1))
HM_halG1_v_RS7G4_DEG <- filter(HM_halG1_v_RS7G4_DEG,
                               !(scaffold %in% HM_hal_lowC_scaffolds$V1))

### lyrata side

HM_lyrG1_v_RS7G1_DEG <- filter(HM_lyrG1_v_RS7G1_DEG,
                               !(scaffold %in% HM_lyr_lowC_scaffolds$V1))
HM_lyrG1_v_RS7G4_DEG <- filter(HM_lyrG1_v_RS7G4_DEG,
                               !(scaffold %in% HM_lyr_lowC_scaffolds$V1))

# Hot conditions (LL)

### halleri side

LL_halG1_v_RS7G1_DEG <- filter(LL_halG1_v_RS7G1_DEG,
                               !(scaffold %in% LL_hal_lowC_scaffolds$V1))
LL_halG1_v_RS7G4_DEG <- filter(LL_halG1_v_RS7G4_DEG,
                               !(scaffold %in% LL_hal_lowC_scaffolds$V1))

### lyrata side

LL_lyrG1_v_RS7G1_DEG <- filter(LL_lyrG1_v_RS7G1_DEG,
                               !(scaffold %in% LL_lyr_lowC_scaffolds$V1))
LL_lyrG1_v_RS7G4_DEG <- filter(LL_lyrG1_v_RS7G4_DEG,
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


## lyrata side

HM_lyrG1_v_RS7G1_DMG <- consolidate(HM_lyrG1_v_RS7G1_DMG_CG,
                                    HM_lyrG1_v_RS7G1_DMG_CHG,
                                    HM_lyrG1_v_RS7G1_DMG_CHH)

HM_lyrG1_v_RS7G4_DMG <- consolidate(HM_lyrG1_v_RS7G4_DMG_CG,
                                    HM_lyrG1_v_RS7G4_DMG_CHG,
                                    HM_lyrG1_v_RS7G4_DMG_CHH)

### Hot conditions

## halleri side

LL_halG1_v_RS7G1_DMG <- consolidate(LL_halG1_v_RS7G1_DMG_CG,
                                    LL_halG1_v_RS7G1_DMG_CHG,
                                    LL_halG1_v_RS7G1_DMG_CHH)

LL_halG1_v_RS7G4_DMG <- consolidate(LL_halG1_v_RS7G4_DMG_CG,
                                    LL_halG1_v_RS7G4_DMG_CHG,
                                    LL_halG1_v_RS7G4_DMG_CHH)

## lyrata side

LL_lyrG1_v_RS7G1_DMG <- consolidate(LL_lyrG1_v_RS7G1_DMG_CG,
                                    LL_lyrG1_v_RS7G1_DMG_CHG,
                                    LL_lyrG1_v_RS7G1_DMG_CHH)

LL_lyrG1_v_RS7G4_DMG <- consolidate(LL_lyrG1_v_RS7G4_DMG_CG,
                                    LL_lyrG1_v_RS7G4_DMG_CHG,
                                    LL_lyrG1_v_RS7G4_DMG_CHH)


## Checking overlap across generations for DEGs and DMGs

## First look at diploid vs synthetic overlap

HM_pro1Vsyn1_halDEG <- sum(!is.na(match(HM_halG1_v_RS7G1_DEG$geneID, HM_halG1_v_RS7G4_DEG$geneID)))
paste0(round(HM_pro1Vsyn1_halDEG/length(HM_halG1_v_RS7G1_DEG$geneID)*100, 2), "%")

HM_pro1Vsyn1_halDMG <- sum(!is.na(match(HM_halG1_v_RS7G1_DMG$geneID, HM_halG1_v_RS7G4_DMG$geneID)))
paste0(round(HM_pro1Vsyn1_halDMG/length(HM_halG1_v_RS7G1_DMG$geneID)*100, 2), "%")

HM_pro1Vsyn1_lyrDEG <- sum(!is.na(match(HM_lyrG1_v_RS7G1_DEG$geneID, HM_lyrG1_v_RS7G4_DEG$geneID)))
paste0(round(HM_pro1Vsyn1_lyrDEG/length(HM_lyrG1_v_RS7G1_DEG$geneID)*100, 2), "%")

HM_pro1Vsyn1_lyrDMG <- sum(!is.na(match(HM_lyrG1_v_RS7G1_DMG$geneID, HM_lyrG1_v_RS7G4_DMG$geneID)))
paste0(round(HM_pro1Vsyn1_lyrDMG/length(HM_lyrG1_v_RS7G1_DMG$geneID)*100, 2), "%")

LL_pro1Vsyn1_halDEG <- sum(!is.na(match(LL_halG1_v_RS7G1_DEG$geneID, LL_halG1_v_RS7G4_DEG$geneID)))
paste0(round(LL_pro1Vsyn1_halDEG/length(LL_halG1_v_RS7G1_DEG$geneID)*100, 2), "%")

LL_pro1Vsyn1_halDMG <- sum(!is.na(match(LL_halG1_v_RS7G1_DMG$geneID, LL_halG1_v_RS7G4_DMG$geneID)))
paste0(round(LL_pro1Vsyn1_halDMG/length(LL_halG1_v_RS7G1_DMG$geneID)*100, 2), "%")

LL_pro1Vsyn1_lyrDEG <- sum(!is.na(match(LL_lyrG1_v_RS7G1_DEG$geneID, LL_lyrG1_v_RS7G4_DEG$geneID)))
paste0(round(LL_pro1Vsyn1_lyrDEG/length(LL_lyrG1_v_RS7G1_DEG$geneID)*100, 2), "%")

LL_pro1Vsyn1_lyrDMG <- sum(!is.na(match(LL_lyrG1_v_RS7G1_DMG$geneID, LL_lyrG1_v_RS7G4_DMG$geneID)))
paste0(round(LL_pro1Vsyn1_lyrDMG/length(LL_lyrG1_v_RS7G1_DMG$geneID)*100, 2), "%")


### We can also do a simple contingency
### table and test if the DMGs and DEGs are independent (H0)
### or dependent (H1)

# create contingency table for HM hal side

## G1
HM_table_halG1_v_RS7G1 <- matrix(c(length(intersect(HM_halG1_v_RS7G1_DEG$geneID,
                                                  HM_halG1_v_RS7G1_DMG$geneID)),
                                 nrow(HM_halG1_v_RS7G1_DEG) - length(intersect(HM_halG1_v_RS7G1_DEG$geneID,
                                                                                 HM_halG1_v_RS7G1_DMG$geneID)),
                                 nrow(HM_halG1_v_RS7G1_DMG) - length(intersect(HM_halG1_v_RS7G1_DEG$geneID,
                                                                               HM_halG1_v_RS7G1_DMG$geneID)),
                                 21518 - nrow(HM_halG1_v_RS7G1_DEG) - nrow(HM_halG1_v_RS7G1_DMG) + length(intersect(HM_halG1_v_RS7G1_DEG$geneID,
                                                                                                                    HM_halG1_v_RS7G1_DMG$geneID))),
                                 ncol = 2)
colnames(HM_table_halG1_v_RS7G1) <- c("DEG", "NotDEG")
rownames(HM_table_halG1_v_RS7G1) <- c("DMG", "NotDMG")

chisq.test(HM_table_halG1_v_RS7G1, correct = F)

## G4
HM_table_halG1_v_RS7G4 <- matrix(c(length(intersect(HM_halG1_v_RS7G4_DEG$geneID,
                                                    HM_halG1_v_RS7G4_DMG$geneID)),
                                   nrow(HM_halG1_v_RS7G4_DEG) - length(intersect(HM_halG1_v_RS7G4_DEG$geneID,
                                                                                 HM_halG1_v_RS7G4_DMG$geneID)),
                                   nrow(HM_halG1_v_RS7G4_DMG) - length(intersect(HM_halG1_v_RS7G4_DEG$geneID,
                                                                                 HM_halG1_v_RS7G4_DMG$geneID)),
                                   21518 - nrow(HM_halG1_v_RS7G4_DEG) - nrow(HM_halG1_v_RS7G4_DMG) + length(intersect(HM_halG1_v_RS7G4_DEG$geneID,
                                                                                                                      HM_halG1_v_RS7G4_DMG$geneID))),
                                 ncol = 2)
colnames(HM_table_halG1_v_RS7G4) <- c("DEG", "NotDEG")
rownames(HM_table_halG1_v_RS7G4) <- c("DMG", "NotDMG")

chisq.test(HM_table_halG1_v_RS7G4, correct = F)

# create contingency table for HM lyr side

## G1
HM_table_lyrG1_v_RS7G1 <- matrix(c(length(intersect(HM_lyrG1_v_RS7G1_DEG$geneID,
                                                    HM_lyrG1_v_RS7G1_DMG$geneID)),
                                   nrow(HM_lyrG1_v_RS7G1_DEG) - length(intersect(HM_lyrG1_v_RS7G1_DEG$geneID,
                                                                                 HM_lyrG1_v_RS7G1_DMG$geneID)),
                                   nrow(HM_lyrG1_v_RS7G1_DMG) - length(intersect(HM_lyrG1_v_RS7G1_DEG$geneID,
                                                                                 HM_lyrG1_v_RS7G1_DMG$geneID)),
                                   19595 - nrow(HM_lyrG1_v_RS7G1_DEG) - nrow(HM_lyrG1_v_RS7G1_DMG) + length(intersect(HM_lyrG1_v_RS7G1_DEG$geneID,
                                                                                                                      HM_lyrG1_v_RS7G1_DMG$geneID))),
                                 ncol = 2)
colnames(HM_table_lyrG1_v_RS7G1) <- c("DEG", "NotDEG")
rownames(HM_table_lyrG1_v_RS7G1) <- c("DMG", "NotDMG")

chisq.test(HM_table_lyrG1_v_RS7G1)

## G4
HM_table_lyrG1_v_RS7G4 <- matrix(c(length(intersect(HM_lyrG1_v_RS7G4_DEG$geneID,
                                                    HM_lyrG1_v_RS7G4_DMG$geneID)),
                                   nrow(HM_lyrG1_v_RS7G4_DEG) - length(intersect(HM_lyrG1_v_RS7G4_DEG$geneID,
                                                                                 HM_lyrG1_v_RS7G4_DMG$geneID)),
                                   nrow(HM_lyrG1_v_RS7G4_DMG) - length(intersect(HM_lyrG1_v_RS7G4_DEG$geneID,
                                                                                 HM_lyrG1_v_RS7G4_DMG$geneID)),
                                   19595 - nrow(HM_lyrG1_v_RS7G4_DEG) - nrow(HM_lyrG1_v_RS7G4_DMG) + length(intersect(HM_lyrG1_v_RS7G4_DEG$geneID,
                                                                                                                      HM_lyrG1_v_RS7G4_DMG$geneID))),
                                 ncol = 2)
colnames(HM_table_lyrG1_v_RS7G4) <- c("DEG", "NotDEG")
rownames(HM_table_lyrG1_v_RS7G4) <- c("DMG", "NotDMG")

chisq.test(HM_table_lyrG1_v_RS7G4)

# create contingency table for LL hal side

## G1
LL_table_halG1_v_RS7G1 <- matrix(c(length(intersect(LL_halG1_v_RS7G1_DEG$geneID,
                                                    LL_halG1_v_RS7G1_DMG$geneID)),
                                   nrow(LL_halG1_v_RS7G1_DEG) - length(intersect(LL_halG1_v_RS7G1_DEG$geneID,
                                                                                 LL_halG1_v_RS7G1_DMG$geneID)),
                                   nrow(LL_halG1_v_RS7G1_DMG) - length(intersect(LL_halG1_v_RS7G1_DEG$geneID,
                                                                                 LL_halG1_v_RS7G1_DMG$geneID)),
                                   21518 - nrow(LL_halG1_v_RS7G1_DEG) - nrow(LL_halG1_v_RS7G1_DMG) + length(intersect(LL_halG1_v_RS7G1_DEG$geneID,
                                                                                                                      LL_halG1_v_RS7G1_DMG$geneID))),
                                 ncol = 2)
colnames(LL_table_halG1_v_RS7G1) <- c("DEG", "NotDEG")
rownames(LL_table_halG1_v_RS7G1) <- c("DMG", "NotDMG")

chisq.test(LL_table_halG1_v_RS7G1, correct = F)

## G4
LL_table_halG1_v_RS7G4 <- matrix(c(length(intersect(LL_halG1_v_RS7G4_DEG$geneID,
                                                    LL_halG1_v_RS7G4_DMG$geneID)),
                                   nrow(LL_halG1_v_RS7G4_DEG) - length(intersect(LL_halG1_v_RS7G4_DEG$geneID,
                                                                                 LL_halG1_v_RS7G4_DMG$geneID)),
                                   nrow(LL_halG1_v_RS7G4_DMG) - length(intersect(LL_halG1_v_RS7G4_DEG$geneID,
                                                                                 LL_halG1_v_RS7G4_DMG$geneID)),
                                   21518 - nrow(LL_halG1_v_RS7G4_DEG) - nrow(LL_halG1_v_RS7G4_DMG) + length(intersect(LL_halG1_v_RS7G4_DEG$geneID,
                                                                                                                      LL_halG1_v_RS7G4_DMG$geneID))),
                                 ncol = 2)
colnames(LL_table_halG1_v_RS7G4) <- c("DEG", "NotDEG")
rownames(LL_table_halG1_v_RS7G4) <- c("DMG", "NotDMG")

chisq.test(LL_table_halG1_v_RS7G4, correct = F)

# create contingency table for LL lyr side

## G1
LL_table_lyrG1_v_RS7G1 <- matrix(c(length(intersect(LL_lyrG1_v_RS7G1_DEG$geneID,
                                                    LL_lyrG1_v_RS7G1_DMG$geneID)),
                                   nrow(LL_lyrG1_v_RS7G1_DEG) - length(intersect(LL_lyrG1_v_RS7G1_DEG$geneID,
                                                                                 LL_lyrG1_v_RS7G1_DMG$geneID)),
                                   nrow(LL_lyrG1_v_RS7G1_DMG) - length(intersect(LL_lyrG1_v_RS7G1_DEG$geneID,
                                                                                 LL_lyrG1_v_RS7G1_DMG$geneID)),
                                   19595 - nrow(LL_lyrG1_v_RS7G1_DEG) - nrow(LL_lyrG1_v_RS7G1_DMG) + length(intersect(LL_lyrG1_v_RS7G1_DEG$geneID,
                                                                                                                      LL_lyrG1_v_RS7G1_DMG$geneID))),
                                 ncol = 2)
colnames(LL_table_lyrG1_v_RS7G1) <- c("DEG", "NotDEG")
rownames(LL_table_lyrG1_v_RS7G1) <- c("DMG", "NotDMG")

chisq.test(LL_table_lyrG1_v_RS7G1)

## G4
LL_table_lyrG1_v_RS7G4 <- matrix(c(length(intersect(LL_lyrG1_v_RS7G4_DEG$geneID,
                                                    LL_lyrG1_v_RS7G4_DMG$geneID)),
                                   nrow(LL_lyrG1_v_RS7G4_DEG) - length(intersect(LL_lyrG1_v_RS7G4_DEG$geneID,
                                                                                 LL_lyrG1_v_RS7G4_DMG$geneID)),
                                   nrow(LL_lyrG1_v_RS7G4_DMG) - length(intersect(LL_lyrG1_v_RS7G4_DEG$geneID,
                                                                                 LL_lyrG1_v_RS7G4_DMG$geneID)),
                                   19595 - nrow(LL_lyrG1_v_RS7G4_DEG) - nrow(LL_lyrG1_v_RS7G4_DMG) + length(intersect(LL_lyrG1_v_RS7G4_DEG$geneID,
                                                                                                                      LL_lyrG1_v_RS7G4_DMG$geneID))),
                                 ncol = 2)
colnames(LL_table_lyrG1_v_RS7G4) <- c("DEG", "NotDEG")
rownames(LL_table_lyrG1_v_RS7G4) <- c("DMG", "NotDMG")

chisq.test(LL_table_lyrG1_v_RS7G4)
