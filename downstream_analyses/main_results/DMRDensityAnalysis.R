# This script will compute the density of DMRs
# in homeologous exchanged regions and normal regions

### import libraries

library(ggplot2)
library(GenomicRanges)
library(tidyverse)

### Import files

homeo_HM_overlap_continuous <- read.csv("~/OneDrive/PhD/Project/Chapter_3/coverage_analysis/homeo_regions/homeo_HM_overlap_continuous.txt", sep="")

homeo_LL_overlap_continuous <- read.csv("~/OneDrive/PhD/Project/Chapter_3/coverage_analysis/homeo_regions/homeo_LL_overlap_continuous.txt", sep="")

nonhomeo_HM_overlap_continuous <- read.csv("~/OneDrive/PhD/Project/Chapter_3/coverage_analysis/homeo_regions/nonhomeo_HM_overlap_continuous.txt", sep="")

nonhomeo_LL_overlap_continuous <- read.csv("~/OneDrive/PhD/Project/Chapter_3/coverage_analysis/homeo_regions/nonhomeo_LL_overlap_continuous.txt", sep="")

folder_path <- c("~/OneDrive/PhD/Project/Chapter_3/DMR_results/ARPEGGIO_results_syn4Vpro1_HM/dmrseq/")

HM_RS7_G4_CG <- read.csv(paste0(folder_path, "CG_context/parent1_v_allo.txt"))
HM_RS7_G4_CHG <- read.csv(paste0(folder_path, "CHG_context/parent1_v_allo.txt"))
HM_RS7_G4_CHH <- read.csv(paste0(folder_path, "CHH_context/parent1_v_allo.txt"))

folder_path <- c("~/OneDrive/PhD/Project/Chapter_3/DMR_results/ARPEGGIO_results_syn4Vpro1_LL/dmrseq/")

LL_RS7_G4_CG <- read.csv(paste0(folder_path, "CG_context/parent1_v_allo.txt"))
LL_RS7_G4_CHG <- read.csv(paste0(folder_path, "CHG_context/parent1_v_allo.txt"))
LL_RS7_G4_CHH <- read.csv(paste0(folder_path, "CHH_context/parent1_v_allo.txt"))

### Coordinate converter function to be able to select DMRs

coordinate_converter <- function(dmrseq_output, agp){

      ### Match old scaffolds to new chromosomes

      scaffold_match <- match(dmrseq_output$seqnames, agp$old_scaffold)

      ### Create vector with new scaffolds and new coordinates

      new_scaffolds <- agp$new_scaffold[scaffold_match]
      new_start <- agp$start[scaffold_match] + dmrseq_output$start
      new_end <- agp$start[scaffold_match] + dmrseq_output$end

      ### Modify dmrseq output with new coordinates

      new_dmrseq_output <- dmrseq_output
      new_dmrseq_output$seqnames <- new_scaffolds
      new_dmrseq_output$start <- new_start
      new_dmrseq_output$end <- new_end

      # Selecting only 8 chromosomes

      new_dmrseq_output <- new_dmrseq_output %>%
        dplyr::filter(seqnames=="scaffold_1_RagTag" |
                        seqnames=="scaffold_2_RagTag" |
                        seqnames=="scaffold_3_RagTag" |
                        seqnames=="scaffold_4_RagTag" |
                        seqnames=="scaffold_5_RagTag" |
                        seqnames=="scaffold_6_RagTag" |
                        seqnames=="scaffold_7_RagTag" |
                        seqnames=="scaffold_8_RagTag")

      ### return file with new coordinates

      return(new_dmrseq_output)
}

agp <- read.delim("~/OneDrive/PhD/Project/Chapter2/AssemblyMapping/Ahal_mm2_tiling.agp", 
                  header=FALSE, comment.char="#")

names(agp) <- c("new_scaffold", 
                "start", 
                "end", 
                "n", 
                "WU", 
                "old_scaffold", 
                "type", 
                "length", 
                "orientiation")

new_HM_RS7_G4_CG <- coordinate_converter(HM_RS7_G4_CG, agp)
new_HM_RS7_G4_CHG <- coordinate_converter(HM_RS7_G4_CHG, agp)
new_HM_RS7_G4_CHH <- coordinate_converter(HM_RS7_G4_CHH, agp)
new_LL_RS7_G4_CG <- coordinate_converter(LL_RS7_G4_CG, agp)
new_LL_RS7_G4_CHG <- coordinate_converter(LL_RS7_G4_CHG, agp)
new_LL_RS7_G4_CHH <- coordinate_converter(LL_RS7_G4_CHH, agp)

## merge all contexts and select significant regions only

all_HM_RS7_G4 <- rbind(new_HM_RS7_G4_CG,
                       new_HM_RS7_G4_CHG,
                       new_HM_RS7_G4_CHH)

all_LL_RS7_G4 <- rbind(new_LL_RS7_G4_CG,
                       new_LL_RS7_G4_CHG,
                       new_LL_RS7_G4_CHH)

all_HM_RS7_G4_sig <- filter(all_HM_RS7_G4, qval < 0.05)
all_LL_RS7_G4_sig <- filter(all_LL_RS7_G4, qval < 0.05)

## convert regions to Granges objects to find overlap with
## homeo or non homeo regions

all_HM_RS7_G4_sig_G <- GRanges(seqnames = all_HM_RS7_G4_sig$seqnames,
                               ranges = IRanges(start = all_HM_RS7_G4_sig$start,
                                                end = all_HM_RS7_G4_sig$end))

all_LL_RS7_G4_sig_G <- GRanges(seqnames = all_LL_RS7_G4_sig$seqnames,
                               ranges = IRanges(start = all_LL_RS7_G4_sig$start,
                                                end = all_LL_RS7_G4_sig$end))

homeo_HM_overlap_continuous_G <- GRanges(seqnames = homeo_HM_overlap_continuous$seqnames,
                                             ranges = IRanges(start = homeo_HM_overlap_continuous$start,
                                                              end = homeo_HM_overlap_continuous$end))


homeo_LL_overlap_continuous_G <- GRanges(seqnames = homeo_LL_overlap_continuous$seqnames,
                                             ranges = IRanges(start = homeo_LL_overlap_continuous$start,
                                                              end = homeo_LL_overlap_continuous$end))

nonhomeo_HM_overlap_continuous_G <- GRanges(seqnames = nonhomeo_HM_overlap_continuous$start_scaffold,
                                             ranges = IRanges(start = nonhomeo_HM_overlap_continuous$start_coordinate,
                                                              end = nonhomeo_HM_overlap_continuous$end_coordinate))


nonhomeo_LL_overlap_continuous_G <- GRanges(seqnames = nonhomeo_LL_overlap_continuous$start_scaffold,
                                             ranges = IRanges(start = nonhomeo_LL_overlap_continuous$start_coordinate,
                                                              end = nonhomeo_LL_overlap_continuous$end_coordinate))


# Check overlaps

DMRs_homeo_HM <- findOverlaps(all_HM_RS7_G4_sig_G, homeo_HM_overlap_continuous_G)

DMRhomeo_HM <- all_HM_RS7_G4_sig[DMRs_homeo_HM@from,]

DMRs_nonhomeo_HM <- findOverlaps(all_HM_RS7_G4_sig_G, 
                                 nonhomeo_HM_overlap_continuous_G)

DMRnonhomeo_HM <- all_HM_RS7_G4_sig[DMRs_nonhomeo_HM@from,]


DMRs_homeo_LL <- findOverlaps(all_LL_RS7_G4_sig_G, homeo_LL_overlap_continuous_G)

DMRhomeo_LL <- all_LL_RS7_G4_sig[DMRs_homeo_LL@from,]

DMRs_nonhomeo_LL <- findOverlaps(all_LL_RS7_G4_sig_G, 
                                 nonhomeo_LL_overlap_continuous_G)

DMRnonhomeo_LL <- all_LL_RS7_G4_sig[DMRs_nonhomeo_LL@from,]

# Let's first look at DMR density over the whole range

homeo_range_HM <- sum(homeo_HM_overlap_continuous$end -
                        homeo_HM_overlap_continuous$start)

homeo_range_LL <- sum(homeo_LL_overlap_continuous$end -
                        homeo_LL_overlap_continuous$start)

nonhomeo_range_HM <- sum(nonhomeo_HM_overlap_continuous$end_coordinate -
                        nonhomeo_HM_overlap_continuous$start_coordinate)

nonhomeo_range_LL <- sum(nonhomeo_LL_overlap_continuous$end_coordinate -
                        nonhomeo_LL_overlap_continuous$start_coordinate)

global_density_homeo_HM <- nrow(DMRhomeo_HM) / homeo_range_HM

global_density_homeo_LL <- nrow(DMRhomeo_LL) / homeo_range_LL

global_density_nonhomeo_HM <- nrow(DMRnonhomeo_HM) / nonhomeo_range_HM

global_density_nonhomeo_LL <- nrow(DMRnonhomeo_LL) / nonhomeo_range_LL

global_density_homeo_HM / global_density_nonhomeo_HM

global_density_homeo_LL / global_density_nonhomeo_LL

## It looks like in HM the density of DMRs is 1.5 times higher in exchanged
## regions than non-exchanged, and this is even higher in hot conditions 

## We take a look at the proportion of hypo and hypermethylated regions to
## see if the proportion shows any interesting pattern

# Hypermethylated DMRs in HE in cold conditions
sum(DMRhomeo_HM$stat < 0)

# Hypomethylated DMRs in HE in cold conditions
sum(DMRhomeo_HM$stat > 0)

# Hypermethylated DMRs in normal regions in cold conditions
sum(DMRnonhomeo_HM$stat < 0)

# Hypomethylated DMRs in normal regions in cold conditions
sum(DMRnonhomeo_HM$stat > 0)

# Hypermethylated DMRs in HE in hot conditions
sum(DMRhomeo_LL$stat < 0)

# Hypomethylated DMRs in HE in hot conditions
sum(DMRhomeo_LL$stat > 0)

# Hypermethylated DMRs in normal regions in hot conditions
sum(DMRnonhomeo_LL$stat < 0)

# Hypomethylated DMRs in normal regions in hot conditions
sum(DMRnonhomeo_LL$stat > 0)


## The exchanged regions have lower density of decreased methylation
## compared to the non exchanged regions
