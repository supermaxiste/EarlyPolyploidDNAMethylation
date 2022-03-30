### RNA spatial distribution

### Import libraries

library(tidyverse)
library(karyoploteR)

## Path to density data folder

folder_path <- c("coverage_data/")

## Import all densities

HM_hal_G1_1_dens <- read.csv(paste0(folder_path, "HM_hal_G1_1_dens.txt"))

HM_lyr_G1_1_dens <- read.csv(paste0(folder_path, "HM_lyr_G1_1_dens.txt"))

## polyploid data

HM_RS7K_G1_1_hal_dens <- read.csv(paste0(folder_path, "HM_RS7_G1_1_hal_dens.txt"))
HM_RS7K_G1_2_hal_dens <- read.csv(paste0(folder_path, "HM_RS7_G1_2_hal_dens.txt"))
HM_RS7K_G1_3_hal_dens <- read.csv(paste0(folder_path, "HM_RS7_G1_3_hal_dens.txt"))

HM_RS7K_G1_1_lyr_dens <- read.csv(paste0(folder_path, "HM_RS7_G1_1_lyr_dens.txt"))
HM_RS7K_G1_2_lyr_dens <- read.csv(paste0(folder_path, "HM_RS7_G1_2_lyr_dens.txt"))
HM_RS7K_G1_3_lyr_dens <- read.csv(paste0(folder_path, "HM_RS7_G1_3_lyr_dens.txt"))

HM_RS7_G4_1_hal_dens <- read.csv(paste0(folder_path, "HM_RS7_G4_1_hal_dens.txt"))
HM_RS7_G4_2_hal_dens <- read.csv(paste0(folder_path, "HM_RS7_G4_2_hal_dens.txt"))
HM_RS7_G4_3_hal_dens <- read.csv(paste0(folder_path, "HM_RS7_G4_3_hal_dens.txt"))

HM_RS7_G4_1_lyr_dens <- read.csv(paste0(folder_path, "HM_RS7_G4_1_lyr_dens.txt"))
HM_RS7_G4_2_lyr_dens <- read.csv(paste0(folder_path, "HM_RS7_G4_2_lyr_dens.txt"))
HM_RS7_G4_3_lyr_dens <- read.csv(paste0(folder_path, "HM_RS7_G4_3_lyr_dens.txt"))

HM_RS7_G4_1_hal_dens_old <- read.csv(paste0(folder_path, "old/RS7_G4_1_hal_dens.txt"))
HM_RS7_G4_2_hal_dens_old <- read.csv(paste0(folder_path, "old/RS7_G4_2_hal_dens.txt"))
HM_RS7_G4_3_hal_dens_old <- read.csv(paste0(folder_path, "old/RS7_G4_3_hal_dens.txt"))

HM_RS7_G4_1_lyr_dens_old <- read.csv(paste0(folder_path, "old/RS7_G4_1_lyr_dens.txt"))
HM_RS7_G4_2_lyr_dens_old <- read.csv(paste0(folder_path, "old/RS7_G4_2_lyr_dens.txt"))
HM_RS7_G4_3_lyr_dens_old <- read.csv(paste0(folder_path, "old/RS7_G4_3_lyr_dens.txt"))

LL_RS7_G1_1_hal_dens <- read.csv(paste0(folder_path, "LL_RS7_G1_1_hal_dens.txt"))
LL_RS7_G1_2_hal_dens <- read.csv(paste0(folder_path, "LL_RS7_G1_2_hal_dens.txt"))
LL_RS7_G1_3_hal_dens <- read.csv(paste0(folder_path, "LL_RS7_G1_3_hal_dens.txt"))

LL_RS7_G1_1_lyr_dens <- read.csv(paste0(folder_path, "LL_RS7_G1_1_lyr_dens.txt"))
LL_RS7_G1_2_lyr_dens <- read.csv(paste0(folder_path, "LL_RS7_G1_2_lyr_dens.txt"))
LL_RS7_G1_3_lyr_dens <- read.csv(paste0(folder_path, "LL_RS7_G1_3_lyr_dens.txt"))

LL_RS7K_G4_1_hal_dens <- read.csv(paste0(folder_path, "LL_RS7K_G4_1_hal_dens.txt"))
LL_RS7K_G4_2_hal_dens <- read.csv(paste0(folder_path, "LL_RS7K_G4_2_hal_dens.txt"))
LL_RS7K_G4_3_hal_dens <- read.csv(paste0(folder_path, "LL_RS7K_G4_3_hal_dens.txt"))

LL_RS7K_G4_1_lyr_dens <- read.csv(paste0(folder_path, "LL_RS7K_G4_1_lyr_dens.txt"))
LL_RS7K_G4_2_lyr_dens <- read.csv(paste0(folder_path, "LL_RS7K_G4_2_lyr_dens.txt"))
LL_RS7K_G4_3_lyr_dens <- read.csv(paste0(folder_path, "LL_RS7K_G4_3_lyr_dens.txt"))

## total density

hal_totC_dens <- read.delim("coverage_data/hal_totC_dens.txt")
lyr_totC_dens <- read.delim("coverage_data/lyr_totC_dens.txt")

### We first compute the real cytosine coverage per window for all samples

## polyploid data

HM_RS7K_G1_1_hal_real_coverage <- HM_RS7K_G1_1_hal_dens$window_scores / hal_totC_dens$dens
HM_RS7K_G1_2_hal_real_coverage <- HM_RS7K_G1_2_hal_dens$window_scores / hal_totC_dens$dens
HM_RS7K_G1_3_hal_real_coverage <- HM_RS7K_G1_3_hal_dens$window_scores / hal_totC_dens$dens

HM_RS7K_G1_1_lyr_real_coverage <- HM_RS7K_G1_1_lyr_dens$window_scores / lyr_totC_dens$dens
HM_RS7K_G1_2_lyr_real_coverage <- HM_RS7K_G1_2_lyr_dens$window_scores / lyr_totC_dens$dens
HM_RS7K_G1_3_lyr_real_coverage <- HM_RS7K_G1_3_lyr_dens$window_scores / lyr_totC_dens$dens

HM_RS7_G4_1_hal_real_coverage <- HM_RS7_G4_1_hal_dens$window_scores / hal_totC_dens$dens
HM_RS7_G4_2_hal_real_coverage <- HM_RS7_G4_2_hal_dens$window_scores / hal_totC_dens$dens
HM_RS7_G4_3_hal_real_coverage <- HM_RS7_G4_3_hal_dens$window_scores / hal_totC_dens$dens

HM_RS7_G4_1_lyr_real_coverage <- HM_RS7_G4_1_lyr_dens$window_scores / lyr_totC_dens$dens
HM_RS7_G4_2_lyr_real_coverage <- HM_RS7_G4_2_lyr_dens$window_scores / lyr_totC_dens$dens
HM_RS7_G4_3_lyr_real_coverage <- HM_RS7_G4_3_lyr_dens$window_scores / lyr_totC_dens$dens

HM_RS7_G4_1_hal_real_coverage_old <- HM_RS7_G4_1_hal_dens_old$window_scores / hal_totC_dens$dens
HM_RS7_G4_2_hal_real_coverage_old <- HM_RS7_G4_2_hal_dens_old$window_scores / hal_totC_dens$dens
HM_RS7_G4_3_hal_real_coverage_old <- HM_RS7_G4_3_hal_dens_old$window_scores / hal_totC_dens$dens

HM_RS7_G4_1_lyr_real_coverage_old <- HM_RS7_G4_1_lyr_dens_old$window_scores / lyr_totC_dens$dens
HM_RS7_G4_2_lyr_real_coverage_old <- HM_RS7_G4_2_lyr_dens_old$window_scores / lyr_totC_dens$dens
HM_RS7_G4_3_lyr_real_coverage_old <- HM_RS7_G4_3_lyr_dens_old$window_scores / lyr_totC_dens$dens

LL_RS7_G1_1_hal_real_coverage <- LL_RS7_G1_1_hal_dens$window_scores / hal_totC_dens$dens
LL_RS7_G1_2_hal_real_coverage <- LL_RS7_G1_2_hal_dens$window_scores / hal_totC_dens$dens
LL_RS7_G1_3_hal_real_coverage <- LL_RS7_G1_3_hal_dens$window_scores / hal_totC_dens$dens

LL_RS7_G1_1_lyr_real_coverage <- LL_RS7_G1_1_lyr_dens$window_scores / lyr_totC_dens$dens
LL_RS7_G1_2_lyr_real_coverage <- LL_RS7_G1_2_lyr_dens$window_scores / lyr_totC_dens$dens
LL_RS7_G1_3_lyr_real_coverage <- LL_RS7_G1_3_lyr_dens$window_scores / lyr_totC_dens$dens

LL_RS7K_G4_1_hal_real_coverage <- LL_RS7K_G4_1_hal_dens$window_scores / hal_totC_dens$dens
LL_RS7K_G4_2_hal_real_coverage <- LL_RS7K_G4_2_hal_dens$window_scores / hal_totC_dens$dens
LL_RS7K_G4_3_hal_real_coverage <- LL_RS7K_G4_3_hal_dens$window_scores / hal_totC_dens$dens

LL_RS7K_G4_1_lyr_real_coverage <- LL_RS7K_G4_1_lyr_dens$window_scores / lyr_totC_dens$dens
LL_RS7K_G4_2_lyr_real_coverage <- LL_RS7K_G4_2_lyr_dens$window_scores / lyr_totC_dens$dens
LL_RS7K_G4_3_lyr_real_coverage <- LL_RS7K_G4_3_lyr_dens$window_scores / lyr_totC_dens$dens

# Compute average for each sample

HM_RS7K_G1_hal_avg_coverage <- c((HM_RS7K_G1_1_hal_real_coverage +
                                    HM_RS7K_G1_2_hal_real_coverage +
                                    HM_RS7K_G1_3_hal_real_coverage) / 3)

HM_RS7K_G1_lyr_avg_coverage <- c((HM_RS7K_G1_1_lyr_real_coverage +
                                    HM_RS7K_G1_2_lyr_real_coverage +
                                    HM_RS7K_G1_3_lyr_real_coverage) / 3)

HM_RS7_G4_hal_avg_coverage <- c((HM_RS7_G4_1_hal_real_coverage +
                                   HM_RS7_G4_2_hal_real_coverage +
                                   HM_RS7_G4_3_hal_real_coverage) / 3)

HM_RS7_G4_lyr_avg_coverage <- c((HM_RS7_G4_1_lyr_real_coverage +
                                   HM_RS7_G4_2_lyr_real_coverage +
                                   HM_RS7_G4_3_lyr_real_coverage) / 3)

HM_RS7_G4_hal_avg_coverage_old <- c((HM_RS7_G4_1_hal_real_coverage_old +
                                       HM_RS7_G4_2_hal_real_coverage_old +
                                       HM_RS7_G4_3_hal_real_coverage_old) / 3)

HM_RS7_G4_lyr_avg_coverage_old <- c((HM_RS7_G4_1_lyr_real_coverage_old +
                                       HM_RS7_G4_2_lyr_real_coverage_old +
                                       HM_RS7_G4_3_lyr_real_coverage_old) / 3)

LL_RS7_G1_hal_avg_coverage <- c((LL_RS7_G1_1_hal_real_coverage +
                                   LL_RS7_G1_2_hal_real_coverage +
                                   LL_RS7_G1_3_hal_real_coverage) / 3)

LL_RS7_G1_lyr_avg_coverage <- c((LL_RS7_G1_1_lyr_real_coverage +
                                   LL_RS7_G1_2_lyr_real_coverage +
                                   LL_RS7_G1_3_lyr_real_coverage) / 3)

LL_RS7K_G4_hal_avg_coverage <- c((LL_RS7K_G4_1_hal_real_coverage +
                                    LL_RS7K_G4_2_hal_real_coverage +
                                    LL_RS7K_G4_3_hal_real_coverage) / 3)

LL_RS7K_G4_lyr_avg_coverage <- c((LL_RS7K_G4_1_lyr_real_coverage +
                                    LL_RS7K_G4_2_lyr_real_coverage +
                                    LL_RS7K_G4_3_lyr_real_coverage) / 3)

### Compute total windows per chromosome

chrom_length_hal <- c(1:8)
chrom_length_lyr <- c(1:8)
counter <- 1
for (i in c("scaffold_1_RagTag", "scaffold_2_RagTag", "scaffold_3_RagTag",
            "scaffold_4_RagTag", "scaffold_5_RagTag", "scaffold_6_RagTag",
            "scaffold_7_RagTag", "scaffold_8_RagTag")){
  chrom_length_hal[counter] <- sum(HM_hal_G1_1_dens$seqnames == i)
  chrom_length_lyr[counter] <- sum(HM_lyr_G1_1_dens$seqnames == i)
  counter <- counter + 1
}

average_chromosomes_hal <- c()
average_chromosomes_lyr <- c()

for (i in c(1:8)){
  average_chromosomes_hal <- append(average_chromosomes_hal, seq(from = 50000, 
                                 to = chrom_length_hal[i]*100000, 
                                 by = 100000))
  average_chromosomes_lyr <- append(average_chromosomes_lyr, seq(from = 50000, 
                                        to = chrom_length_lyr[i]*100000, 
                                        by = 100000))
}

## Create GRanges object for plotting coverage over chromosomes

#hal_coverage_G <- GRanges(seqnames = HM_RS7K_G1_1_hal_dens$seqnames,
                                    #ranges = IRanges(start = average_chromosomes_hal, end = average_chromosomes_hal))

### We define a function to plot DEG 
### on chromosomes with karyoploteR


DEGs_karyoploter <- function(deg_output, side, coverage, chrom_length, average_chromosomes){

  direction <- ifelse(deg_output$dir>0, "increase", "decrease")

  deg_output <- mutate(deg_output, direction = direction)
  
  # We import the annotations to add genomic coordinates
  # for each DEG based on their geneID
  
  annotation_import <- function(side){
    if(side=="hal"){
      anno <- read.delim("Ahal_v2_2.gff",
                         header=FALSE,
                         col.names = c("scaffold", 
                                       "tool", 
                                       "context", 
                                       "start", 
                                       "end", 
                                       "number",
                                       "strand", 
                                       "dot", 
                                       "extra"))
    }
    else {
      anno <- read.delim("Alyr_v2_2_renamed.gff",
                         header=FALSE,
                         col.names = c("scaffold", 
                                       "tool", 
                                       "context", 
                                       "start", 
                                       "end", 
                                       "number",
                                       "strand", 
                                       "dot", 
                                       "extra"))
    }
    return(anno)
  }
  
  anno <- annotation_import(side)
  
  anno <- filter(anno, context=="gene")
  
  genes_deg <- as.numeric(as.data.frame(str_split_fixed(deg_output$geneID, "g", n=Inf))$V2)
  
  deg_output <- mutate(deg_output, 
                       scaffold = anno$scaffold[genes_deg],
                       start = anno$start[genes_deg],
                       end = anno$end[genes_deg])
  
  ## We convert the new coordinates to chromosome ones

  coordinate_converter <- function(deg_output, agp){

    ### Match old scaffolds to new chromosomes

    scaffold_match <- match(deg_output$seqnames, agp$old_scaffold)

    ### Create vector with new scaffolds and new coordinates

    new_scaffolds <- agp$new_scaffold[scaffold_match]
    new_start <- agp$start[scaffold_match] + deg_output$start
    new_end <- agp$start[scaffold_match] + deg_output$end

    ### Modify deg output with new coordinates

    new_deg_output <- deg_output
    new_deg_output$seqnames <- new_scaffolds
    new_deg_output$start <- new_start
    new_deg_output$end <- new_end

    # Selecting only 8 chromosomes

    new_deg_output <- new_deg_output %>%
      dplyr::filter(seqnames=="scaffold_1_RagTag" |
                      seqnames=="scaffold_2_RagTag" |
                      seqnames=="scaffold_3_RagTag" |
                      seqnames=="scaffold_4_RagTag" |
                      seqnames=="scaffold_5_RagTag" |
                      seqnames=="scaffold_6_RagTag" |
                      seqnames=="scaffold_7_RagTag" |
                      seqnames=="scaffold_8_RagTag")

    ### return file with new coordinates

    return(new_deg_output)
  }

  # We define function to import agp file
  
  import_agp <- function(side){
    
    if (side=="hal"){
      
      ### Import agp file of the corresponding assembly
      
      agp <- read.delim("agp_files/Ahal_mum_tiling.agp", 
                        header=FALSE,
                        comment.char="#")
      
      names(agp) <- c("new_scaffold", 
                      "start", 
                      "end", 
                      "n", 
                      "WU", 
                      "old_scaffold", 
                      "type", 
                      "length", 
                      "orientiation")
  }
  else{

    agp <- read.delim("agp_files/Alyr_mum_tiling.agp", 
                          header=FALSE, 
                          comment.char="#")
    
    names(agp) <- c("new_scaffold", 
                        "start", 
                        "end", 
                        "n", 
                        "WU", 
                        "old_scaffold", 
                        "type", 
                        "length", 
                        "orientiation")
    
    agp_scaffold_renamer <- function(agp){
      
      # Get old scaffold numbers
      old_scaffolds <- agp$old_scaffold
      
      # Split the string to just get the number
      scaffold_number <- as.integer(str_split_fixed(old_scaffolds, "_", 2)[,2])
      
      # Add 2249 for the new numbering
      new_scaffold_number <- scaffold_number + 2239
      
      # Create new set of strings to put in the agp file
      new_scaffolds <- paste0("scaffold_", new_scaffold_number)
      
      # Replace old names with new names
      agp$old_scaffold <- new_scaffolds
      
      return(agp)
    }
    agp <- agp_scaffold_renamer(agp)
  }
    return(agp)
  }
  
  # We import the agp file
  
  agp <- import_agp(side)
  
  # We create function for converting coordinates
  # to chromosome coordinates
  
  chromosome_converter <- function(deg_output, agp){
    
    ### Match old scaffolds to new chromosomes
    
    scaffold_match <- match(deg_output$scaffold,
                            agp$old_scaffold)
    
    ### Create vector with new scaffolds 
    ### and new coordinates

    new_scaffolds <- agp$new_scaffold[scaffold_match]
    new_start <- agp$start[scaffold_match] + 
      deg_output$start
    new_end <- agp$start[scaffold_match] +
      deg_output$end

    ### Modify dmrseq output with new coordinates

    new_deg_output <- deg_output
    new_deg_output$scaffold <- new_scaffolds
    new_deg_output$start <- new_start
    new_deg_output$end <- new_end
    
    # Selecting only 8 chromosomes

    new_deg_output <- new_deg_output %>%
      dplyr::filter(scaffold=="scaffold_1_RagTag" |
                      scaffold=="scaffold_2_RagTag" |
                      scaffold=="scaffold_3_RagTag" |
                      scaffold=="scaffold_4_RagTag" |
                      scaffold=="scaffold_5_RagTag" |
                      scaffold=="scaffold_6_RagTag" |
                      scaffold=="scaffold_7_RagTag" |
                      scaffold=="scaffold_8_RagTag")

    ### return file with new coordinates
    
    return(new_deg_output)
    
  }
  
  new_deg_output <- chromosome_converter(deg_output,
                                         agp)
  
  ### Define centromere coordinates
  
  centromere_coord <-   c(18589977,
                          5071832,
                          15417662,
                          7248175,
                          7610403,
                          14782549,
                          15405871,
                          7623002)

  # Only black for DMRs

  new_deg_output_col <- data.frame(Chrom = as.character(new_deg_output$scaffold),
                                   Start = as.integer(new_deg_output$start),
                                      End = as.integer(new_deg_output$end),
                                      Name = as.character(new_deg_output$scaffold),
                                      gieStain = "black")


  # Define chromosome sizes
  
  import_chromSizes <- function(side){
    if(side=="hal"){
      chrom_sizes <- data.frame(Chrom = c(1:8),
                                Start = c(rep(1, 8)),
                                End = c(29853912, 
                                        19430572, 
                                        22405798,
                                        22332524, 
                                        17489651, 
                                        21293240,
                                        25127945, 
                                        21784696),
                                Name = rep(NA, 8),
                                Colors = rep("lightgrey", 
                                             8))
      
    }
    
    else{
      # Define chromosome sizes lyrata
      chrom_sizes <- data.frame(Chrom = c(1:8),
                                Start = c(rep(1, 8)),
                                End = c(26235489, 
                                        17728226, 
                                        22610115,
                                        18792104, 
                                        21074452, 
                                        18886653,
                                        22122096, 
                                        18984421),
                                Name = rep(NA, 8),
                                Colors = rep("lightgrey", 
                                             8))
    }
    return(chrom_sizes)
  }
  
  # Import chromosome sizes
  
  chrom_sizes <- import_chromSizes(side)

  # Define vector with chromosome names
  chrom_names <- chrom_sizes$End
  names(chrom_names) <- chrom_sizes$Chrom

  chrom_sizes_G <- GRanges(seqnames = chrom_sizes$Chrom,
                           ranges = IRanges(start = chrom_sizes$Start,
                                            end = chrom_sizes$End),
                           seqlengths = chrom_names)

  # Turn chromosome scaffolds into simple numbers
  
  new_deg_output$scaffold[new_deg_output$scaffold == "scaffold_1_RagTag"] <- 1
  new_deg_output$scaffold[new_deg_output$scaffold == "scaffold_2_RagTag"] <- 2
  new_deg_output$scaffold[new_deg_output$scaffold == "scaffold_3_RagTag"] <- 3
  new_deg_output$scaffold[new_deg_output$scaffold == "scaffold_4_RagTag"] <- 4
  new_deg_output$scaffold[new_deg_output$scaffold == "scaffold_5_RagTag"] <- 5
  new_deg_output$scaffold[new_deg_output$scaffold == "scaffold_6_RagTag"] <- 6
  new_deg_output$scaffold[new_deg_output$scaffold == "scaffold_7_RagTag"] <- 7
  new_deg_output$scaffold[new_deg_output$scaffold == "scaffold_8_RagTag"] <- 8
  
  # Turn deg output into GRanges object
  
  new_deg_output_G <- GRanges(seqnames = new_deg_output$scaffold,
                          ranges = IRanges(start = new_deg_output$start,
                                           end = new_deg_output$end))
  
  # Plot with karyoploter

  # First set parameters

  plot.params <- getDefaultPlotParams(plot.type=2)
  plot.params$data1height <- 50
  plot.params$data2height <- 50
  plot.params$ideogramheight <- 50

  kp <- plotKaryotype(genome = chrom_sizes_G, plot.type = 2, plot.params = plot.params)
  kpDataBackground(kp, data.panel = 1, color = "grey")
  kpDataBackground(kp, data.panel = 2, color = "white")
  kpDataBackground(kp, data.panel = "ideogram", color = "white")
  
  kpLines(kp, data.panel = 1, chr = rep(c(1:8), chrom_length), x = average_chromosomes, y = coverage, col = "black", ymax = 20)
  kpAxis(kp, data.panel = 1, side = 2, numticks = 2, tick.pos = c(0, 1), labels = c(0, 20), cex = 0.5)
  kpPlotRegions(kp, data = new_deg_output_G, col = c("black"), border = NA, data.panel = "ideogram")
  
}

HM_halG1_v_RS7G1_DEG <- read.delim("DEGs/HM_halG1_v_RS7G1_DEG.txt")
HM_halG1_v_RS7G4_DEG <- read.delim("DEGs/HM_halG1_v_RS7G4_DEG.txt")
HM_lyrG1_v_RS7G1_DEG <- read.delim("DEGs/HM_lyrG1_v_RS7G1_DEG.txt")
HM_lyrG1_v_RS7G4_DEG <- read.delim("DEGs/HM_lyrG1_v_RS7G4_DEG.txt")

DEGs_karyoploter(HM_halG1_v_RS7G1_DEG, 
                 "hal", 
                 HM_RS7K_G1_1_hal_real_coverage,
                 chrom_length_hal,
                 average_chromosomes_hal)
DEGs_karyoploter(HM_halG1_v_RS7G4_DEG, 
                 "hal",
                 HM_RS7_G4_1_hal_real_coverage,
                 chrom_length_hal,
                 average_chromosomes_hal)
DEGs_karyoploter(HM_lyrG1_v_RS7G1_DEG, 
                 "lyr",
                 HM_RS7K_G1_1_lyr_real_coverage,
                 chrom_length_lyr,
                 average_chromosomes_lyr)
DEGs_karyoploter(HM_lyrG1_v_RS7G4_DEG, 
                 "lyr",
                 HM_RS7_G4_1_lyr_real_coverage,
                 chrom_length_lyr,
                 average_chromosomes_lyr)

LL_halG1_v_RS7G1_DEG <- read.delim("DEGs/LL_halG1_v_RS7G1_DEG.txt")
LL_halG1_v_RS7G4_DEG <- read.delim("DEGs/LL_halG1_v_RS7G4_DEG.txt")
LL_lyrG1_v_RS7G1_DEG <- read.delim("DEGs/LL_lyrG1_v_RS7G1_DEG.txt")
LL_lyrG1_v_RS7G4_DEG <- read.delim("DEGs/LL_lyrG1_v_RS7G4_DEG.txt")

DEGs_karyoploter(LL_halG1_v_RS7G1_DEG, 
                 "hal", 
                 LL_RS7_G1_1_hal_real_coverage,
                 chrom_length_hal,
                 average_chromosomes_hal)
DEGs_karyoploter(LL_halG1_v_RS7G4_DEG, 
                 "hal", 
                 LL_RS7K_G4_1_hal_real_coverage,
                 chrom_length_hal,
                 average_chromosomes_hal)
DEGs_karyoploter(LL_lyrG1_v_RS7G1_DEG, 
                 "lyr",
                 LL_RS7_G1_1_lyr_real_coverage,
                 chrom_length_lyr,
                 average_chromosomes_lyr)
DEGs_karyoploter(LL_lyrG1_v_RS7G4_DEG, 
                 "lyr",
                 LL_RS7K_G4_1_lyr_real_coverage,
                 chrom_length_lyr,
                 average_chromosomes_lyr)

HM_RS7G1_v_G4_hal_DEG <- read.delim("DEGs/HM_RS7G1_v_G4_hal_DEG.txt")

DEGs_karyoploter(HM_RS7G1_v_G4_hal_DEG, 
                 "hal", 
                 HM_RS7_G4_1_hal_real_coverage,
                 chrom_length_hal,
                 average_chromosomes_hal)
