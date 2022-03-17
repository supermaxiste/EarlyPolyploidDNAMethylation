### DMRs spatial distribution along chromosomes

# Import libraries

library(tidyverse)
library(karyoploteR)

#We use karyoploteR to visualize spatial information of DMRs
# The following function takes as an input the output from
# dmrseq, the side we're considering ("halleri" or "lyrata")
# and if the reference is fine ("ok" for any comparison besides
# natural vs diploids)

DMRs_karyoploter <- function(dmrseq_output, side, ref){

  ## Filter only significant regions

  dmrseq_output <- filter(dmrseq_output, qval < 0.05)

  ## Count regions showing increase or decrease in methylation depending on ref
  ## if ref is ok it means the reference used was correct, otherwise we need to
  ## switch the interpretation the other way round

  if (ref == "ok"){
    direction <- ifelse(dmrseq_output$stat>0, "decrease", "increase")
  }
  else {
    direction <- ifelse(dmrseq_output$stat>0, "increase", "decrease")
  }



  dmrseq_output <- mutate(dmrseq_output, direction = direction)

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

  ### Define centromere coordinates

  centromere_coord <-   c(18589977,
                          5071832,
                          15417662,
                          7248175,
                          7610403,
                          14782549,
                          15405871,
                          7623002)

  if (side=="halleri"){

    ### Import agp file of the corresponding assembly

  agp <- read.delim("Ahal_mm2_tiling.agp", header=FALSE, comment.char="#")

  names(agp) <- c("new_scaffold", "start", "end", "n", "WU", "old_scaffold", "type", "length", "orientiation")

  new_dmrseq_output <- coordinate_converter(dmrseq_output, agp)

  # Splitting the scaffold column to get just numbers

  new_dmrseq_output <- separate(new_dmrseq_output, seqnames, into = c("final_scaffold", "number"), sep = "_" )
  new_dmrseq_output <- mutate(new_dmrseq_output, scaffold = new_dmrseq_output$number)


  # Increase start and end coordinates for thicker bands

  new_dmrseq_output2 <- new_dmrseq_output

  new_dmrseq_output$start <- new_dmrseq_output$end - 1000

  new_dmrseq_output$end <- new_dmrseq_output$end + 1000

  # Only black for DMRs

  new_dmrseq_output_col <- data.frame(Chrom = as.character(new_dmrseq_output$scaffold),
                                      Start = as.integer(new_dmrseq_output$start),
                                      End = as.integer(new_dmrseq_output$end),
                                      Name = as.character(new_dmrseq_output$scaffold),
                                      gieStain = "black")

  new_dmrseq_output_col2 <- data.frame(Chrom = as.character(new_dmrseq_output2$scaffold),
                                       Start = as.integer(new_dmrseq_output2$start),
                                       End = as.integer(new_dmrseq_output2$end),
                                       Name = as.character(new_dmrseq_output2$scaffold),
                                       gieStain = "black")

  # Red or blue depending on increase decrease (not expanded)

  new_dmrseq_output_col_hh2 <- data.frame(Chrom = as.character(new_dmrseq_output2$scaffold),
                                         Start = as.integer(new_dmrseq_output2$start),
                                         End = as.integer(new_dmrseq_output2$end),
                                         Name = as.character(new_dmrseq_output2$scaffold),
                                         gieStain = ifelse(new_dmrseq_output2$direction=="increase", "red", "blue"))

  new_dmrseq_output_col_red2 <- filter(new_dmrseq_output_col_hh2, gieStain=="red")
  new_dmrseq_output_col_blue2 <- filter(new_dmrseq_output_col_hh2, gieStain=="blue")

  # Red or blue depending on increase decrease (expanded)

  new_dmrseq_output_col_hh <- data.frame(Chrom = as.character(new_dmrseq_output$scaffold),
                                         Start = as.integer(new_dmrseq_output$start),
                                         End = as.integer(new_dmrseq_output$end),
                                         Name = as.character(new_dmrseq_output$scaffold),
                                         gieStain = ifelse(new_dmrseq_output$direction=="increase", "red", "blue"))

  new_dmrseq_output_col_red <- filter(new_dmrseq_output_col_hh, gieStain=="red")
  new_dmrseq_output_col_blue <- filter(new_dmrseq_output_col_hh, gieStain=="blue")

  # Define chromosome sizes halleri
  chrom_sizes <- data.frame(Chrom = c(1:8),
                            Start = c(rep(1, 8)),
                            End = c(29639309, 19005547, 26189827,
                                    24286108, 18220511, 26999834,
                                    20312336, 21625522),
                            Name = rep(NA, 8),
                            Colors = rep("lightgrey", 8))

  # Define vector with chromosome names
  chrom_names <- chrom_sizes$End
  names(chrom_names) <- chrom_sizes$Chrom

  chrom_sizes_G <- GRanges(seqnames = chrom_sizes$Chrom,
                           ranges = IRanges(start = chrom_sizes$Start,
                                            end = chrom_sizes$End),
                           seqlengths = chrom_names)


  windows <- tileGenome(chrom_names,
                        tilewidth = 1e5, cut.last.tile.in.chrom = TRUE)

  new_dmrseq_output_col_blue_G <- toGRanges(new_dmrseq_output_col_blue)
  new_dmrseq_output_col_red_G <- toGRanges(new_dmrseq_output_col_red)

  dens_red <- GenomicRanges::countOverlaps(windows, new_dmrseq_output_col_red_G)
  dens_blue <- GenomicRanges::countOverlaps(windows, new_dmrseq_output_col_blue_G)

  axis_limit <- max(c(dens_red), (dens_blue))

  complete_dens <- dens_blue + dens_red

  overall_dens <- dens_red - dens_blue

  positive_dens <- overall_dens
  positive_dens[positive_dens < 0] <- 0

  negative_dens <- overall_dens
  negative_dens[negative_dens > 0] <- 0

  no_dens <- rep(0, 1867)

  # Plot with karyoploter

  # First set parameters

  plot.params <- getDefaultPlotParams(plot.type=2)
  plot.params$data1height <- 200
  plot.params$data2height <- 200
  plot.params$ideogramheight <- 50
  plot.params$rightmargin <- 0.1

  kp <- plotKaryotype(genome = chrom_sizes_G, plot.type = 2, plot.params = plot.params)
  kpDataBackground(kp, data.panel = 1, color = "white")
  kpAxis(kp, ymax = -axis_limit, ymin = axis_limit, side = 2, data.panel = 2)
  kpDataBackground(kp, data.panel = 2, color = "white")

  new_dmrseq_output_col_red_G <- toGRanges(new_dmrseq_output_col_red)

  kpPlotRibbon(kp, data = windows, y1 = dens_red, ymin = -axis_limit, ymax = axis_limit, col = c("red"), border = NA, data.panel = 1)

  kpPlotRibbon(kp, data = windows, y1 = -dens_blue, ymin = -axis_limit, ymax = axis_limit, col = c("blue"), border = NA, data.panel = 1)

  kpPlotRibbon(kp, data = windows, y1 = no_dens, ymin = -axis_limit, ymax = axis_limit, col = c("black"), data.panel = 1)

  kpPlotRibbon(kp, data = windows, y1 = -positive_dens, ymin = -axis_limit, ymax = axis_limit, col = c("red"), border = NA, data.panel = 2)

  kpPlotRibbon(kp, data = windows, y1 = -negative_dens, ymin = -axis_limit, ymax = axis_limit, col = c("blue"), border = NA, data.panel = 2)

  kpPlotRibbon(kp, data = windows, y1 = no_dens, ymin = -axis_limit, ymax = axis_limit, col = c("black"), data.panel = 2)



  kpRect(kp, data = chrom_sizes_G, y0=0, y1=1, data.panel = "ideogram", col="white")

  kpRect(kp, data=toGRanges(new_dmrseq_output_col_blue), y0=0, y1=1, data.panel = "ideogram", col="black", border = NA)

  kpRect(kp, data=toGRanges(new_dmrseq_output_col_red), y0=0, y1=1, data.panel = "ideogram", col="black", border = NA)

  kpAbline(kp, chr = chrom_sizes$Chrom, v = centromere_coord, data.panel = 1)
  kpAbline(kp, chr = chrom_sizes$Chrom, v = centromere_coord, data.panel = 2)

  }

  else {

    agp_lyr <- read.delim("Alyr_mm2_tiling.agp", header=FALSE, comment.char="#")

    names(agp_lyr) <- c("new_scaffold", "start", "end", "n", "WU", "old_scaffold", "type", "length", "orientiation")

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
    agp <- agp_scaffold_renamer(agp_lyr)

    new_dmrseq_output <- coordinate_converter(dmrseq_output, agp)

  # Splitting the scaffold column to get just numbers

  new_dmrseq_output <- separate(new_dmrseq_output, seqnames, into = c("final_scaffold", "number"), sep = "_" )
  new_dmrseq_output <- mutate(new_dmrseq_output, scaffold = new_dmrseq_output$number)


  # Increase start and end coordinates for thicker bands

  new_dmrseq_output2 <- new_dmrseq_output

  new_dmrseq_output$start <- new_dmrseq_output$end - 1000

  new_dmrseq_output$end <- new_dmrseq_output$end + 1000

  # Only black for DMRs

  new_dmrseq_output_col <- data.frame(Chrom = as.character(new_dmrseq_output$scaffold),
                                      Start = as.integer(new_dmrseq_output$start),
                                      End = as.integer(new_dmrseq_output$end),
                                      Name = as.character(new_dmrseq_output$scaffold),
                                      gieStain = "black")

  new_dmrseq_output_col2 <- data.frame(Chrom = as.character(new_dmrseq_output2$scaffold),
                                       Start = as.integer(new_dmrseq_output2$start),
                                       End = as.integer(new_dmrseq_output2$end),
                                       Name = as.character(new_dmrseq_output2$scaffold),
                                       gieStain = "black")

  # Red or blue depending on increase decrease (not expanded)

  new_dmrseq_output_col_hh2 <- data.frame(Chrom = as.character(new_dmrseq_output2$scaffold),
                                         Start = as.integer(new_dmrseq_output2$start),
                                         End = as.integer(new_dmrseq_output2$end),
                                         Name = as.character(new_dmrseq_output2$scaffold),
                                         gieStain = ifelse(new_dmrseq_output2$direction=="increase", "red", "blue"))

  new_dmrseq_output_col_red2 <- filter(new_dmrseq_output_col_hh2, gieStain=="red")
  new_dmrseq_output_col_blue2 <- filter(new_dmrseq_output_col_hh2, gieStain=="blue")

  # Red or blue depending on increase decrease (expanded)

  new_dmrseq_output_col_hh <- data.frame(Chrom = as.character(new_dmrseq_output$scaffold),
                                         Start = as.integer(new_dmrseq_output$start),
                                         End = as.integer(new_dmrseq_output$end),
                                         Name = as.character(new_dmrseq_output$scaffold),
                                         gieStain = ifelse(new_dmrseq_output$direction=="increase", "red", "blue"))

  new_dmrseq_output_col_red <- filter(new_dmrseq_output_col_hh, gieStain=="red")
  new_dmrseq_output_col_blue <- filter(new_dmrseq_output_col_hh, gieStain=="blue")

  # Define chromosome sizes lyrata
  chrom_sizes <- data.frame(Chrom = c(1:8),
                           Start = c(rep(1, 8)),
                           End = c(25296738, 17552197, 23171854,
                                   20542238, 20888863, 19620591,
                                   22072404, 19203912),
                           Name = rep(NA, 8),
                           Colors = rep("lightgrey", 8))

  # Define vector with chromosome names
  chrom_names <- chrom_sizes$End
  names(chrom_names) <- chrom_sizes$Chrom

  chrom_sizes_G <- GRanges(seqnames = chrom_sizes$Chrom,
                           ranges = IRanges(start = chrom_sizes$Start,
                                            end = chrom_sizes$End),
                           seqlengths = chrom_names)


  windows <- tileGenome(chrom_names,
                        tilewidth = 1e5, cut.last.tile.in.chrom = TRUE)

  new_dmrseq_output_col_blue_G <- toGRanges(new_dmrseq_output_col_blue)
  new_dmrseq_output_col_red_G <- toGRanges(new_dmrseq_output_col_red)

  dens_red <- GenomicRanges::countOverlaps(windows, new_dmrseq_output_col_red_G)
  dens_blue <- GenomicRanges::countOverlaps(windows, new_dmrseq_output_col_blue_G)

  axis_limit <- max(c(dens_red), (dens_blue))

  complete_dens <- dens_blue + dens_red

  overall_dens <- dens_red - dens_blue

  positive_dens <- overall_dens
  positive_dens[positive_dens < 0] <- 0

  negative_dens <- overall_dens
  negative_dens[negative_dens > 0] <- 0

  no_dens <- rep(0, 1867)

  # Plot with karyoploter

  # First set parameters

  plot.params <- getDefaultPlotParams(plot.type=2)
  plot.params$data1height <- 200
  plot.params$data2height <- 200
  plot.params$ideogramheight <- 50
  plot.params$rightmargin <- 0.1

  kp <- plotKaryotype(genome = chrom_sizes_G, plot.type = 2, plot.params = plot.params)
  kpDataBackground(kp, data.panel = 1, color = "white")
  kpAxis(kp, ymax = -axis_limit, ymin = axis_limit, side = 2, data.panel = 2)
  kpDataBackground(kp, data.panel = 2, color = "white")

  new_dmrseq_output_col_red_G <- toGRanges(new_dmrseq_output_col_red)

  kpPlotRibbon(kp, data = windows, y1 = dens_red, ymin = -axis_limit, ymax = axis_limit, col = c("red"), border = NA, data.panel = 1)

  kpPlotRibbon(kp, data = windows, y1 = -dens_blue, ymin = -axis_limit, ymax = axis_limit, col = c("blue"), border = NA, data.panel = 1)

  kpPlotRibbon(kp, data = windows, y1 = no_dens, ymin = -axis_limit, ymax = axis_limit, col = c("black"), data.panel = 1)

  kpPlotRibbon(kp, data = windows, y1 = -positive_dens, ymin = -axis_limit, ymax = axis_limit, col = c("red"), border = NA, data.panel = 2)

  kpPlotRibbon(kp, data = windows, y1 = -negative_dens, ymin = -axis_limit, ymax = axis_limit, col = c("blue"), border = NA, data.panel = 2)

  kpPlotRibbon(kp, data = windows, y1 = no_dens, ymin = -axis_limit, ymax = axis_limit, col = c("black"), data.panel = 2)



  kpRect(kp, data = chrom_sizes_G, y0=0, y1=1, data.panel = "ideogram", col="white")

  kpRect(kp, data=toGRanges(new_dmrseq_output_col_blue), y0=0, y1=1, data.panel = "ideogram", col="black", border = NA)

  kpRect(kp, data=toGRanges(new_dmrseq_output_col_red), y0=0, y1=1, data.panel = "ideogram", col="black", border = NA)

  kpAbline(kp, chr = chrom_sizes$Chrom, v = centromere_coord, data.panel = 1)
  kpAbline(kp, chr = chrom_sizes$Chrom, v = centromere_coord, data.panel = 2)

  }
}

DMRs_karyoploter_special <- function(dmrseq_output_CG_AvB,
                                     dmrseq_output_CHG_AvB,
                                     dmrseq_output_CHH_AvB, ref){
  scaffold_number_CG <- as.integer(str_split_fixed(dmrseq_output_CG_AvB$seqnames, "_", 2)[,2])
  scaffold_number_CHG <- as.integer(str_split_fixed(dmrseq_output_CHG_AvB$seqnames, "_", 2)[,2])
  scaffold_number_CHH <- as.integer(str_split_fixed(dmrseq_output_CHH_AvB$seqnames, "_", 2)[,2])
  dmrseq_output_CG_hal <- dmrseq_output_CG_AvB[scaffold_number_CG < 2240]
  dmrseq_output_CG_lyr <- dmrseq_output_CG_AvB[scaffold_number_CG > 2239]
  dmrseq_output_CHG_hal <- dmrseq_output_CHG_AvB[scaffold_number_CHG < 2240]
  dmrseq_output_CHG_lyr <- dmrseq_output_CHG_AvB[scaffold_number_CHG > 2240]
  dmrseq_output_CHH_hal <- dmrseq_output_CHH_AvB[scaffold_number_CHH < 2240]
  dmrseq_output_CHH_lyr <- dmrseq_output_CHH_AvB[scaffold_number_CHH > 2240]
  dmrseq_output <- rbind(dmrseq_output_CG_hal, dmrseq_output_CHG_hal, dmrseq_output_CHH_hal)
  DMRs_karyoploter(dmrseq_output = dmrseq_output, side = "halleri", ref)
  dmrseq_output <- rbind(dmrseq_output_CG_lyr, dmrseq_output_CHG_lyr, dmrseq_output_CHH_lyr)
  DMRs_karyoploter(dmrseq_output = dmrseq_output, side = "lyrata", ref)
}

# Import datasets and plot 

# HM_syn1Vpro1

setwd("ARPEGGIO_results_syn1Vpro1_HM/dmrseq")

dmrseq_output_CG_hal <- fread("CG_context/parent1_v_allo.txt")
dmrseq_output_CHG_hal <- fread("CHG_context/parent1_v_allo.txt")
dmrseq_output_CHH_hal <- fread("CHH_context/parent1_v_allo.txt")

dmrseq_output_all_hal <- rbind(dmrseq_output_CG_hal,
                           dmrseq_output_CHG_hal,
                           dmrseq_output_CHH_hal)

dmrseq_output_CG_lyr <- fread("CG_context/parent2_v_allo.txt")
dmrseq_output_CHG_lyr <- fread("CHG_context/parent2_v_allo.txt")
dmrseq_output_CHH_lyr <- fread("CHH_context/parent2_v_allo.txt")

dmrseq_output_all_lyr <- rbind(dmrseq_output_CG_lyr,
                               dmrseq_output_CHG_lyr,
                               dmrseq_output_CHH_lyr)

HM_syn1_v_hal1 <- DMRs_karyoploter(dmrseq_output = dmrseq_output_all_hal,
                 side = "halleri",
                 ref = "ok")

HM_syn1_v_lyr1 <- DMRs_karyoploter(dmrseq_output = dmrseq_output_all_lyr,
                                   side = "lyrata",
                                   ref = "ok")

# HM_syn4Vpro1
setwd("ARPEGGIO_results_syn4Vpro1_HM/dmrseq")

dmrseq_output_CG_hal <- fread("CG_context/parent1_v_allo.txt")
dmrseq_output_CHG_hal <- fread("CHG_context/parent1_v_allo.txt")
dmrseq_output_CHH_hal <- fread("CHH_context/parent1_v_allo.txt")

dmrseq_output_all_hal <- rbind(dmrseq_output_CG_hal,
                               dmrseq_output_CHG_hal,
                               dmrseq_output_CHH_hal)

dmrseq_output_CG_lyr <- fread("CG_context/parent2_v_allo.txt")
dmrseq_output_CHG_lyr <- fread("CHG_context/parent2_v_allo.txt")
dmrseq_output_CHH_lyr <- fread("CHH_context/parent2_v_allo.txt")

dmrseq_output_all_lyr <- rbind(dmrseq_output_CG_lyr,
                               dmrseq_output_CHG_lyr,
                               dmrseq_output_CHH_lyr)

HM_syn4_v_hal1 <- DMRs_karyoploter(dmrseq_output = dmrseq_output_all_hal,
                                   side = "halleri",
                                   ref = "ok")

HM_syn4_v_lyr1 <- DMRs_karyoploter(dmrseq_output = dmrseq_output_all_lyr,
                                   side = "lyrata",
                                   ref = "ok")


# LL_syn1Vpro1

setwd("ARPEGGIO_results_syn1Vpro1_LL/dmrseq")

dmrseq_output_CG_hal <- fread("CG_context/parent1_v_allo.txt")
dmrseq_output_CHG_hal <- fread("CHG_context/parent1_v_allo.txt")
dmrseq_output_CHH_hal <- fread("CHH_context/parent1_v_allo.txt")

dmrseq_output_all_hal <- rbind(dmrseq_output_CG_hal,
                               dmrseq_output_CHG_hal,
                               dmrseq_output_CHH_hal)

dmrseq_output_CG_lyr <- fread("CG_context/parent2_v_allo.txt")
dmrseq_output_CHG_lyr <- fread("CHG_context/parent2_v_allo.txt")
dmrseq_output_CHH_lyr <- fread("CHH_context/parent2_v_allo.txt")

dmrseq_output_all_lyr <- rbind(dmrseq_output_CG_lyr,
                               dmrseq_output_CHG_lyr,
                               dmrseq_output_CHH_lyr)

LL_syn1_v_hal1 <- DMRs_karyoploter(dmrseq_output = dmrseq_output_all_hal,
                                   side = "halleri",
                                   ref = "ok")

LL_syn1_v_lyr1 <- DMRs_karyoploter(dmrseq_output = dmrseq_output_all_lyr,
                                   side = "lyrata",
                                   ref = "ok")

# LL_syn4Vpro1
setwd("ARPEGGIO_results_syn4Vpro1_LL/dmrseq")

dmrseq_output_CG_hal <- fread("CG_context/parent1_v_allo.txt")
dmrseq_output_CHG_hal <- fread("CHG_context/parent1_v_allo.txt")
dmrseq_output_CHH_hal <- fread("CHH_context/parent1_v_allo.txt")

dmrseq_output_all_hal <- rbind(dmrseq_output_CG_hal,
                               dmrseq_output_CHG_hal,
                               dmrseq_output_CHH_hal)

dmrseq_output_CG_lyr <- fread("CG_context/parent2_v_allo.txt")
dmrseq_output_CHG_lyr <- fread("CHG_context/parent2_v_allo.txt")
dmrseq_output_CHH_lyr <- fread("CHH_context/parent2_v_allo.txt")

dmrseq_output_all_lyr <- rbind(dmrseq_output_CG_lyr,
                               dmrseq_output_CHG_lyr,
                               dmrseq_output_CHH_lyr)

LL_syn4_v_hal1 <- DMRs_karyoploter(dmrseq_output = dmrseq_output_all_hal,
                                   side = "halleri",
                                   ref = "ok")

LL_syn4_v_lyr1 <- DMRs_karyoploter(dmrseq_output = dmrseq_output_all_lyr,
                                   side = "lyrata",
                                   ref = "ok")
