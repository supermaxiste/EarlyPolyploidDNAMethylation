## In this script we will try to take a look at
## the results from remapping to see if the 
## three methods we used have similar results
## or not

## Trials for ggbio and karyoploteR
## too difficult for doing what I'd like
## Packages tried without success: karyoploteR, ggbio.
## Successful packages: chromPlot and circlize
## Dotplots successfully done with: D-GENIES (http://dgenies.toulouse.inra.fr/)
## Alternative packages/software not tried: 
## MCScanX (https://github.com/wyp1125/MCScanX) and minidot (https://github.com/thackl/minidot)

library(GenomicRanges)
library(tidyverse)


####### HALLERI #######

###### MUMMER ########

Ahal_mum_tiling <- read.delim("Ahal_mum_tiling.agp", header=FALSE, comment.char="#")
colnames(Ahal_mum_tiling) <- c("final_scaffold", "start", "end", "n", "contig_or_gap", "query_scaffold", "query_scaffold_start", "query_scaffold_end", "orientation")

# Formatting the data to have the same format as hg_CytoBandIdeo data
# 5 columns with the following variables: Chrom, Start, End, Name, gieStain

# Removing gaps

Ahal_mum_tiling <- filter(Ahal_mum_tiling, query_scaffold!=100)

# Selecting only 8 chromosomes

Ahal_mum_tiling_8 <- Ahal_mum_tiling %>%  dplyr::filter(final_scaffold=="scaffold_1_RagTag" | 
                                                          final_scaffold=="scaffold_2_RagTag" | final_scaffold=="scaffold_3_RagTag" |
                                                          final_scaffold=="scaffold_4_RagTag" | final_scaffold=="scaffold_5_RagTag" |
                                                          final_scaffold=="scaffold_6_RagTag" | final_scaffold=="scaffold_7_RagTag" |
                                                          final_scaffold=="scaffold_8_RagTag")

# Splitting the scaffold column to get just numbers

Ahal_mum_tiling_8_clean <- separate(Ahal_mum_tiling_8, final_scaffold, into = c("final_scaffold", "number"), sep = "_" )
Ahal_mum_tiling_8_clean <- mutate(Ahal_mum_tiling_8_clean, scaffold = Ahal_mum_tiling_8_clean$number)
Ahal_mum_tiling_8_clean2 <- select(Ahal_mum_tiling_8_clean, scaffold, start, end, query_scaffold)

# Create final data frame with all the info as hg_CytoBandIdeo data (see above)

Ahal_mum_tiling_8G <- data.frame(Chrom = as.character(Ahal_mum_tiling_8_clean2$scaffold), 
                                 Start = as.integer(Ahal_mum_tiling_8_clean2$start),
                                 End = as.integer(Ahal_mum_tiling_8_clean2$end),
                                 Name = as.character(Ahal_mum_tiling_8_clean2$query_scaffold),
                                 gieStain = as.character(Ahal_mum_tiling_8_clean2$query_scaffold))

# Now we colour the chromosomes with the results from the mappers
# First we colour in 8 different colours for mummer mappings

chrom_colors_mum <- c(rep("red", sum(Ahal_mum_tiling_8G$Chrom=="1")), 
                  rep("green", sum(Ahal_mum_tiling_8G$Chrom=="2")),
                  rep("blue", sum(Ahal_mum_tiling_8G$Chrom=="3")),
                  rep("yellow", sum(Ahal_mum_tiling_8G$Chrom=="4")),
                  rep("brown", sum(Ahal_mum_tiling_8G$Chrom=="5")),
                  rep("purple", sum(Ahal_mum_tiling_8G$Chrom=="6")),
                  rep("turquoise", sum(Ahal_mum_tiling_8G$Chrom=="7")),
                  rep("pink", sum(Ahal_mum_tiling_8G$Chrom=="8")))

Ahal_mum_tiling_8G_col <- data.frame(Chrom = as.character(Ahal_mum_tiling_8_clean2$scaffold), 
                                     Start = as.integer(Ahal_mum_tiling_8_clean2$start),
                                     End = as.integer(Ahal_mum_tiling_8_clean2$end),
                                     Name = as.character(Ahal_mum_tiling_8_clean2$query_scaffold),
                                     Colors = chrom_colors_mum)

chrom_sizes <- data.frame(Chrom = c(1:8), 
                          Start = c(rep(1, 8)), 
                          End = c(33132539, 19320864, 24464547,
                                  23328337, 21221946, 25113588,
                                  24649197, 22951293),
                          Name = rep(NA, 8),
                          Colors = rep("lightgrey", 8))

Ahal_mum_tiling_8G_col <- rbind(chrom_sizes, Ahal_mum_tiling_8G_col)


Ahal_mm2_tiling <- read.delim("Ahal_mm2_tiling.agp", header=FALSE, comment.char="#")
colnames(Ahal_mm2_tiling) <- c("final_scaffold", "start", "end", "n", "contig_or_gap", "query_scaffold", "query_scaffold_start", "query_scaffold_end", "orientation")

# Formatting the data to have the same format as hg_CytoBandIdeo data
# 5 columns with the following variables: Chrom, Start, End, Name, gieStain

# Removing gaps

Ahal_mm2_tiling <- filter(Ahal_mm2_tiling, query_scaffold!=100)

# Selecting only 8 chromosomes

Ahal_mm2_tiling_8 <- Ahal_mm2_tiling %>%  dplyr::filter(final_scaffold=="scaffold_1_RagTag" | 
                                                          final_scaffold=="scaffold_2_RagTag" | final_scaffold=="scaffold_3_RagTag" |
                                                          final_scaffold=="scaffold_4_RagTag" | final_scaffold=="scaffold_5_RagTag" |
                                                          final_scaffold=="scaffold_6_RagTag" | final_scaffold=="scaffold_7_RagTag" |
                                                          final_scaffold=="scaffold_8_RagTag")

# Splitting the scaffold column to get just numbers
Ahal_mm2_tiling_8_clean <- separate(Ahal_mm2_tiling_8, final_scaffold, into = c("final_scaffold", "number"), sep = "_" )
Ahal_mm2_tiling_8_clean <- mutate(Ahal_mm2_tiling_8_clean, scaffold = Ahal_mm2_tiling_8_clean$number)
Ahal_mm2_tiling_8_clean2 <- select(Ahal_mm2_tiling_8_clean, scaffold, start, end, query_scaffold)

# Create final data frame with all the info as hg_CytoBandIdeo data (see above)

Ahal_mm2_tiling_8G <- data.frame(Chrom = as.character(Ahal_mm2_tiling_8_clean2$scaffold), 
                                 Start = as.integer(Ahal_mm2_tiling_8_clean2$start),
                                 End = as.integer(Ahal_mm2_tiling_8_clean2$end),
                                 Name = as.character(Ahal_mm2_tiling_8_clean2$query_scaffold),
                                 gieStain = as.character(Ahal_mm2_tiling_8_clean2$query_scaffold))

Ahal_mm2_tiling_8G_col <- select(Ahal_mm2_tiling_8G, Chrom, Start, End, Name)

# Now we match the colors set for mummer to minimap's

matching_colors_mm2 <- match(Ahal_mm2_tiling_8_clean2$query_scaffold, Ahal_mum_tiling_8_clean2$query_scaffold)

chrom_colors_mm2 <- Ahal_mum_tiling_8G_col[matching_colors_mm2, "Colors"]


Ahal_mm2_tiling_8G_col <- mutate(Ahal_mm2_tiling_8G_col, Colors = chrom_colors_mm2)

Ahal_mm2_tiling_8G_col$Colors[is.na(Ahal_mm2_tiling_8G_col$Colors)] <- "lightgrey"

## Lastal

Ahal_lst_tiling <- read.delim("Ahal_last_tiling.agp", header=FALSE, comment.char="#")
colnames(Ahal_lst_tiling) <- c("final_scaffold", "start", "end", "n", "contig_or_gap", "query_scaffold", "query_scaffold_start", "query_scaffold_end", "orientation")

# Formatting the data to have the same format as hg_CytoBandIdeo data
# 5 columns with the following variables: Chrom, Start, End, Name, gieStain

# Removing gaps

Ahal_lst_tiling <- filter(Ahal_lst_tiling, query_scaffold!=100)

# Selecting only 8 chromosomes

Ahal_lst_tiling_8 <- Ahal_lst_tiling %>%  dplyr::filter(final_scaffold=="scaffold_1_RagTag" |
                                                          final_scaffold=="scaffold_2_RagTag" | final_scaffold=="scaffold_3_RagTag" |
                                                          final_scaffold=="scaffold_4_RagTag" | final_scaffold=="scaffold_5_RagTag" |
                                                          final_scaffold=="scaffold_6_RagTag" | final_scaffold=="scaffold_7_RagTag" |
                                                          final_scaffold=="scaffold_8_RagTag")

# Splitting the scaffold column to get just numbers

Ahal_lst_tiling_8_clean <- separate(Ahal_lst_tiling_8, final_scaffold, into = c("final_scaffold", "number"), sep = "_" )
Ahal_lst_tiling_8_clean <- mutate(Ahal_lst_tiling_8_clean, scaffold = Ahal_lst_tiling_8_clean$number)
Ahal_lst_tiling_8_clean2 <- select(Ahal_lst_tiling_8_clean, scaffold, start, end, query_scaffold)

# Now we match the colors set for minimap2 to mummer's

matching_colors_lst <- match(Ahal_lst_tiling_8_clean2$query_scaffold, Ahal_mm2_tiling_8_clean2$query_scaffold)

original_colors <- Ahal_mm2_tiling_8G_col$Colors

new_colors <- na.omit(Ahal_mm2_tiling_8G_col[matching_colors_lst, "Colors"])

# We create a new vector from the old colors and move some
# around depending on their mapping position

original_colors_lst <- original_colors

original_colors_lst[na.omit(matching_colors_lst)] <- new_colors

# Turn into grey the ones not matching

original_colors_lst[which(is.na(matching_colors_lst))] <- "lightgrey"

# Turn into grey the unmapped too

unmapped_colors_lst <- which(is.na(match(Ahal_mm2_tiling_8_clean2$query_scaffold, Ahal_lst_tiling_8_clean2$query_scaffold)))
original_colors_lst[unmapped_colors_lst] <- "lightgrey"

Ahal_lst_tiling_8G_col <- data.frame(Chrom = as.character(Ahal_mm2_tiling_8_clean2$scaffold),
                                     Start = as.integer(Ahal_mm2_tiling_8_clean2$start),
                                     End = as.integer(Ahal_mm2_tiling_8_clean2$end),
                                     Name = as.character(Ahal_mm2_tiling_8_clean2$query_scaffold),
                                     Colors = original_colors_lst)

# Now we plot synteny with circular plots

library(circlize)

############ halleri side ###############


chromosomes <- as.character(c(1:16))
scaffold_length <- as.matrix(cbind(c(1:16), c(33132539, 19320864, 24464547,
                     23328337, 21221946, 25113588, 24649197, 22951293, 
                     33132539, 19320864, 24464547, 23328337, 21221946, 
                     25113588, 24649197, 22951293)))

# We inialize the plot with chromosome names and scaffold lengths

circos.initialize(factors = chromosomes, xlim = scaffold_length)

# We show the outcome

circos.track(chromosomes, ylim = c(0, 1))

# We thin the tracks

circos.par(track.height = 0.1)

# circos.par(RESET = TRUE)

circos.initialize(factors = chromosomes, xlim = scaffold_length)

circos.track(chromosomes, ylim = c(0, 1))


# To draw links:
# bed1 = generateRandomBed(nr = 100)
# bed1 = bed1[sample(nrow(bed1), 20), ]
# bed2 = generateRandomBed(nr = 100)
# bed2 = bed2[sample(nrow(bed2), 20), ]
# circos.genomicLink(bed1, bed2, col = rand_color(nrow(bed1), transparency = 0.5),  border = NA)

bed1 <- select(Ahal_mum_tiling_8_clean2, chr = scaffold, start, end)

bed1 <- cbind(bed1, value1 = rep(0, nrow(bed1)))

bed2 <- select(Ahal_mm2_tiling_8_clean2, chr = scaffold, start, end)

bed2 <- cbind(bed2, value1 = rep(0, nrow(bed2)))

bed3 <- select(Ahal_lst_tiling_8_clean2, chr = scaffold, start, end)

bed3 <- cbind(bed3, value1 = rep(0, nrow(bed3)))

## To plot mm2 vs the other two tools, we add 8 to the scaffold # for
## coordinate matching. mum and lst will have the first 8 "bands" while 
## mm2 will have the spots from 9 to 16

bed2$chr <- as.character(as.integer(bed2$chr) + 8)

## We match scaffolds to find out the links between them

match_mum_mm2 <- match(Ahal_mum_tiling_8_clean2$query_scaffold, Ahal_mm2_tiling_8_clean2$query_scaffold)

match_lst_mm2 <- match(Ahal_lst_tiling_8_clean2$query_scaffold, Ahal_mm2_tiling_8_clean2$query_scaffold)

## mm2 vs mum

bed2_filter <- bed2[na.omit(match_mum_mm2),]

bed1_filter <- bed1[!is.na(match_mum_mm2),]

link_colors <- c(rep("red", sum(bed1_filter$chr=="1")),
                 rep("green", sum(bed1_filter$chr=="2")),
                 rep("blue", sum(bed1_filter$chr=="3")),
                 rep("yellow", sum(bed1_filter$chr=="4")),
                 rep("brown", sum(bed1_filter$chr=="5")),
                 rep("purple", sum(bed1_filter$chr=="6")),
                 rep("turquoise", sum(bed1_filter$chr=="7")),
                 rep("pink", sum(bed1_filter$chr=="8")))

circos.genomicLink(bed1_filter, bed2_filter, col = link_colors,  border = NA)

circos.text(scaffold_length[1,2]/2, 0.5, "1", 1, 1, facing = "bending.inside")
circos.text(scaffold_length[2,2]/2, 0.5, "2", 2, 1, facing = "bending.inside")
circos.text(scaffold_length[3,2]/2, 0.5, "3", 3, 1, facing = "bending.inside")
circos.text(scaffold_length[4,2]/2, 0.5, "4", 4, 1, facing = "bending.inside")
circos.text(scaffold_length[5,2]/2, 0.5, "5", 5, 1, facing = "bending.inside")
circos.text(scaffold_length[6,2]/2, 0.5, "6", 6, 1, facing = "bending.inside")
circos.text(scaffold_length[7,2]/2, 0.5, "7", 7, 1, facing = "bending.inside")
circos.text(scaffold_length[8,2]/2, 0.5, "8", 8, 1, facing = "bending.inside")

circos.text(scaffold_length[1,2]/2, 0.5, "1", 9, 1, facing = "bending.inside")
circos.text(scaffold_length[2,2]/2, 0.5, "2", 10, 1, facing = "bending.inside")
circos.text(scaffold_length[3,2]/2, 0.5, "3", 11, 1, facing = "bending.inside")
circos.text(scaffold_length[4,2]/2, 0.5, "4", 12, 1, facing = "bending.inside")
circos.text(scaffold_length[5,2]/2, 0.5, "5", 13, 1, facing = "bending.inside")
circos.text(scaffold_length[6,2]/2, 0.5, "6", 14, 1, facing = "bending.inside")
circos.text(scaffold_length[7,2]/2, 0.5, "7", 15, 1, facing = "bending.inside")
circos.text(scaffold_length[8,2]/2, 0.5, "8", 16, 1, facing = "bending.inside")

circos.text(scaffold_length[4,2], 2.2, "mummer", 4, 1, facing = "bending.inside", cex = 2)
circos.text(scaffold_length[12,2], 2.2, "minimap2", 12, 1, facing = "bending.inside", cex = 2)

## mm2 vs lst

circos.initialize(factors = chromosomes, xlim = scaffold_length)

circos.track(chromosomes, ylim = c(0, 1))

bed2_filter <- bed2[na.omit(match_lst_mm2),]

bed3_filter <- bed3[!is.na(match_lst_mm2),]

link_colors <- c(rep("red", sum(bed3_filter$chr=="1")),
                 rep("green", sum(bed3_filter$chr=="2")),
                 rep("blue", sum(bed3_filter$chr=="3")),
                 rep("yellow", sum(bed3_filter$chr=="4")),
                 rep("brown", sum(bed3_filter$chr=="5")),
                 rep("purple", sum(bed3_filter$chr=="6")),
                 rep("turquoise", sum(bed3_filter$chr=="7")),
                 rep("pink", sum(bed3_filter$chr=="8")))

circos.genomicLink(bed3_filter, bed2_filter, col = link_colors,  border = NA)

circos.text(scaffold_length[1,2]/2, 0.5, "1", 1, 1, facing = "bending.inside")
circos.text(scaffold_length[2,2]/2, 0.5, "2", 2, 1, facing = "bending.inside")
circos.text(scaffold_length[3,2]/2, 0.5, "3", 3, 1, facing = "bending.inside")
circos.text(scaffold_length[4,2]/2, 0.5, "4", 4, 1, facing = "bending.inside")
circos.text(scaffold_length[5,2]/2, 0.5, "5", 5, 1, facing = "bending.inside")
circos.text(scaffold_length[6,2]/2, 0.5, "6", 6, 1, facing = "bending.inside")
circos.text(scaffold_length[7,2]/2, 0.5, "7", 7, 1, facing = "bending.inside")
circos.text(scaffold_length[8,2]/2, 0.5, "8", 8, 1, facing = "bending.inside")

circos.text(scaffold_length[1,2]/2, 0.5, "1", 9, 1, facing = "bending.inside")
circos.text(scaffold_length[2,2]/2, 0.5, "2", 10, 1, facing = "bending.inside")
circos.text(scaffold_length[3,2]/2, 0.5, "3", 11, 1, facing = "bending.inside")
circos.text(scaffold_length[4,2]/2, 0.5, "4", 12, 1, facing = "bending.inside")
circos.text(scaffold_length[5,2]/2, 0.5, "5", 13, 1, facing = "bending.inside")
circos.text(scaffold_length[6,2]/2, 0.5, "6", 14, 1, facing = "bending.inside")
circos.text(scaffold_length[7,2]/2, 0.5, "7", 15, 1, facing = "bending.inside")
circos.text(scaffold_length[8,2]/2, 0.5, "8", 16, 1, facing = "bending.inside")

circos.text(scaffold_length[4,2], 2.2, "last", 4, 1, facing = "bending.inside", cex = 2)
circos.text(scaffold_length[12,2], 2.2, "minimap2", 12, 1, facing = "bending.inside", cex = 2)

## mum vs lst

bed1$chr <- as.character(as.integer(bed1$chr) + 8)

## We match scaffolds to find out the links between them

match_lst_mum <- match(Ahal_lst_tiling_8_clean2$query_scaffold, Ahal_mum_tiling_8_clean2$query_scaffold)

circos.initialize(factors = chromosomes, xlim = scaffold_length)

circos.track(chromosomes, ylim = c(0, 1))

bed1_filter <- bed1[na.omit(match_lst_mum),]

bed3_filter <- bed3[!is.na(match_lst_mum),]

link_colors <- c(rep("red", sum(bed3_filter$chr=="1")),
                 rep("green", sum(bed3_filter$chr=="2")),
                 rep("blue", sum(bed3_filter$chr=="3")),
                 rep("yellow", sum(bed3_filter$chr=="4")),
                 rep("brown", sum(bed3_filter$chr=="5")),
                 rep("purple", sum(bed3_filter$chr=="6")),
                 rep("turquoise", sum(bed3_filter$chr=="7")),
                 rep("pink", sum(bed3_filter$chr=="8")))

circos.genomicLink(bed3_filter, bed1_filter, col = link_colors,  border = NA)

circos.text(scaffold_length[1,2]/2, 0.5, "1", 1, 1, facing = "bending.inside")
circos.text(scaffold_length[2,2]/2, 0.5, "2", 2, 1, facing = "bending.inside")
circos.text(scaffold_length[3,2]/2, 0.5, "3", 3, 1, facing = "bending.inside")
circos.text(scaffold_length[4,2]/2, 0.5, "4", 4, 1, facing = "bending.inside")
circos.text(scaffold_length[5,2]/2, 0.5, "5", 5, 1, facing = "bending.inside")
circos.text(scaffold_length[6,2]/2, 0.5, "6", 6, 1, facing = "bending.inside")
circos.text(scaffold_length[7,2]/2, 0.5, "7", 7, 1, facing = "bending.inside")
circos.text(scaffold_length[8,2]/2, 0.5, "8", 8, 1, facing = "bending.inside")

circos.text(scaffold_length[1,2]/2, 0.5, "1", 9, 1, facing = "bending.inside")
circos.text(scaffold_length[2,2]/2, 0.5, "2", 10, 1, facing = "bending.inside")
circos.text(scaffold_length[3,2]/2, 0.5, "3", 11, 1, facing = "bending.inside")
circos.text(scaffold_length[4,2]/2, 0.5, "4", 12, 1, facing = "bending.inside")
circos.text(scaffold_length[5,2]/2, 0.5, "5", 13, 1, facing = "bending.inside")
circos.text(scaffold_length[6,2]/2, 0.5, "6", 14, 1, facing = "bending.inside")
circos.text(scaffold_length[7,2]/2, 0.5, "7", 15, 1, facing = "bending.inside")
circos.text(scaffold_length[8,2]/2, 0.5, "8", 16, 1, facing = "bending.inside")

circos.text(scaffold_length[4,2], 2.2, "last", 4, 1, facing = "bending.inside", cex = 2)
circos.text(scaffold_length[12,2], 2.2, "mummer", 12, 1, facing = "bending.inside", cex = 2)


########## lyrata side #############

# Read files

Alyr_mm2_tiling <- read.delim("Alyr_mm2_tiling.agp", header=FALSE, comment.char="#")
colnames(Alyr_mm2_tiling) <- c("final_scaffold", "start", "end", "n", "contig_or_gap", "query_scaffold", "query_scaffold_start", "query_scaffold_end", "orientation")

# Formatting the data to have the same format as hg_CytoBandIdeo data
# 5 columns with the following variables: Chrom, Start, End, Name, gieStain

# Removing gaps

Alyr_mm2_tiling <- filter(Alyr_mm2_tiling, query_scaffold!=100)

# Selecting only 8 chromosomes

Alyr_mm2_tiling_8 <- Alyr_mm2_tiling %>%  dplyr::filter(final_scaffold=="scaffold_1_RagTag" | 
                                                          final_scaffold=="scaffold_2_RagTag" | final_scaffold=="scaffold_3_RagTag" |
                                                          final_scaffold=="scaffold_4_RagTag" | final_scaffold=="scaffold_5_RagTag" |
                                                          final_scaffold=="scaffold_6_RagTag" | final_scaffold=="scaffold_7_RagTag" |
                                                          final_scaffold=="scaffold_8_RagTag")

# Splitting the scaffold column to get just numbers
Alyr_mm2_tiling_8_clean <- separate(Alyr_mm2_tiling_8, final_scaffold, into = c("final_scaffold", "number"), sep = "_" )
Alyr_mm2_tiling_8_clean <- mutate(Alyr_mm2_tiling_8_clean, scaffold = Alyr_mm2_tiling_8_clean$number)
Alyr_mm2_tiling_8_clean2 <- select(Alyr_mm2_tiling_8_clean, scaffold, start, end, query_scaffold)

# Create final data frame with all the info as hg_CytoBandIdeo data (see above)

Alyr_mm2_tiling_8G <- data.frame(Chrom = as.character(Alyr_mm2_tiling_8_clean2$scaffold), 
                                 Start = as.integer(Alyr_mm2_tiling_8_clean2$start),
                                 End = as.integer(Alyr_mm2_tiling_8_clean2$end),
                                 Name = as.character(Alyr_mm2_tiling_8_clean2$query_scaffold),
                                 gieStain = as.character(Alyr_mm2_tiling_8_clean2$query_scaffold))

# Plotting the chromosomes

#chromPlot(gaps = gap_file, bands = Alyr_mm2_tiling_8G)

# Now we colour the chromosomes with the results from the mappers
# First we colour in 8 different colours for minimap2 mappings

Alyr_mm2_tiling_8G_col <- dplyr::select(Alyr_mm2_tiling_8G, Chrom, Start, End, Name)

chrom_colors <- c(rep("red", sum(Alyr_mm2_tiling_8G$Chrom=="1")), 
                  rep("green", sum(Alyr_mm2_tiling_8G$Chrom=="2")),
                  rep("blue", sum(Alyr_mm2_tiling_8G$Chrom=="3")),
                  rep("yellow", sum(Alyr_mm2_tiling_8G$Chrom=="4")),
                  rep("brown", sum(Alyr_mm2_tiling_8G$Chrom=="5")),
                  rep("purple", sum(Alyr_mm2_tiling_8G$Chrom=="6")),
                  rep("turquoise", sum(Alyr_mm2_tiling_8G$Chrom=="7")),
                  rep("pink", sum(Alyr_mm2_tiling_8G$Chrom=="8")))

Alyr_mm2_tiling_8G_col <- mutate(Alyr_mm2_tiling_8G_col, Colors = chrom_colors)

chromPlot(bands = Alyr_mm2_tiling_8G_col)

## The first version of the plot will keep all the colours but I can also make
## one plot for each chromosome

## Mummer

Alyr_mum_tiling <- read.delim("Alyr_mum_tiling.agp", header=FALSE, comment.char="#")
colnames(Alyr_mum_tiling) <- c("final_scaffold", "start", "end", "n", "contig_or_gap", "query_scaffold", "query_scaffold_start", "query_scaffold_end", "orientation")

# Formatting the data to have the same format as hg_CytoBandIdeo data
# 5 columns with the following variables: Chrom, Start, End, Name, gieStain

# Removing gaps

Alyr_mum_tiling <- filter(Alyr_mum_tiling, query_scaffold!=100)

# Selecting only 8 chromosomes

Alyr_mum_tiling_8 <- Alyr_mum_tiling %>%  dplyr::filter(final_scaffold=="scaffold_1_RagTag" | 
                                                          final_scaffold=="scaffold_2_RagTag" | final_scaffold=="scaffold_3_RagTag" |
                                                          final_scaffold=="scaffold_4_RagTag" | final_scaffold=="scaffold_5_RagTag" |
                                                          final_scaffold=="scaffold_6_RagTag" | final_scaffold=="scaffold_7_RagTag" |
                                                          final_scaffold=="scaffold_8_RagTag")

# Splitting the scaffold column to get just numbers

Alyr_mum_tiling_8_clean <- separate(Alyr_mum_tiling_8, final_scaffold, into = c("final_scaffold", "number"), sep = "_" )
Alyr_mum_tiling_8_clean <- mutate(Alyr_mum_tiling_8_clean, scaffold = Alyr_mum_tiling_8_clean$number)
Alyr_mum_tiling_8_clean2 <- dplyr::select(Alyr_mum_tiling_8_clean, scaffold, start, end, query_scaffold)

# Now we match the colors set for minimap2 to mummer's

matching_colors_mum <- match(Alyr_mum_tiling_8_clean2$query_scaffold, Alyr_mm2_tiling_8_clean2$query_scaffold)

original_colors <- Alyr_mm2_tiling_8G_col$Colors

new_colors <- na.omit(Alyr_mm2_tiling_8G_col[matching_colors_mum, "Colors"])

# We create a new vector from the old colors and move some
# around depending on their mapping position

original_colors_mum <- original_colors

original_colors_mum[na.omit(matching_colors_mum)] <- new_colors

# Turn into grey the ones not matching

original_colors_mum[which(is.na(matching_colors_mum))] <- "lightgrey"

# Turn into grey the unmapped too

unmapped_colors_mum <- which(is.na(match(Alyr_mm2_tiling_8_clean2$query_scaffold, Alyr_mum_tiling_8_clean2$query_scaffold)))
original_colors_mum[unmapped_colors_mum] <- "lightgrey"

Alyr_mum_tiling_8G_col <- data.frame(Chrom = as.character(Alyr_mm2_tiling_8_clean2$scaffold), 
                                     Start = as.integer(Alyr_mm2_tiling_8_clean2$start),
                                     End = as.integer(Alyr_mm2_tiling_8_clean2$end),
                                     Name = as.character(Alyr_mm2_tiling_8_clean2$query_scaffold),
                                     Colors = original_colors_mum)

chromPlot(bands = Alyr_mum_tiling_8G_col)

## Lastal

Alyr_lst_tiling <- read.delim("Alyr_last_tiling.agp", header=FALSE, comment.char="#")
colnames(Alyr_lst_tiling) <- c("final_scaffold", "start", "end", "n", "contig_or_gap", "query_scaffold", "query_scaffold_start", "query_scaffold_end", "orientation")

# Formatting the data to have the same format as hg_CytoBandIdeo data
# 5 columns with the following variables: Chrom, Start, End, Name, gieStain

# Removing gaps

Alyr_lst_tiling <- filter(Alyr_lst_tiling, query_scaffold!=100)

# Selecting only 8 chromosomes

Alyr_lst_tiling_8 <- Alyr_lst_tiling %>%  dplyr::filter(final_scaffold=="scaffold_1_RagTag" | 
                                                          final_scaffold=="scaffold_2_RagTag" | final_scaffold=="scaffold_3_RagTag" |
                                                          final_scaffold=="scaffold_4_RagTag" | final_scaffold=="scaffold_5_RagTag" |
                                                          final_scaffold=="scaffold_6_RagTag" | final_scaffold=="scaffold_7_RagTag" |
                                                          final_scaffold=="scaffold_8_RagTag")

# Splitting the scaffold column to get just numbers

Alyr_lst_tiling_8_clean <- separate(Alyr_lst_tiling_8, final_scaffold, into = c("final_scaffold", "number"), sep = "_" )
Alyr_lst_tiling_8_clean <- mutate(Alyr_lst_tiling_8_clean, scaffold = Alyr_lst_tiling_8_clean$number)
Alyr_lst_tiling_8_clean2 <- dplyr::select(Alyr_lst_tiling_8_clean, scaffold, start, end, query_scaffold)

# Now we match the colors set for minimap2 to mummer's

matching_colors_lst <- match(Alyr_lst_tiling_8_clean2$query_scaffold, Alyr_mm2_tiling_8_clean2$query_scaffold)

original_colors <- Alyr_mm2_tiling_8G_col$Colors

new_colors <- na.omit(Alyr_mm2_tiling_8G_col[matching_colors_lst, "Colors"])

# We create a new vector from the old colors and move some
# around depending on their mapping position

original_colors_lst <- original_colors

original_colors_lst[na.omit(matching_colors_lst)] <- new_colors

# Turn into grey the ones not matching

original_colors_lst[which(is.na(matching_colors_lst))] <- "lightgrey"

# Turn into grey the unmapped too

unmapped_colors_lst <- which(is.na(match(Alyr_mm2_tiling_8_clean2$query_scaffold, Alyr_lst_tiling_8_clean2$query_scaffold)))
original_colors_lst[unmapped_colors_lst] <- "lightgrey"

Alyr_lst_tiling_8G_col <- data.frame(Chrom = as.character(Alyr_mm2_tiling_8_clean2$scaffold), 
                                     Start = as.integer(Alyr_mm2_tiling_8_clean2$start),
                                     End = as.integer(Alyr_mm2_tiling_8_clean2$end),
                                     Name = as.character(Alyr_mm2_tiling_8_clean2$query_scaffold),
                                     Colors = original_colors_lst)

chromPlot(bands = Alyr_lst_tiling_8G_col)

# Now we plot synteny with circular plots

chromosomes <- as.character(c(1:16))
scaffold_length <- as.matrix(cbind(c(1:16), c(33132539, 19320864, 24464547,
                                              23328337, 21221946, 25113588, 24649197, 22951293,
                                              33132539, 19320864, 24464547, 23328337, 21221946,
                                              25113588, 24649197, 22951293)))

# We inialize the plot with chromosome names and scaffold lengths

circos.initialize(factors = chromosomes, xlim = scaffold_length)

# We show the outcome

circos.track(chromosomes, ylim = c(0, 1))

# We thin the tracks

circos.par(track.height = 0.1)

# circos.par(RESET = TRUE)

circos.initialize(factors = chromosomes, xlim = scaffold_length)

circos.track(chromosomes, ylim = c(0, 1))


# To draw links:
# bed1 = generateRandomBed(nr = 100)
# bed1 = bed1[sample(nrow(bed1), 20), ]
# bed2 = generateRandomBed(nr = 100)
# bed2 = bed2[sample(nrow(bed2), 20), ]
# circos.genomicLink(bed1, bed2, col = rand_color(nrow(bed1), transparency = 0.5),  border = NA)

bed1 <- select(Alyr_mum_tiling_8_clean2, chr = scaffold, start, end)

bed1 <- cbind(bed1, value1 = rep(0, nrow(bed1)))

bed2 <- select(Alyr_mm2_tiling_8_clean2, chr = scaffold, start, end)

bed2 <- cbind(bed2, value1 = rep(0, nrow(bed2)))

bed3 <- select(Alyr_lst_tiling_8_clean2, chr = scaffold, start, end)

bed3 <- cbind(bed3, value1 = rep(0, nrow(bed3)))

## To plot mm2 vs the other two tools, we add 8 to the scaffold # for
## coordinate matching. mum and lst will have the first 8 "bands" while
## mm2 will have the spots from 9 to 16

bed2$chr <- as.character(as.integer(bed2$chr) + 8)

## We match scaffolds to find out the links between them

match_mum_mm2 <- match(Alyr_mum_tiling_8_clean2$query_scaffold, Alyr_mm2_tiling_8_clean2$query_scaffold)

match_lst_mm2 <- match(Alyr_lst_tiling_8_clean2$query_scaffold, Alyr_mm2_tiling_8_clean2$query_scaffold)

## mm2 vs mum

bed2_filter <- bed2[na.omit(match_mum_mm2),]

bed1_filter <- bed1[!is.na(match_mum_mm2),]

link_colors <- c(rep("red", sum(bed1_filter$chr=="1")),
                 rep("green", sum(bed1_filter$chr=="2")),
                 rep("blue", sum(bed1_filter$chr=="3")),
                 rep("yellow", sum(bed1_filter$chr=="4")),
                 rep("brown", sum(bed1_filter$chr=="5")),
                 rep("purple", sum(bed1_filter$chr=="6")),
                 rep("turquoise", sum(bed1_filter$chr=="7")),
                 rep("pink", sum(bed1_filter$chr=="8")))

circos.genomicLink(bed1_filter, bed2_filter, col = link_colors,  border = NA)

circos.text(scaffold_length[1,2]/2, 0.5, "1", 1, 1, facing = "bending.inside")
circos.text(scaffold_length[2,2]/2, 0.5, "2", 2, 1, facing = "bending.inside")
circos.text(scaffold_length[3,2]/2, 0.5, "3", 3, 1, facing = "bending.inside")
circos.text(scaffold_length[4,2]/2, 0.5, "4", 4, 1, facing = "bending.inside")
circos.text(scaffold_length[5,2]/2, 0.5, "5", 5, 1, facing = "bending.inside")
circos.text(scaffold_length[6,2]/2, 0.5, "6", 6, 1, facing = "bending.inside")
circos.text(scaffold_length[7,2]/2, 0.5, "7", 7, 1, facing = "bending.inside")
circos.text(scaffold_length[8,2]/2, 0.5, "8", 8, 1, facing = "bending.inside")

circos.text(scaffold_length[1,2]/2, 0.5, "1", 9, 1, facing = "bending.inside")
circos.text(scaffold_length[2,2]/2, 0.5, "2", 10, 1, facing = "bending.inside")
circos.text(scaffold_length[3,2]/2, 0.5, "3", 11, 1, facing = "bending.inside")
circos.text(scaffold_length[4,2]/2, 0.5, "4", 12, 1, facing = "bending.inside")
circos.text(scaffold_length[5,2]/2, 0.5, "5", 13, 1, facing = "bending.inside")
circos.text(scaffold_length[6,2]/2, 0.5, "6", 14, 1, facing = "bending.inside")
circos.text(scaffold_length[7,2]/2, 0.5, "7", 15, 1, facing = "bending.inside")
circos.text(scaffold_length[8,2]/2, 0.5, "8", 16, 1, facing = "bending.inside")

circos.text(scaffold_length[4,2], 2.2, "mummer", 4, 1, facing = "bending.inside", cex = 2)
circos.text(scaffold_length[12,2], 2.2, "minimap2", 12, 1, facing = "bending.inside", cex = 2)

## mm2 vs lst

circos.initialize(factors = chromosomes, xlim = scaffold_length)

circos.track(chromosomes, ylim = c(0, 1))

bed2_filter <- bed2[na.omit(match_lst_mm2),]

bed3_filter <- bed3[!is.na(match_lst_mm2),]

link_colors <- c(rep("red", sum(bed3_filter$chr=="1")),
                 rep("green", sum(bed3_filter$chr=="2")),
                 rep("blue", sum(bed3_filter$chr=="3")),
                 rep("yellow", sum(bed3_filter$chr=="4")),
                 rep("brown", sum(bed3_filter$chr=="5")),
                 rep("purple", sum(bed3_filter$chr=="6")),
                 rep("turquoise", sum(bed3_filter$chr=="7")),
                 rep("pink", sum(bed3_filter$chr=="8")))

circos.genomicLink(bed3_filter, bed2_filter, col = link_colors,  border = NA)

circos.text(scaffold_length[1,2]/2, 0.5, "1", 1, 1, facing = "bending.inside")
circos.text(scaffold_length[2,2]/2, 0.5, "2", 2, 1, facing = "bending.inside")
circos.text(scaffold_length[3,2]/2, 0.5, "3", 3, 1, facing = "bending.inside")
circos.text(scaffold_length[4,2]/2, 0.5, "4", 4, 1, facing = "bending.inside")
circos.text(scaffold_length[5,2]/2, 0.5, "5", 5, 1, facing = "bending.inside")
circos.text(scaffold_length[6,2]/2, 0.5, "6", 6, 1, facing = "bending.inside")
circos.text(scaffold_length[7,2]/2, 0.5, "7", 7, 1, facing = "bending.inside")
circos.text(scaffold_length[8,2]/2, 0.5, "8", 8, 1, facing = "bending.inside")

circos.text(scaffold_length[1,2]/2, 0.5, "1", 9, 1, facing = "bending.inside")
circos.text(scaffold_length[2,2]/2, 0.5, "2", 10, 1, facing = "bending.inside")
circos.text(scaffold_length[3,2]/2, 0.5, "3", 11, 1, facing = "bending.inside")
circos.text(scaffold_length[4,2]/2, 0.5, "4", 12, 1, facing = "bending.inside")
circos.text(scaffold_length[5,2]/2, 0.5, "5", 13, 1, facing = "bending.inside")
circos.text(scaffold_length[6,2]/2, 0.5, "6", 14, 1, facing = "bending.inside")
circos.text(scaffold_length[7,2]/2, 0.5, "7", 15, 1, facing = "bending.inside")
circos.text(scaffold_length[8,2]/2, 0.5, "8", 16, 1, facing = "bending.inside")

circos.text(scaffold_length[4,2], 2.2, "last", 4, 1, facing = "bending.inside", cex = 2)
circos.text(scaffold_length[12,2], 2.2, "minimap2", 12, 1, facing = "bending.inside", cex = 2)

## mum vs lst

bed1$chr <- as.character(as.integer(bed1$chr) + 8)

## We match scaffolds to find out the links between them

match_lst_mum <- match(Alyr_lst_tiling_8_clean2$query_scaffold, Alyr_mum_tiling_8_clean2$query_scaffold)

circos.initialize(factors = chromosomes, xlim = scaffold_length)

circos.track(chromosomes, ylim = c(0, 1))

bed1_filter <- bed1[na.omit(match_lst_mum),]

bed3_filter <- bed3[!is.na(match_lst_mum),]

link_colors <- c(rep("red", sum(bed3_filter$chr=="1")),
                 rep("green", sum(bed3_filter$chr=="2")),
                 rep("blue", sum(bed3_filter$chr=="3")),
                 rep("yellow", sum(bed3_filter$chr=="4")),
                 rep("brown", sum(bed3_filter$chr=="5")),
                 rep("purple", sum(bed3_filter$chr=="6")),
                 rep("turquoise", sum(bed3_filter$chr=="7")),
                 rep("pink", sum(bed3_filter$chr=="8")))

circos.genomicLink(bed3_filter, bed1_filter, col = link_colors,  border = NA)

circos.text(scaffold_length[1,2]/2, 0.5, "1", 1, 1, facing = "bending.inside")
circos.text(scaffold_length[2,2]/2, 0.5, "2", 2, 1, facing = "bending.inside")
circos.text(scaffold_length[3,2]/2, 0.5, "3", 3, 1, facing = "bending.inside")
circos.text(scaffold_length[4,2]/2, 0.5, "4", 4, 1, facing = "bending.inside")
circos.text(scaffold_length[5,2]/2, 0.5, "5", 5, 1, facing = "bending.inside")
circos.text(scaffold_length[6,2]/2, 0.5, "6", 6, 1, facing = "bending.inside")
circos.text(scaffold_length[7,2]/2, 0.5, "7", 7, 1, facing = "bending.inside")
circos.text(scaffold_length[8,2]/2, 0.5, "8", 8, 1, facing = "bending.inside")

circos.text(scaffold_length[1,2]/2, 0.5, "1", 9, 1, facing = "bending.inside")
circos.text(scaffold_length[2,2]/2, 0.5, "2", 10, 1, facing = "bending.inside")
circos.text(scaffold_length[3,2]/2, 0.5, "3", 11, 1, facing = "bending.inside")
circos.text(scaffold_length[4,2]/2, 0.5, "4", 12, 1, facing = "bending.inside")
circos.text(scaffold_length[5,2]/2, 0.5, "5", 13, 1, facing = "bending.inside")
circos.text(scaffold_length[6,2]/2, 0.5, "6", 14, 1, facing = "bending.inside")
circos.text(scaffold_length[7,2]/2, 0.5, "7", 15, 1, facing = "bending.inside")
circos.text(scaffold_length[8,2]/2, 0.5, "8", 16, 1, facing = "bending.inside")

circos.text(scaffold_length[4,2], 2.2, "last", 4, 1, facing = "bending.inside", cex = 2)
circos.text(scaffold_length[12,2], 2.2, "mummer", 12, 1, facing = "bending.inside", cex = 2)



## let's plot one chromosome at a time

chromosome1 <- as.character(c(1,9))
scaffold_length1 <- as.matrix(cbind(c(1,9), c(33132539, 33132539)))

circos.initialize(factors = chromosome1, xlim = scaffold_length1)

circos.track(chromosome1, ylim = c(0, 1))

Ahal_mum_1 <- filter(Ahal_mum_tiling_8_clean2, scaffold == 1)

bed11 <- select(Ahal_mum_1, chr = scaffold, start, end)

bed11 <- cbind(bed11, value1 = rep(0, nrow(bed11)))

Ahal_mm2_1 <- filter(Ahal_mum_tiling_8_clean2, scaffold == 1)

bed21 <- select(Ahal_mm2_1, chr = scaffold, start, end)

bed21 <- cbind(bed21, value1 = rep(0, nrow(bed21)))

bed21$chr <- as.character(as.integer(bed21$chr) + 8)

match_mum_mm2_1 <- match(Ahal_mum_1$query_scaffold, Ahal_mm2_1$query_scaffold)

match_mum_mm2_1na <- is.na(match_mum_mm2_1)

bed21 <- bed21[na.omit(match_mum_mm2_1),]

bed11 <- bed11[!is.na(match_mum_mm2_1),]

circos.genomicLink(bed11, bed21, col = rand_color(nrow(bed11), transparency = 0.5),  border = NA)

