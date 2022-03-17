### Script to plot Barplot of DMRs over generations

### Import libraries

library(tidyverse)
library(data.table)
library(patchwork)

### Import first set of files, natural ALK vs progenitors (HM)

setwd("ARPEGGIO_results_alk1Vpro1_HM/dmrseq")
dmrseq_output_CG_hal <- fread("CG_context/parent1_v_allo.txt")
dmrseq_output_CHG_hal <- fread("CHG_context/parent1_v_allo.txt")
dmrseq_output_CHH_hal <- fread("CHH_context/parent1_v_allo.txt")
dmrseq_output_CG_lyr <- fread("CG_context/parent2_v_allo.txt")
dmrseq_output_CHG_lyr <- fread("CHG_context/parent2_v_allo.txt")
dmrseq_output_CHH_lyr <- fread("CHH_context/parent2_v_allo.txt")


## Filter only significant regions

dmrseq_output_CG_hal_sig <- dmrseq_output_CG_hal[dmrseq_output_CG_hal$qval<0.05]
dmrseq_output_CHG_hal_sig <- dmrseq_output_CHG_hal[dmrseq_output_CHG_hal$qval<0.05]
dmrseq_output_CHH_hal_sig <- dmrseq_output_CHH_hal[dmrseq_output_CHH_hal$qval<0.05]
dmrseq_output_CG_lyr_sig <- dmrseq_output_CG_lyr[dmrseq_output_CG_lyr$qval<0.05]
dmrseq_output_CHG_lyr_sig <- dmrseq_output_CHG_lyr[dmrseq_output_CHG_lyr$qval<0.05]
dmrseq_output_CHH_lyr_sig <- dmrseq_output_CHH_lyr[dmrseq_output_CHH_lyr$qval<0.05]

## Filter low coverage scaffolds

HM_hal_lowC_scaffolds  <- read.table("HM_RS7_G4_hal_LowCovScaffolds.txt", quote="\"", comment.char="")
HM_lyr_lowC_scaffolds  <- read.table("HM_RS7_G4_lyr_LowCovScaffolds.txt", quote="\"", comment.char="")

matchCG <- dmrseq_output_CG_hal_sig$seqnames %in% HM_hal_lowC_scaffolds$V1
dmrseq_output_CG_hal_sig <- dmrseq_output_CG_hal_sig[!matchCG,]
matchCHG <- dmrseq_output_CHG_hal_sig$seqnames %in% HM_hal_lowC_scaffolds$V1
dmrseq_output_CHG_hal_sig <- dmrseq_output_CHG_hal_sig[!matchCHG,]
matchCHH <- dmrseq_output_CHH_hal_sig$seqnames %in% HM_hal_lowC_scaffolds$V1
dmrseq_output_CHH_hal_sig <- dmrseq_output_CHH_hal_sig[!matchCHH,]

## Count regions showing increase or decrease in methylation

CG_hal <- ifelse(dmrseq_output_CG_hal_sig$stat>0, "decrease", "increase")
CHG_hal <- ifelse(dmrseq_output_CHG_hal_sig$stat>0, "decrease", "increase")
CHH_hal <- ifelse(dmrseq_output_CHH_hal_sig$stat>0, "decrease", "increase")
CG_lyr <- ifelse(dmrseq_output_CG_lyr_sig$stat>0, "decrease", "increase")
CHG_lyr <- ifelse(dmrseq_output_CHG_lyr_sig$stat>0, "decrease", "increase")
CHH_lyr <- ifelse(dmrseq_output_CHH_lyr_sig$stat>0, "decrease", "increase")

## To create the barplot, we first generate a dataframe with all
## the information we need about changes from both sides

context <- rep(c("CG", "CHG", "CHH"), 4)

species <- c(rep("halleri", 6), rep("lyrata", 6))

hypo_hyper <- c(rep("hyper", 3), rep("hypo", 3), rep("hyper", 3), rep("hypo", 3))

value <- c(table(CG_hal)[2], table(CHG_hal)[2], table(CHH_hal)[2],
           table(CG_hal)[1], table(CHG_hal)[1], table(CHH_hal)[1],
           table(CG_lyr)[2], table(CHG_lyr)[2], table(CHH_lyr)[2],
           table(CG_lyr)[1], table(CHG_lyr)[1], table(CHH_lyr)[1])

comparison <- c(rep("progenitors", 12))

dmrs_context_total <- data.frame(context, species, hypo_hyper, value, comparison)

### Now we import the second set of files and repeat the same operations
### At the end we will merge everything with the previous data frame

## Import data

setwd("ARPEGGIO_results_alk1Vsyn1_HM/dmrseq")
dmrseq_output_CG_AvB <- fread("CG_context/A_v_B_polyploid.txt")
dmrseq_output_CHG_AvB <- fread("CHG_context/A_v_B_polyploid.txt")
dmrseq_output_CHH_AvB <- fread("CHH_context/A_v_B_polyploid.txt")

## Split the files on two sides

scaffold_number_CG <- as.integer(str_split_fixed(dmrseq_output_CG_AvB$seqnames, "_", 2)[,2])
scaffold_number_CHG <- as.integer(str_split_fixed(dmrseq_output_CHG_AvB$seqnames, "_", 2)[,2])
scaffold_number_CHH <- as.integer(str_split_fixed(dmrseq_output_CHH_AvB$seqnames, "_", 2)[,2])

dmrseq_output_CG_hal <- dmrseq_output_CG_AvB[scaffold_number_CG < 2240]
dmrseq_output_CG_lyr <- dmrseq_output_CG_AvB[scaffold_number_CG > 2239]
dmrseq_output_CHG_hal <- dmrseq_output_CHG_AvB[scaffold_number_CHG < 2240]
dmrseq_output_CHG_lyr <- dmrseq_output_CHG_AvB[scaffold_number_CHG > 2240]
dmrseq_output_CHH_hal <- dmrseq_output_CHH_AvB[scaffold_number_CHH < 2240]
dmrseq_output_CHH_lyr <- dmrseq_output_CHH_AvB[scaffold_number_CHH > 2240]


## Filter only significant regions

dmrseq_output_CG_hal_sig <- dmrseq_output_CG_hal[dmrseq_output_CG_hal$qval<0.05]
dmrseq_output_CHG_hal_sig <- dmrseq_output_CHG_hal[dmrseq_output_CHG_hal$qval<0.05]
dmrseq_output_CHH_hal_sig <- dmrseq_output_CHH_hal[dmrseq_output_CHH_hal$qval<0.05]
dmrseq_output_CG_lyr_sig <- dmrseq_output_CG_lyr[dmrseq_output_CG_lyr$qval<0.05]
dmrseq_output_CHG_lyr_sig <- dmrseq_output_CHG_lyr[dmrseq_output_CHG_lyr$qval<0.05]
dmrseq_output_CHH_lyr_sig <- dmrseq_output_CHH_lyr[dmrseq_output_CHH_lyr$qval<0.05]

## Filter low coverage scaffolds

matchCG <- dmrseq_output_CG_hal_sig$seqnames %in% HM_hal_lowC_scaffolds$V1
dmrseq_output_CG_hal_sig <- dmrseq_output_CG_hal_sig[!matchCG,]
matchCHG <- dmrseq_output_CHG_hal_sig$seqnames %in% HM_hal_lowC_scaffolds$V1
dmrseq_output_CHG_hal_sig <- dmrseq_output_CHG_hal_sig[!matchCHG,]
matchCHH <- dmrseq_output_CHH_hal_sig$seqnames %in% HM_hal_lowC_scaffolds$V1
dmrseq_output_CHH_hal_sig <- dmrseq_output_CHH_hal_sig[!matchCHH,]

## Count regions showing increase or decrease in methylation

CG_hal <- ifelse(dmrseq_output_CG_hal_sig$stat>0, "decrease", "increase")
CHG_hal <- ifelse(dmrseq_output_CHG_hal_sig$stat>0, "decrease", "increase")
CHH_hal <- ifelse(dmrseq_output_CHH_hal_sig$stat>0, "decrease", "increase")
CG_lyr <- ifelse(dmrseq_output_CG_lyr_sig$stat>0, "decrease", "increase")
CHG_lyr <- ifelse(dmrseq_output_CHG_lyr_sig$stat>0, "decrease", "increase")
CHH_lyr <- ifelse(dmrseq_output_CHH_lyr_sig$stat>0, "decrease", "increase")

## To create the barplot, we first generate a dataframe with all
## the information we need about changes from both sides

context <- rep(c("CG", "CHG", "CHH"), 4)

species <- c(rep("halleri", 6), rep("lyrata", 6))

hypo_hyper <- c(rep("hyper", 3), rep("hypo", 3), rep("hyper", 3), rep("hypo", 3))

value <- c(table(CG_hal)[2], table(CHG_hal)[2], table(CHH_hal)[2],
           table(CG_hal)[1], table(CHG_hal)[1], table(CHH_hal)[1],
           table(CG_lyr)[2], table(CHG_lyr)[2], table(CHH_lyr)[2],
           table(CG_lyr)[1], table(CHG_lyr)[1], table(CHH_lyr)[1])

comparison <- c(rep("synthetic_1", 12))

dmrs_context_total_new <- data.frame(context, species, hypo_hyper, value, comparison)

dmrs_context_total <- rbind(dmrs_context_total, dmrs_context_total_new)

### We import the data one last time
### natural vs synthetic generation 4

setwd("ARPEGGIO_results_alk1Vsyn4_HM/dmrseq")
dmrseq_output_CG_AvB <- fread("CG_context/A_v_B_polyploid.txt")
dmrseq_output_CHG_AvB <- fread("CHG_context/A_v_B_polyploid.txt")
dmrseq_output_CHH_AvB <- fread("CHH_context/A_v_B_polyploid.txt")

## Split the files on two sides

scaffold_number_CG <- as.integer(str_split_fixed(dmrseq_output_CG_AvB$seqnames, "_", 2)[,2])
scaffold_number_CHG <- as.integer(str_split_fixed(dmrseq_output_CHG_AvB$seqnames, "_", 2)[,2])
scaffold_number_CHH <- as.integer(str_split_fixed(dmrseq_output_CHH_AvB$seqnames, "_", 2)[,2])

dmrseq_output_CG_hal <- dmrseq_output_CG_AvB[scaffold_number_CG < 2240]
dmrseq_output_CG_lyr <- dmrseq_output_CG_AvB[scaffold_number_CG > 2239]
dmrseq_output_CHG_hal <- dmrseq_output_CHG_AvB[scaffold_number_CHG < 2240]
dmrseq_output_CHG_lyr <- dmrseq_output_CHG_AvB[scaffold_number_CHG > 2240]
dmrseq_output_CHH_hal <- dmrseq_output_CHH_AvB[scaffold_number_CHH < 2240]
dmrseq_output_CHH_lyr <- dmrseq_output_CHH_AvB[scaffold_number_CHH > 2240]


## Filter only significant regions

dmrseq_output_CG_hal_sig <- dmrseq_output_CG_hal[dmrseq_output_CG_hal$qval<0.05]
dmrseq_output_CHG_hal_sig <- dmrseq_output_CHG_hal[dmrseq_output_CHG_hal$qval<0.05]
dmrseq_output_CHH_hal_sig <- dmrseq_output_CHH_hal[dmrseq_output_CHH_hal$qval<0.05]
dmrseq_output_CG_lyr_sig <- dmrseq_output_CG_lyr[dmrseq_output_CG_lyr$qval<0.05]
dmrseq_output_CHG_lyr_sig <- dmrseq_output_CHG_lyr[dmrseq_output_CHG_lyr$qval<0.05]
dmrseq_output_CHH_lyr_sig <- dmrseq_output_CHH_lyr[dmrseq_output_CHH_lyr$qval<0.05]

## Filter low coverage scaffolds

matchCG <- dmrseq_output_CG_hal_sig$seqnames %in% HM_hal_lowC_scaffolds$V1
dmrseq_output_CG_hal_sig <- dmrseq_output_CG_hal_sig[!matchCG,]
matchCHG <- dmrseq_output_CHG_hal_sig$seqnames %in% HM_hal_lowC_scaffolds$V1
dmrseq_output_CHG_hal_sig <- dmrseq_output_CHG_hal_sig[!matchCHG,]
matchCHH <- dmrseq_output_CHH_hal_sig$seqnames %in% HM_hal_lowC_scaffolds$V1
dmrseq_output_CHH_hal_sig <- dmrseq_output_CHH_hal_sig[!matchCHH,]

## Count regions showing increase or decrease in methylation

CG_hal <- ifelse(dmrseq_output_CG_hal_sig$stat>0, "decrease", "increase")
CHG_hal <- ifelse(dmrseq_output_CHG_hal_sig$stat>0, "decrease", "increase")
CHH_hal <- ifelse(dmrseq_output_CHH_hal_sig$stat>0, "decrease", "increase")
CG_lyr <- ifelse(dmrseq_output_CG_lyr_sig$stat>0, "decrease", "increase")
CHG_lyr <- ifelse(dmrseq_output_CHG_lyr_sig$stat>0, "decrease", "increase")
CHH_lyr <- ifelse(dmrseq_output_CHH_lyr_sig$stat>0, "decrease", "increase")

## To create the barplot, we first generate a dataframe with all
## the information we need about changes from both sides

context <- rep(c("CG", "CHG", "CHH"), 4)

species <- c(rep("halleri", 6), rep("lyrata", 6))

hypo_hyper <- c(rep("hyper", 3), rep("hypo", 3), rep("hyper", 3), rep("hypo", 3))

value <- c(table(CG_hal)[2], table(CHG_hal)[2], table(CHH_hal)[2],
           table(CG_hal)[1], table(CHG_hal)[1], table(CHH_hal)[1],
           table(CG_lyr)[2], table(CHG_lyr)[2], table(CHH_lyr)[2],
           table(CG_lyr)[1], table(CHG_lyr)[1], table(CHH_lyr)[1])

comparison <- c(rep("synthetic_4", 12))

dmrs_context_total_new <- data.frame(context, species, hypo_hyper, value, comparison)

dmrs_context_total <- rbind(dmrs_context_total, dmrs_context_total_new)

## To check the dataframe used with ggplot
## View(dmrs_context_total)
## Barplot code

## We create a vector to change the labels on top of the plot
## i.e. from halleri to halleri-side

species.labs <- c(paste0("halleri-side"),
                  paste0("lyrata-side"))

names(species.labs) <- c("halleri", "lyrata")

## Second plot showing the proportion of hyper and hypo methylated regions

cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

ggplot(dmrs_context_total, aes(x=factor(interaction(context, comparison), levels = c("CG.progenitors",
                                                                                     "CG.synthetic_1",
                                                                                     "CG.synthetic_4",
                                                                                     "CHG.progenitors",
                                                                                     "CHG.synthetic_1",
                                                                                     "CHG.synthetic_4",
                                                                                     "CHH.progenitors",
                                                                                     "CHH.synthetic_1",
                                                                                     "CHH.synthetic_4")), y=value, fill = comparison, col = factor(hypo_hyper, levels = c("hyper", "hypo")))) +
  geom_col(width = 0.5, show.legend = TRUE, position = "stack") +
  scale_fill_manual(values=cbPalette) +
  scale_color_manual(labels = c("hypo", "hyper"), values=c("white", "black")) +
  facet_grid(~species, labeller = labeller(species = species.labs)) +
  ylab("Number of DMRs") +
  theme_bw() +
  theme(legend.title = element_blank(), text = element_text(size=25),
        legend.text = element_text(size = 15),
        strip.background = element_rect(color="white", fill="white", size=1.5, linetype="solid"),
        axis.ticks = element_blank()) +
  scale_x_discrete(name = "DMR context", breaks = ,
                   labels = c("", "CG", "", "", "CHG", "", "", "CHH", "")) +
  ylim(0, 35000) +
  geom_vline(xintercept = c(3.5, 6.5), linetype = "longdash") +
  ggtitle("ALK (Cold Conditions)")

# Simple line plot for trend over comparisons

# We modify the original dataframe so that comparisons are integers instead of strings

dmrs_context_total_new <- data.frame(context = rep(c("CG", "CHG", "CHH"), 6),
           side = c(rep(c(rep("halleri", 3),
                    rep("lyrata", 3)), 3)),
           value = c(dmrs_context_total$value[1] + dmrs_context_total$value[4],
                     dmrs_context_total$value[2] + dmrs_context_total$value[5],
                     dmrs_context_total$value[3] + dmrs_context_total$value[6],
                     dmrs_context_total$value[7] + dmrs_context_total$value[10],
                     dmrs_context_total$value[8] + dmrs_context_total$value[11],
                     dmrs_context_total$value[9] + dmrs_context_total$value[12],
                     dmrs_context_total$value[13] + dmrs_context_total$value[16],
                     dmrs_context_total$value[14] + dmrs_context_total$value[17],
                     dmrs_context_total$value[15] + dmrs_context_total$value[18],
                     dmrs_context_total$value[19] + dmrs_context_total$value[22],
                     dmrs_context_total$value[20] + dmrs_context_total$value[23],
                     dmrs_context_total$value[21] + dmrs_context_total$value[24],
                     dmrs_context_total$value[25] + dmrs_context_total$value[28],
                     dmrs_context_total$value[26] + dmrs_context_total$value[29],
                     dmrs_context_total$value[27] + dmrs_context_total$value[30],
                     dmrs_context_total$value[31] + dmrs_context_total$value[34],
                     dmrs_context_total$value[32] + dmrs_context_total$value[35],
                     dmrs_context_total$value[33] + dmrs_context_total$value[36]),
           comparison = c(rep("NvsP1", 6), rep("NvsS1", 6), rep("NvsS4",6)))


HM_ALK <- ggplot(dmrs_context_total_new, aes(x=comparison, y=value, group = context, colour = context)) +
  geom_line(show.legend = TRUE, size = 2) +
  #scale_x_continuous(labels = c("pro", "G1", "G2", "G3", "G4")) +
  #scale_fill_manual(values=cbPalette) +
  #scale_color_manual(labels = c("hypo", "hyper"), values=c("white", "black")) +
  facet_grid(~side) +
  ylab("Number of DMRs") +
  xlab("Generation") +
  theme_bw() +
  theme(legend.title = element_blank(), text = element_text(size=15),
        legend.text = element_text(size = 10),
        #strip.background = element_rect(color="white", fill="white", size=1.5, linetype="solid"),
        axis.ticks = element_blank()) +
  ggtitle("ALK (Cold Conditions)") +
  ylim(0, 35000)

#####################################################################
# Repeat the whole procedure from above to compare TKS vs all (HM) ##
#####################################################################

### Import first set of files, natural TKS vs progenitors (HM)

setwd("ARPEGGIO_results_tks1Vpro1_HM/dmrseq")
dmrseq_output_CG_hal <- fread("CG_context/parent1_v_allo.txt")
dmrseq_output_CHG_hal <- fread("CHG_context/parent1_v_allo.txt")
dmrseq_output_CHH_hal <- fread("CHH_context/parent1_v_allo.txt")
dmrseq_output_CG_lyr <- fread("CG_context/parent2_v_allo.txt")
dmrseq_output_CHG_lyr <- fread("CHG_context/parent2_v_allo.txt")
dmrseq_output_CHH_lyr <- fread("CHH_context/parent2_v_allo.txt")


## Filter only significant regions

dmrseq_output_CG_hal_sig <- dmrseq_output_CG_hal[dmrseq_output_CG_hal$qval<0.05]
dmrseq_output_CHG_hal_sig <- dmrseq_output_CHG_hal[dmrseq_output_CHG_hal$qval<0.05]
dmrseq_output_CHH_hal_sig <- dmrseq_output_CHH_hal[dmrseq_output_CHH_hal$qval<0.05]
dmrseq_output_CG_lyr_sig <- dmrseq_output_CG_lyr[dmrseq_output_CG_lyr$qval<0.05]
dmrseq_output_CHG_lyr_sig <- dmrseq_output_CHG_lyr[dmrseq_output_CHG_lyr$qval<0.05]
dmrseq_output_CHH_lyr_sig <- dmrseq_output_CHH_lyr[dmrseq_output_CHH_lyr$qval<0.05]

## Filter low coverage scaffolds

HM_hal_lowC_scaffolds  <- read.table("HM_RS7_G4_hal_LowCovScaffolds.txt", quote="\"", comment.char="")
HM_lyr_lowC_scaffolds  <- read.table("HM_RS7_G4_lyr_LowCovScaffolds.txt", quote="\"", comment.char="")

matchCG <- dmrseq_output_CG_hal_sig$seqnames %in% HM_hal_lowC_scaffolds$V1
dmrseq_output_CG_hal_sig <- dmrseq_output_CG_hal_sig[!matchCG,]
matchCHG <- dmrseq_output_CHG_hal_sig$seqnames %in% HM_hal_lowC_scaffolds$V1
dmrseq_output_CHG_hal_sig <- dmrseq_output_CHG_hal_sig[!matchCHG,]
matchCHH <- dmrseq_output_CHH_hal_sig$seqnames %in% HM_hal_lowC_scaffolds$V1
dmrseq_output_CHH_hal_sig <- dmrseq_output_CHH_hal_sig[!matchCHH,]

## Count regions showing increase or decrease in methylation

CG_hal <- ifelse(dmrseq_output_CG_hal_sig$stat>0, "decrease", "increase")
CHG_hal <- ifelse(dmrseq_output_CHG_hal_sig$stat>0, "decrease", "increase")
CHH_hal <- ifelse(dmrseq_output_CHH_hal_sig$stat>0, "decrease", "increase")
CG_lyr <- ifelse(dmrseq_output_CG_lyr_sig$stat>0, "decrease", "increase")
CHG_lyr <- ifelse(dmrseq_output_CHG_lyr_sig$stat>0, "decrease", "increase")
CHH_lyr <- ifelse(dmrseq_output_CHH_lyr_sig$stat>0, "decrease", "increase")

## To create the barplot, we first generate a dataframe with all
## the information we need about changes from both sides

context <- rep(c("CG", "CHG", "CHH"), 4)

species <- c(rep("halleri", 6), rep("lyrata", 6))

hypo_hyper <- c(rep("hyper", 3), rep("hypo", 3), rep("hyper", 3), rep("hypo", 3))

value <- c(table(CG_hal)[2], table(CHG_hal)[2], table(CHH_hal)[2],
           table(CG_hal)[1], table(CHG_hal)[1], table(CHH_hal)[1],
           table(CG_lyr)[2], table(CHG_lyr)[2], table(CHH_lyr)[2],
           table(CG_lyr)[1], table(CHG_lyr)[1], table(CHH_lyr)[1])

comparison <- c(rep("progenitors", 12))

dmrs_context_total <- data.frame(context, species, hypo_hyper, value, comparison)

### Now we import the second set of files and repeat the same operations
### At the end we will merge everything with the previous data frame

## Import data

setwd("ARPEGGIO_results_tks1Vsyn1_HM/dmrseq")
dmrseq_output_CG_AvB <- fread("CG_context/A_v_B_polyploid.txt")
dmrseq_output_CHG_AvB <- fread("CHG_context/A_v_B_polyploid.txt")
dmrseq_output_CHH_AvB <- fread("CHH_context/A_v_B_polyploid.txt")

## Split the files on two sides

scaffold_number_CG <- as.integer(str_split_fixed(dmrseq_output_CG_AvB$seqnames, "_", 2)[,2])
scaffold_number_CHG <- as.integer(str_split_fixed(dmrseq_output_CHG_AvB$seqnames, "_", 2)[,2])
scaffold_number_CHH <- as.integer(str_split_fixed(dmrseq_output_CHH_AvB$seqnames, "_", 2)[,2])

dmrseq_output_CG_hal <- dmrseq_output_CG_AvB[scaffold_number_CG < 2240]
dmrseq_output_CG_lyr <- dmrseq_output_CG_AvB[scaffold_number_CG > 2239]
dmrseq_output_CHG_hal <- dmrseq_output_CHG_AvB[scaffold_number_CHG < 2240]
dmrseq_output_CHG_lyr <- dmrseq_output_CHG_AvB[scaffold_number_CHG > 2240]
dmrseq_output_CHH_hal <- dmrseq_output_CHH_AvB[scaffold_number_CHH < 2240]
dmrseq_output_CHH_lyr <- dmrseq_output_CHH_AvB[scaffold_number_CHH > 2240]


## Filter only significant regions

dmrseq_output_CG_hal_sig <- dmrseq_output_CG_hal[dmrseq_output_CG_hal$qval<0.05]
dmrseq_output_CHG_hal_sig <- dmrseq_output_CHG_hal[dmrseq_output_CHG_hal$qval<0.05]
dmrseq_output_CHH_hal_sig <- dmrseq_output_CHH_hal[dmrseq_output_CHH_hal$qval<0.05]
dmrseq_output_CG_lyr_sig <- dmrseq_output_CG_lyr[dmrseq_output_CG_lyr$qval<0.05]
dmrseq_output_CHG_lyr_sig <- dmrseq_output_CHG_lyr[dmrseq_output_CHG_lyr$qval<0.05]
dmrseq_output_CHH_lyr_sig <- dmrseq_output_CHH_lyr[dmrseq_output_CHH_lyr$qval<0.05]

## Filter low coverage scaffolds

matchCG <- dmrseq_output_CG_hal_sig$seqnames %in% HM_hal_lowC_scaffolds$V1
dmrseq_output_CG_hal_sig <- dmrseq_output_CG_hal_sig[!matchCG,]
matchCHG <- dmrseq_output_CHG_hal_sig$seqnames %in% HM_hal_lowC_scaffolds$V1
dmrseq_output_CHG_hal_sig <- dmrseq_output_CHG_hal_sig[!matchCHG,]
matchCHH <- dmrseq_output_CHH_hal_sig$seqnames %in% HM_hal_lowC_scaffolds$V1
dmrseq_output_CHH_hal_sig <- dmrseq_output_CHH_hal_sig[!matchCHH,]

## Count regions showing increase or decrease in methylation

CG_hal <- ifelse(dmrseq_output_CG_hal_sig$stat>0, "decrease", "increase")
CHG_hal <- ifelse(dmrseq_output_CHG_hal_sig$stat>0, "decrease", "increase")
CHH_hal <- ifelse(dmrseq_output_CHH_hal_sig$stat>0, "decrease", "increase")
CG_lyr <- ifelse(dmrseq_output_CG_lyr_sig$stat>0, "decrease", "increase")
CHG_lyr <- ifelse(dmrseq_output_CHG_lyr_sig$stat>0, "decrease", "increase")
CHH_lyr <- ifelse(dmrseq_output_CHH_lyr_sig$stat>0, "decrease", "increase")

## To create the barplot, we first generate a dataframe with all
## the information we need about changes from both sides

context <- rep(c("CG", "CHG", "CHH"), 4)

species <- c(rep("halleri", 6), rep("lyrata", 6))

hypo_hyper <- c(rep("hyper", 3), rep("hypo", 3), rep("hyper", 3), rep("hypo", 3))

value <- c(table(CG_hal)[2], table(CHG_hal)[2], table(CHH_hal)[2],
           table(CG_hal)[1], table(CHG_hal)[1], table(CHH_hal)[1],
           table(CG_lyr)[2], table(CHG_lyr)[2], table(CHH_lyr)[2],
           table(CG_lyr)[1], table(CHG_lyr)[1], table(CHH_lyr)[1])

comparison <- c(rep("synthetic_1", 12))

dmrs_context_total_new <- data.frame(context, species, hypo_hyper, value, comparison)

dmrs_context_total <- rbind(dmrs_context_total, dmrs_context_total_new)

### We import the data one last time
### natural vs synthetic generation 4

setwd("ARPEGGIO_results_tks1Vsyn4_HM/dmrseq")
dmrseq_output_CG_AvB <- fread("CG_context/A_v_B_polyploid.txt")
dmrseq_output_CHG_AvB <- fread("CHG_context/A_v_B_polyploid.txt")
dmrseq_output_CHH_AvB <- fread("CHH_context/A_v_B_polyploid.txt")

## Split the files on two sides

scaffold_number_CG <- as.integer(str_split_fixed(dmrseq_output_CG_AvB$seqnames, "_", 2)[,2])
scaffold_number_CHG <- as.integer(str_split_fixed(dmrseq_output_CHG_AvB$seqnames, "_", 2)[,2])
scaffold_number_CHH <- as.integer(str_split_fixed(dmrseq_output_CHH_AvB$seqnames, "_", 2)[,2])

dmrseq_output_CG_hal <- dmrseq_output_CG_AvB[scaffold_number_CG < 2240]
dmrseq_output_CG_lyr <- dmrseq_output_CG_AvB[scaffold_number_CG > 2239]
dmrseq_output_CHG_hal <- dmrseq_output_CHG_AvB[scaffold_number_CHG < 2240]
dmrseq_output_CHG_lyr <- dmrseq_output_CHG_AvB[scaffold_number_CHG > 2240]
dmrseq_output_CHH_hal <- dmrseq_output_CHH_AvB[scaffold_number_CHH < 2240]
dmrseq_output_CHH_lyr <- dmrseq_output_CHH_AvB[scaffold_number_CHH > 2240]


## Filter only significant regions

dmrseq_output_CG_hal_sig <- dmrseq_output_CG_hal[dmrseq_output_CG_hal$qval<0.05]
dmrseq_output_CHG_hal_sig <- dmrseq_output_CHG_hal[dmrseq_output_CHG_hal$qval<0.05]
dmrseq_output_CHH_hal_sig <- dmrseq_output_CHH_hal[dmrseq_output_CHH_hal$qval<0.05]
dmrseq_output_CG_lyr_sig <- dmrseq_output_CG_lyr[dmrseq_output_CG_lyr$qval<0.05]
dmrseq_output_CHG_lyr_sig <- dmrseq_output_CHG_lyr[dmrseq_output_CHG_lyr$qval<0.05]
dmrseq_output_CHH_lyr_sig <- dmrseq_output_CHH_lyr[dmrseq_output_CHH_lyr$qval<0.05]

## Filter low coverage scaffolds

matchCG <- dmrseq_output_CG_hal_sig$seqnames %in% HM_hal_lowC_scaffolds$V1
dmrseq_output_CG_hal_sig <- dmrseq_output_CG_hal_sig[!matchCG,]
matchCHG <- dmrseq_output_CHG_hal_sig$seqnames %in% HM_hal_lowC_scaffolds$V1
dmrseq_output_CHG_hal_sig <- dmrseq_output_CHG_hal_sig[!matchCHG,]
matchCHH <- dmrseq_output_CHH_hal_sig$seqnames %in% HM_hal_lowC_scaffolds$V1
dmrseq_output_CHH_hal_sig <- dmrseq_output_CHH_hal_sig[!matchCHH,]

## Count regions showing increase or decrease in methylation

CG_hal <- ifelse(dmrseq_output_CG_hal_sig$stat>0, "decrease", "increase")
CHG_hal <- ifelse(dmrseq_output_CHG_hal_sig$stat>0, "decrease", "increase")
CHH_hal <- ifelse(dmrseq_output_CHH_hal_sig$stat>0, "decrease", "increase")
CG_lyr <- ifelse(dmrseq_output_CG_lyr_sig$stat>0, "decrease", "increase")
CHG_lyr <- ifelse(dmrseq_output_CHG_lyr_sig$stat>0, "decrease", "increase")
CHH_lyr <- ifelse(dmrseq_output_CHH_lyr_sig$stat>0, "decrease", "increase")

## To create the barplot, we first generate a dataframe with all
## the information we need about changes from both sides

context <- rep(c("CG", "CHG", "CHH"), 4)

species <- c(rep("halleri", 6), rep("lyrata", 6))

hypo_hyper <- c(rep("hyper", 3), rep("hypo", 3), rep("hyper", 3), rep("hypo", 3))

value <- c(table(CG_hal)[2], table(CHG_hal)[2], table(CHH_hal)[2],
           table(CG_hal)[1], table(CHG_hal)[1], table(CHH_hal)[1],
           table(CG_lyr)[2], table(CHG_lyr)[2], table(CHH_lyr)[2],
           table(CG_lyr)[1], table(CHG_lyr)[1], table(CHH_lyr)[1])

comparison <- c(rep("synthetic_4", 12))

dmrs_context_total_new <- data.frame(context, species, hypo_hyper, value, comparison)

dmrs_context_total <- rbind(dmrs_context_total, dmrs_context_total_new)

## To check the dataframe used with ggplot
## View(dmrs_context_total)
## Barplot code

## We create a vector to change the labels on top of the plot
## i.e. from halleri to halleri-side

species.labs <- c(paste0("halleri-side"),
                  paste0("lyrata-side"))

names(species.labs) <- c("halleri", "lyrata")

## Second plot showing the proportion of hyper and hypo methylated regions

cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

ggplot(dmrs_context_total, aes(x=factor(interaction(context, comparison), levels = c("CG.progenitors",
                                                                                     "CG.synthetic_1",
                                                                                     "CG.synthetic_4",
                                                                                     "CHG.progenitors",
                                                                                     "CHG.synthetic_1",
                                                                                     "CHG.synthetic_4",
                                                                                     "CHH.progenitors",
                                                                                     "CHH.synthetic_1",
                                                                                     "CHH.synthetic_4")), y=value, fill = comparison, col = factor(hypo_hyper, levels = c("hyper", "hypo")))) +
  geom_col(width = 0.5, show.legend = TRUE, position = "stack") +
  scale_fill_manual(values=cbPalette) +
  scale_color_manual(labels = c("hypo", "hyper"), values=c("white", "black")) +
  facet_grid(~species, labeller = labeller(species = species.labs)) +
  ylab("Number of DMRs") +
  theme_bw() +
  theme(legend.title = element_blank(), text = element_text(size=25),
        legend.text = element_text(size = 15),
        strip.background = element_rect(color="white", fill="white", size=1.5, linetype="solid"),
        axis.ticks = element_blank()) +
  scale_x_discrete(name = "DMR context", breaks = ,
                   labels = c("", "CG", "", "", "CHG", "", "", "CHH", "")) +
  ylim(0, 35000) +
  geom_vline(xintercept = c(3.5, 6.5), linetype = "longdash") +
  ggtitle("TKS (Cold Conditions)")

# Simple line plot for trend over comparisons

# We modify the original dataframe so that comparisons are integers instead of strings

dmrs_context_total_new <- data.frame(context = rep(c("CG", "CHG", "CHH"), 6),
                                     side = c(rep(c(rep("halleri", 3),
                                                    rep("lyrata", 3)), 3)),
                                     value = c(dmrs_context_total$value[1] + dmrs_context_total$value[4],
                                               dmrs_context_total$value[2] + dmrs_context_total$value[5],
                                               dmrs_context_total$value[3] + dmrs_context_total$value[6],
                                               dmrs_context_total$value[7] + dmrs_context_total$value[10],
                                               dmrs_context_total$value[8] + dmrs_context_total$value[11],
                                               dmrs_context_total$value[9] + dmrs_context_total$value[12],
                                               dmrs_context_total$value[13] + dmrs_context_total$value[16],
                                               dmrs_context_total$value[14] + dmrs_context_total$value[17],
                                               dmrs_context_total$value[15] + dmrs_context_total$value[18],
                                               dmrs_context_total$value[19] + dmrs_context_total$value[22],
                                               dmrs_context_total$value[20] + dmrs_context_total$value[23],
                                               dmrs_context_total$value[21] + dmrs_context_total$value[24],
                                               dmrs_context_total$value[25] + dmrs_context_total$value[28],
                                               dmrs_context_total$value[26] + dmrs_context_total$value[29],
                                               dmrs_context_total$value[27] + dmrs_context_total$value[30],
                                               dmrs_context_total$value[31] + dmrs_context_total$value[34],
                                               dmrs_context_total$value[32] + dmrs_context_total$value[35],
                                               dmrs_context_total$value[33] + dmrs_context_total$value[36]),
                                     comparison = c(rep("NvsP1", 6), rep("NvsS1", 6), rep("NvsS4",6)))


HM_TKS <- ggplot(dmrs_context_total_new, aes(x=comparison, y=value, group = context, colour = context)) +
  geom_line(show.legend = TRUE, size = 2) +
  #scale_x_continuous(labels = c("pro", "G1", "G2", "G3", "G4")) +
  #scale_fill_manual(values=cbPalette) +
  #scale_color_manual(labels = c("hypo", "hyper"), values=c("white", "black")) +
  facet_grid(~side) +
  ylab("Number of DMRs") +
  xlab("Generation") +
  theme_bw() +
  theme(legend.title = element_blank(), text = element_text(size=15),
        legend.text = element_text(size = 10),
        #strip.background = element_rect(color="white", fill="white", size=1.5, linetype="solid"),
        axis.ticks = element_blank()) +
  ggtitle("TKS (Cold Conditions)") +
  ylim(0, 35000)


#####################################################################
# Repeat the whole procedure from above to compare ALK vs all (LL) ##
#####################################################################

### Import first set of files, natural ALK vs progenitors (LL)

setwd("ARPEGGIO_results_alk1Vpro1_LL/dmrseq")
dmrseq_output_CG_hal <- fread("CG_context/parent1_v_allo.txt")
dmrseq_output_CHG_hal <- fread("CHG_context/parent1_v_allo.txt")
dmrseq_output_CHH_hal <- fread("CHH_context/parent1_v_allo.txt")
dmrseq_output_CG_lyr <- fread("CG_context/parent2_v_allo.txt")
dmrseq_output_CHG_lyr <- fread("CHG_context/parent2_v_allo.txt")
dmrseq_output_CHH_lyr <- fread("CHH_context/parent2_v_allo.txt")


## Filter only significant regions

dmrseq_output_CG_hal_sig <- dmrseq_output_CG_hal[dmrseq_output_CG_hal$qval<0.05]
dmrseq_output_CHG_hal_sig <- dmrseq_output_CHG_hal[dmrseq_output_CHG_hal$qval<0.05]
dmrseq_output_CHH_hal_sig <- dmrseq_output_CHH_hal[dmrseq_output_CHH_hal$qval<0.05]
dmrseq_output_CG_lyr_sig <- dmrseq_output_CG_lyr[dmrseq_output_CG_lyr$qval<0.05]
dmrseq_output_CHG_lyr_sig <- dmrseq_output_CHG_lyr[dmrseq_output_CHG_lyr$qval<0.05]
dmrseq_output_CHH_lyr_sig <- dmrseq_output_CHH_lyr[dmrseq_output_CHH_lyr$qval<0.05]

## Filter low coverage scaffolds

LL_hal_lowC_scaffolds  <- read.table("LL_RS7_G4_hal_LowCovScaffolds.txt", quote="\"", comment.char="")
LL_lyr_lowC_scaffolds  <- read.table("LL_RS7_G4_lyr_LowCovScaffolds.txt", quote="\"", comment.char="")

matchCG <- dmrseq_output_CG_hal_sig$seqnames %in% LL_hal_lowC_scaffolds$V1
dmrseq_output_CG_hal_sig <- dmrseq_output_CG_hal_sig[!matchCG,]
matchCHG <- dmrseq_output_CHG_hal_sig$seqnames %in% LL_hal_lowC_scaffolds$V1
dmrseq_output_CHG_hal_sig <- dmrseq_output_CHG_hal_sig[!matchCHG,]
matchCHH <- dmrseq_output_CHH_hal_sig$seqnames %in% LL_hal_lowC_scaffolds$V1
dmrseq_output_CHH_hal_sig <- dmrseq_output_CHH_hal_sig[!matchCHH,]

## Count regions showing increase or decrease in methylation

CG_hal <- ifelse(dmrseq_output_CG_hal_sig$stat>0, "decrease", "increase")
CHG_hal <- ifelse(dmrseq_output_CHG_hal_sig$stat>0, "decrease", "increase")
CHH_hal <- ifelse(dmrseq_output_CHH_hal_sig$stat>0, "decrease", "increase")
CG_lyr <- ifelse(dmrseq_output_CG_lyr_sig$stat>0, "decrease", "increase")
CHG_lyr <- ifelse(dmrseq_output_CHG_lyr_sig$stat>0, "decrease", "increase")
CHH_lyr <- ifelse(dmrseq_output_CHH_lyr_sig$stat>0, "decrease", "increase")

## To create the barplot, we first generate a dataframe with all
## the information we need about changes from both sides

context <- rep(c("CG", "CHG", "CHH"), 4)

species <- c(rep("halleri", 6), rep("lyrata", 6))

hypo_hyper <- c(rep("hyper", 3), rep("hypo", 3), rep("hyper", 3), rep("hypo", 3))

value <- c(table(CG_hal)[2], table(CHG_hal)[2], table(CHH_hal)[2],
           table(CG_hal)[1], table(CHG_hal)[1], table(CHH_hal)[1],
           table(CG_lyr)[2], table(CHG_lyr)[2], table(CHH_lyr)[2],
           table(CG_lyr)[1], table(CHG_lyr)[1], table(CHH_lyr)[1])

comparison <- c(rep("progenitors", 12))

dmrs_context_total <- data.frame(context, species, hypo_hyper, value, comparison)

### Now we import the second set of files and repeat the same operations
### At the end we will merge everything with the previous data frame

## Import data

setwd("ARPEGGIO_results_alk1Vsyn1_LL/dmrseq")
dmrseq_output_CG_AvB <- fread("CG_context/A_v_B_polyploid.txt")
dmrseq_output_CHG_AvB <- fread("CHG_context/A_v_B_polyploid.txt")
dmrseq_output_CHH_AvB <- fread("CHH_context/A_v_B_polyploid.txt")

## Split the files on two sides

scaffold_number_CG <- as.integer(str_split_fixed(dmrseq_output_CG_AvB$seqnames, "_", 2)[,2])
scaffold_number_CHG <- as.integer(str_split_fixed(dmrseq_output_CHG_AvB$seqnames, "_", 2)[,2])
scaffold_number_CHH <- as.integer(str_split_fixed(dmrseq_output_CHH_AvB$seqnames, "_", 2)[,2])

dmrseq_output_CG_hal <- dmrseq_output_CG_AvB[scaffold_number_CG < 2240]
dmrseq_output_CG_lyr <- dmrseq_output_CG_AvB[scaffold_number_CG > 2239]
dmrseq_output_CHG_hal <- dmrseq_output_CHG_AvB[scaffold_number_CHG < 2240]
dmrseq_output_CHG_lyr <- dmrseq_output_CHG_AvB[scaffold_number_CHG > 2240]
dmrseq_output_CHH_hal <- dmrseq_output_CHH_AvB[scaffold_number_CHH < 2240]
dmrseq_output_CHH_lyr <- dmrseq_output_CHH_AvB[scaffold_number_CHH > 2240]


## Filter only significant regions

dmrseq_output_CG_hal_sig <- dmrseq_output_CG_hal[dmrseq_output_CG_hal$qval<0.05]
dmrseq_output_CHG_hal_sig <- dmrseq_output_CHG_hal[dmrseq_output_CHG_hal$qval<0.05]
dmrseq_output_CHH_hal_sig <- dmrseq_output_CHH_hal[dmrseq_output_CHH_hal$qval<0.05]
dmrseq_output_CG_lyr_sig <- dmrseq_output_CG_lyr[dmrseq_output_CG_lyr$qval<0.05]
dmrseq_output_CHG_lyr_sig <- dmrseq_output_CHG_lyr[dmrseq_output_CHG_lyr$qval<0.05]
dmrseq_output_CHH_lyr_sig <- dmrseq_output_CHH_lyr[dmrseq_output_CHH_lyr$qval<0.05]

## Filter low coverage scaffolds

matchCG <- dmrseq_output_CG_hal_sig$seqnames %in% LL_hal_lowC_scaffolds$V1
dmrseq_output_CG_hal_sig <- dmrseq_output_CG_hal_sig[!matchCG,]
matchCHG <- dmrseq_output_CHG_hal_sig$seqnames %in% LL_hal_lowC_scaffolds$V1
dmrseq_output_CHG_hal_sig <- dmrseq_output_CHG_hal_sig[!matchCHG,]
matchCHH <- dmrseq_output_CHH_hal_sig$seqnames %in% LL_hal_lowC_scaffolds$V1
dmrseq_output_CHH_hal_sig <- dmrseq_output_CHH_hal_sig[!matchCHH,]

## Count regions showing increase or decrease in methylation

CG_hal <- ifelse(dmrseq_output_CG_hal_sig$stat>0, "decrease", "increase")
CHG_hal <- ifelse(dmrseq_output_CHG_hal_sig$stat>0, "decrease", "increase")
CHH_hal <- ifelse(dmrseq_output_CHH_hal_sig$stat>0, "decrease", "increase")
CG_lyr <- ifelse(dmrseq_output_CG_lyr_sig$stat>0, "decrease", "increase")
CHG_lyr <- ifelse(dmrseq_output_CHG_lyr_sig$stat>0, "decrease", "increase")
CHH_lyr <- ifelse(dmrseq_output_CHH_lyr_sig$stat>0, "decrease", "increase")

## To create the barplot, we first generate a dataframe with all
## the information we need about changes from both sides

context <- rep(c("CG", "CHG", "CHH"), 4)

species <- c(rep("halleri", 6), rep("lyrata", 6))

hypo_hyper <- c(rep("hyper", 3), rep("hypo", 3), rep("hyper", 3), rep("hypo", 3))

value <- c(table(CG_hal)[2], table(CHG_hal)[2], table(CHH_hal)[2],
           table(CG_hal)[1], table(CHG_hal)[1], table(CHH_hal)[1],
           table(CG_lyr)[2], table(CHG_lyr)[2], table(CHH_lyr)[2],
           table(CG_lyr)[1], table(CHG_lyr)[1], table(CHH_lyr)[1])

comparison <- c(rep("synthetic_1", 12))

dmrs_context_total_new <- data.frame(context, species, hypo_hyper, value, comparison)

dmrs_context_total <- rbind(dmrs_context_total, dmrs_context_total_new)

### We import the data one last time
### natural vs synthetic generation 4

setwd("ARPEGGIO_results_alk1Vsyn4_LL/dmrseq")
dmrseq_output_CG_AvB <- fread("CG_context/A_v_B_polyploid.txt")
dmrseq_output_CHG_AvB <- fread("CHG_context/A_v_B_polyploid.txt")
dmrseq_output_CHH_AvB <- fread("CHH_context/A_v_B_polyploid.txt")

## Split the files on two sides

scaffold_number_CG <- as.integer(str_split_fixed(dmrseq_output_CG_AvB$seqnames, "_", 2)[,2])
scaffold_number_CHG <- as.integer(str_split_fixed(dmrseq_output_CHG_AvB$seqnames, "_", 2)[,2])
scaffold_number_CHH <- as.integer(str_split_fixed(dmrseq_output_CHH_AvB$seqnames, "_", 2)[,2])

dmrseq_output_CG_hal <- dmrseq_output_CG_AvB[scaffold_number_CG < 2240]
dmrseq_output_CG_lyr <- dmrseq_output_CG_AvB[scaffold_number_CG > 2239]
dmrseq_output_CHG_hal <- dmrseq_output_CHG_AvB[scaffold_number_CHG < 2240]
dmrseq_output_CHG_lyr <- dmrseq_output_CHG_AvB[scaffold_number_CHG > 2240]
dmrseq_output_CHH_hal <- dmrseq_output_CHH_AvB[scaffold_number_CHH < 2240]
dmrseq_output_CHH_lyr <- dmrseq_output_CHH_AvB[scaffold_number_CHH > 2240]


## Filter only significant regions

dmrseq_output_CG_hal_sig <- dmrseq_output_CG_hal[dmrseq_output_CG_hal$qval<0.05]
dmrseq_output_CHG_hal_sig <- dmrseq_output_CHG_hal[dmrseq_output_CHG_hal$qval<0.05]
dmrseq_output_CHH_hal_sig <- dmrseq_output_CHH_hal[dmrseq_output_CHH_hal$qval<0.05]
dmrseq_output_CG_lyr_sig <- dmrseq_output_CG_lyr[dmrseq_output_CG_lyr$qval<0.05]
dmrseq_output_CHG_lyr_sig <- dmrseq_output_CHG_lyr[dmrseq_output_CHG_lyr$qval<0.05]
dmrseq_output_CHH_lyr_sig <- dmrseq_output_CHH_lyr[dmrseq_output_CHH_lyr$qval<0.05]

## Filter low coverage scaffolds

matchCG <- dmrseq_output_CG_hal_sig$seqnames %in% LL_hal_lowC_scaffolds$V1
dmrseq_output_CG_hal_sig <- dmrseq_output_CG_hal_sig[!matchCG,]
matchCHG <- dmrseq_output_CHG_hal_sig$seqnames %in% LL_hal_lowC_scaffolds$V1
dmrseq_output_CHG_hal_sig <- dmrseq_output_CHG_hal_sig[!matchCHG,]
matchCHH <- dmrseq_output_CHH_hal_sig$seqnames %in% LL_hal_lowC_scaffolds$V1
dmrseq_output_CHH_hal_sig <- dmrseq_output_CHH_hal_sig[!matchCHH,]

## Count regions showing increase or decrease in methylation

CG_hal <- ifelse(dmrseq_output_CG_hal_sig$stat>0, "decrease", "increase")
CHG_hal <- ifelse(dmrseq_output_CHG_hal_sig$stat>0, "decrease", "increase")
CHH_hal <- ifelse(dmrseq_output_CHH_hal_sig$stat>0, "decrease", "increase")
CG_lyr <- ifelse(dmrseq_output_CG_lyr_sig$stat>0, "decrease", "increase")
CHG_lyr <- ifelse(dmrseq_output_CHG_lyr_sig$stat>0, "decrease", "increase")
CHH_lyr <- ifelse(dmrseq_output_CHH_lyr_sig$stat>0, "decrease", "increase")

## To create the barplot, we first generate a dataframe with all
## the information we need about changes from both sides

context <- rep(c("CG", "CHG", "CHH"), 4)

species <- c(rep("halleri", 6), rep("lyrata", 6))

hypo_hyper <- c(rep("hyper", 3), rep("hypo", 3), rep("hyper", 3), rep("hypo", 3))

value <- c(table(CG_hal)[2], table(CHG_hal)[2], table(CHH_hal)[2],
           table(CG_hal)[1], table(CHG_hal)[1], table(CHH_hal)[1],
           table(CG_lyr)[2], table(CHG_lyr)[2], table(CHH_lyr)[2],
           table(CG_lyr)[1], table(CHG_lyr)[1], table(CHH_lyr)[1])

comparison <- c(rep("synthetic_4", 12))

dmrs_context_total_new <- data.frame(context, species, hypo_hyper, value, comparison)

dmrs_context_total <- rbind(dmrs_context_total, dmrs_context_total_new)

## To check the dataframe used with ggplot
## View(dmrs_context_total)
## Barplot code

## We create a vector to change the labels on top of the plot
## i.e. from halleri to halleri-side

species.labs <- c(paste0("halleri-side"),
                  paste0("lyrata-side"))

names(species.labs) <- c("halleri", "lyrata")

## Second plot showing the proportion of hyper and hypo methylated regions

cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

ggplot(dmrs_context_total, aes(x=factor(interaction(context, comparison), levels = c("CG.progenitors",
                                                                                     "CG.synthetic_1",
                                                                                     "CG.synthetic_4",
                                                                                     "CHG.progenitors",
                                                                                     "CHG.synthetic_1",
                                                                                     "CHG.synthetic_4",
                                                                                     "CHH.progenitors",
                                                                                     "CHH.synthetic_1",
                                                                                     "CHH.synthetic_4")), y=value, fill = comparison, col = factor(hypo_hyper, levels = c("hyper", "hypo")))) +
  geom_col(width = 0.5, show.legend = TRUE, position = "stack") +
  scale_fill_manual(values=cbPalette) +
  scale_color_manual(labels = c("hypo", "hyper"), values=c("white", "black")) +
  facet_grid(~species, labeller = labeller(species = species.labs)) +
  ylab("Number of DMRs") +
  theme_bw() +
  theme(legend.title = element_blank(), text = element_text(size=25),
        legend.text = element_text(size = 15),
        strip.background = element_rect(color="white", fill="white", size=1.5, linetype="solid"),
        axis.ticks = element_blank()) +
  scale_x_discrete(name = "DMR context", breaks = ,
                   labels = c("", "CG", "", "", "CHG", "", "", "CHH", "")) +
  ylim(0, 35000) +
  geom_vline(xintercept = c(3.5, 6.5), linetype = "longdash") +
  ggtitle("ALK (Hot Conditions)")

# Simple line plot for trend over comparisons

# We modify the original dataframe so that comparisons are integers instead of strings

dmrs_context_total_new <- data.frame(context = rep(c("CG", "CHG", "CHH"), 6),
                                     side = c(rep(c(rep("halleri", 3),
                                                    rep("lyrata", 3)), 3)),
                                     value = c(dmrs_context_total$value[1] + dmrs_context_total$value[4],
                                               dmrs_context_total$value[2] + dmrs_context_total$value[5],
                                               dmrs_context_total$value[3] + dmrs_context_total$value[6],
                                               dmrs_context_total$value[7] + dmrs_context_total$value[10],
                                               dmrs_context_total$value[8] + dmrs_context_total$value[11],
                                               dmrs_context_total$value[9] + dmrs_context_total$value[12],
                                               dmrs_context_total$value[13] + dmrs_context_total$value[16],
                                               dmrs_context_total$value[14] + dmrs_context_total$value[17],
                                               dmrs_context_total$value[15] + dmrs_context_total$value[18],
                                               dmrs_context_total$value[19] + dmrs_context_total$value[22],
                                               dmrs_context_total$value[20] + dmrs_context_total$value[23],
                                               dmrs_context_total$value[21] + dmrs_context_total$value[24],
                                               dmrs_context_total$value[25] + dmrs_context_total$value[28],
                                               dmrs_context_total$value[26] + dmrs_context_total$value[29],
                                               dmrs_context_total$value[27] + dmrs_context_total$value[30],
                                               dmrs_context_total$value[31] + dmrs_context_total$value[34],
                                               dmrs_context_total$value[32] + dmrs_context_total$value[35],
                                               dmrs_context_total$value[33] + dmrs_context_total$value[36]),
                                     comparison = c(rep("NvsP1", 6), rep("NvsS1", 6), rep("NvsS4",6)))


LL_ALK <- ggplot(dmrs_context_total_new, aes(x=comparison, y=value, group = context, colour = context)) +
  geom_line(show.legend = TRUE, size = 2) +
  #scale_x_continuous(labels = c("pro", "G1", "G2", "G3", "G4")) +
  #scale_fill_manual(values=cbPalette) +
  #scale_color_manual(labels = c("hypo", "hyper"), values=c("white", "black")) +
  facet_grid(~side) +
  ylab("") +
  xlab("Generation") +
  theme_bw() +
  theme(legend.title = element_blank(), text = element_text(size=15),
        legend.text = element_text(size = 10),
        #strip.background = element_rect(color="white", fill="white", size=1.5, linetype="solid"),
        axis.ticks = element_blank()) +
  ggtitle("ALK (Hot Conditions)") +
  ylim(0, 35000)

#####################################################################
# Repeat the whole procedure from above to compare TKS vs all (LL) ##
#####################################################################

### Import first set of files, natural TKS vs progenitors (LL)

setwd("ARPEGGIO_results_tks1Vpro1_LL/dmrseq")
dmrseq_output_CG_hal <- fread("CG_context/parent1_v_allo.txt")
dmrseq_output_CHG_hal <- fread("CHG_context/parent1_v_allo.txt")
dmrseq_output_CHH_hal <- fread("CHH_context/parent1_v_allo.txt")
dmrseq_output_CG_lyr <- fread("CG_context/parent2_v_allo.txt")
dmrseq_output_CHG_lyr <- fread("CHG_context/parent2_v_allo.txt")
dmrseq_output_CHH_lyr <- fread("CHH_context/parent2_v_allo.txt")


## Filter only significant regions

dmrseq_output_CG_hal_sig <- dmrseq_output_CG_hal[dmrseq_output_CG_hal$qval<0.05]
dmrseq_output_CHG_hal_sig <- dmrseq_output_CHG_hal[dmrseq_output_CHG_hal$qval<0.05]
dmrseq_output_CHH_hal_sig <- dmrseq_output_CHH_hal[dmrseq_output_CHH_hal$qval<0.05]
dmrseq_output_CG_lyr_sig <- dmrseq_output_CG_lyr[dmrseq_output_CG_lyr$qval<0.05]
dmrseq_output_CHG_lyr_sig <- dmrseq_output_CHG_lyr[dmrseq_output_CHG_lyr$qval<0.05]
dmrseq_output_CHH_lyr_sig <- dmrseq_output_CHH_lyr[dmrseq_output_CHH_lyr$qval<0.05]

## Filter low coverage scaffolds

LL_hal_lowC_scaffolds  <- read.table("LL_RS7_G4_hal_LowCovScaffolds.txt", quote="\"", comment.char="")
LL_lyr_lowC_scaffolds  <- read.table("LL_RS7_G4_lyr_LowCovScaffolds.txt", quote="\"", comment.char="")

matchCG <- dmrseq_output_CG_hal_sig$seqnames %in% LL_hal_lowC_scaffolds$V1
dmrseq_output_CG_hal_sig <- dmrseq_output_CG_hal_sig[!matchCG,]
matchCHG <- dmrseq_output_CHG_hal_sig$seqnames %in% LL_hal_lowC_scaffolds$V1
dmrseq_output_CHG_hal_sig <- dmrseq_output_CHG_hal_sig[!matchCHG,]
matchCHH <- dmrseq_output_CHH_hal_sig$seqnames %in% LL_hal_lowC_scaffolds$V1
dmrseq_output_CHH_hal_sig <- dmrseq_output_CHH_hal_sig[!matchCHH,]

## Count regions showing increase or decrease in methylation

CG_hal <- ifelse(dmrseq_output_CG_hal_sig$stat>0, "decrease", "increase")
CHG_hal <- ifelse(dmrseq_output_CHG_hal_sig$stat>0, "decrease", "increase")
CHH_hal <- ifelse(dmrseq_output_CHH_hal_sig$stat>0, "decrease", "increase")
CG_lyr <- ifelse(dmrseq_output_CG_lyr_sig$stat>0, "decrease", "increase")
CHG_lyr <- ifelse(dmrseq_output_CHG_lyr_sig$stat>0, "decrease", "increase")
CHH_lyr <- ifelse(dmrseq_output_CHH_lyr_sig$stat>0, "decrease", "increase")

## To create the barplot, we first generate a dataframe with all
## the information we need about changes from both sides

context <- rep(c("CG", "CHG", "CHH"), 4)

species <- c(rep("halleri", 6), rep("lyrata", 6))

hypo_hyper <- c(rep("hyper", 3), rep("hypo", 3), rep("hyper", 3), rep("hypo", 3))

value <- c(table(CG_hal)[2], table(CHG_hal)[2], table(CHH_hal)[2],
           table(CG_hal)[1], table(CHG_hal)[1], table(CHH_hal)[1],
           table(CG_lyr)[2], table(CHG_lyr)[2], table(CHH_lyr)[2],
           table(CG_lyr)[1], table(CHG_lyr)[1], table(CHH_lyr)[1])

comparison <- c(rep("progenitors", 12))

dmrs_context_total <- data.frame(context, species, hypo_hyper, value, comparison)

### Now we import the second set of files and repeat the same operations
### At the end we will merge everything with the previous data frame

## Import data

setwd("ARPEGGIO_results_tks1Vsyn1_LL/dmrseq")
dmrseq_output_CG_AvB <- fread("CG_context/A_v_B_polyploid.txt")
dmrseq_output_CHG_AvB <- fread("CHG_context/A_v_B_polyploid.txt")
dmrseq_output_CHH_AvB <- fread("CHH_context/A_v_B_polyploid.txt")

## Split the files on two sides

scaffold_number_CG <- as.integer(str_split_fixed(dmrseq_output_CG_AvB$seqnames, "_", 2)[,2])
scaffold_number_CHG <- as.integer(str_split_fixed(dmrseq_output_CHG_AvB$seqnames, "_", 2)[,2])
scaffold_number_CHH <- as.integer(str_split_fixed(dmrseq_output_CHH_AvB$seqnames, "_", 2)[,2])

dmrseq_output_CG_hal <- dmrseq_output_CG_AvB[scaffold_number_CG < 2240]
dmrseq_output_CG_lyr <- dmrseq_output_CG_AvB[scaffold_number_CG > 2239]
dmrseq_output_CHG_hal <- dmrseq_output_CHG_AvB[scaffold_number_CHG < 2240]
dmrseq_output_CHG_lyr <- dmrseq_output_CHG_AvB[scaffold_number_CHG > 2240]
dmrseq_output_CHH_hal <- dmrseq_output_CHH_AvB[scaffold_number_CHH < 2240]
dmrseq_output_CHH_lyr <- dmrseq_output_CHH_AvB[scaffold_number_CHH > 2240]


## Filter only significant regions

dmrseq_output_CG_hal_sig <- dmrseq_output_CG_hal[dmrseq_output_CG_hal$qval<0.05]
dmrseq_output_CHG_hal_sig <- dmrseq_output_CHG_hal[dmrseq_output_CHG_hal$qval<0.05]
dmrseq_output_CHH_hal_sig <- dmrseq_output_CHH_hal[dmrseq_output_CHH_hal$qval<0.05]
dmrseq_output_CG_lyr_sig <- dmrseq_output_CG_lyr[dmrseq_output_CG_lyr$qval<0.05]
dmrseq_output_CHG_lyr_sig <- dmrseq_output_CHG_lyr[dmrseq_output_CHG_lyr$qval<0.05]
dmrseq_output_CHH_lyr_sig <- dmrseq_output_CHH_lyr[dmrseq_output_CHH_lyr$qval<0.05]

## Filter low coverage scaffolds

matchCG <- dmrseq_output_CG_hal_sig$seqnames %in% LL_hal_lowC_scaffolds$V1
dmrseq_output_CG_hal_sig <- dmrseq_output_CG_hal_sig[!matchCG,]
matchCHG <- dmrseq_output_CHG_hal_sig$seqnames %in% LL_hal_lowC_scaffolds$V1
dmrseq_output_CHG_hal_sig <- dmrseq_output_CHG_hal_sig[!matchCHG,]
matchCHH <- dmrseq_output_CHH_hal_sig$seqnames %in% LL_hal_lowC_scaffolds$V1
dmrseq_output_CHH_hal_sig <- dmrseq_output_CHH_hal_sig[!matchCHH,]

## Count regions showing increase or decrease in methylation

CG_hal <- ifelse(dmrseq_output_CG_hal_sig$stat>0, "decrease", "increase")
CHG_hal <- ifelse(dmrseq_output_CHG_hal_sig$stat>0, "decrease", "increase")
CHH_hal <- ifelse(dmrseq_output_CHH_hal_sig$stat>0, "decrease", "increase")
CG_lyr <- ifelse(dmrseq_output_CG_lyr_sig$stat>0, "decrease", "increase")
CHG_lyr <- ifelse(dmrseq_output_CHG_lyr_sig$stat>0, "decrease", "increase")
CHH_lyr <- ifelse(dmrseq_output_CHH_lyr_sig$stat>0, "decrease", "increase")

## To create the barplot, we first generate a dataframe with all
## the information we need about changes from both sides

context <- rep(c("CG", "CHG", "CHH"), 4)

species <- c(rep("halleri", 6), rep("lyrata", 6))

hypo_hyper <- c(rep("hyper", 3), rep("hypo", 3), rep("hyper", 3), rep("hypo", 3))

value <- c(table(CG_hal)[2], table(CHG_hal)[2], table(CHH_hal)[2],
           table(CG_hal)[1], table(CHG_hal)[1], table(CHH_hal)[1],
           table(CG_lyr)[2], table(CHG_lyr)[2], table(CHH_lyr)[2],
           table(CG_lyr)[1], table(CHG_lyr)[1], table(CHH_lyr)[1])

comparison <- c(rep("synthetic_1", 12))

dmrs_context_total_new <- data.frame(context, species, hypo_hyper, value, comparison)

dmrs_context_total <- rbind(dmrs_context_total, dmrs_context_total_new)

### We import the data one last time
### natural vs synthetic generation 4

setwd("ARPEGGIO_results_tks1Vsyn4_LL/dmrseq")
dmrseq_output_CG_AvB <- fread("CG_context/A_v_B_polyploid.txt")
dmrseq_output_CHG_AvB <- fread("CHG_context/A_v_B_polyploid.txt")
dmrseq_output_CHH_AvB <- fread("CHH_context/A_v_B_polyploid.txt")

## Split the files on two sides

scaffold_number_CG <- as.integer(str_split_fixed(dmrseq_output_CG_AvB$seqnames, "_", 2)[,2])
scaffold_number_CHG <- as.integer(str_split_fixed(dmrseq_output_CHG_AvB$seqnames, "_", 2)[,2])
scaffold_number_CHH <- as.integer(str_split_fixed(dmrseq_output_CHH_AvB$seqnames, "_", 2)[,2])

dmrseq_output_CG_hal <- dmrseq_output_CG_AvB[scaffold_number_CG < 2240]
dmrseq_output_CG_lyr <- dmrseq_output_CG_AvB[scaffold_number_CG > 2239]
dmrseq_output_CHG_hal <- dmrseq_output_CHG_AvB[scaffold_number_CHG < 2240]
dmrseq_output_CHG_lyr <- dmrseq_output_CHG_AvB[scaffold_number_CHG > 2240]
dmrseq_output_CHH_hal <- dmrseq_output_CHH_AvB[scaffold_number_CHH < 2240]
dmrseq_output_CHH_lyr <- dmrseq_output_CHH_AvB[scaffold_number_CHH > 2240]


## Filter only significant regions

dmrseq_output_CG_hal_sig <- dmrseq_output_CG_hal[dmrseq_output_CG_hal$qval<0.05]
dmrseq_output_CHG_hal_sig <- dmrseq_output_CHG_hal[dmrseq_output_CHG_hal$qval<0.05]
dmrseq_output_CHH_hal_sig <- dmrseq_output_CHH_hal[dmrseq_output_CHH_hal$qval<0.05]
dmrseq_output_CG_lyr_sig <- dmrseq_output_CG_lyr[dmrseq_output_CG_lyr$qval<0.05]
dmrseq_output_CHG_lyr_sig <- dmrseq_output_CHG_lyr[dmrseq_output_CHG_lyr$qval<0.05]
dmrseq_output_CHH_lyr_sig <- dmrseq_output_CHH_lyr[dmrseq_output_CHH_lyr$qval<0.05]

## Filter low coverage scaffolds

matchCG <- dmrseq_output_CG_hal_sig$seqnames %in% LL_hal_lowC_scaffolds$V1
dmrseq_output_CG_hal_sig <- dmrseq_output_CG_hal_sig[!matchCG,]
matchCHG <- dmrseq_output_CHG_hal_sig$seqnames %in% LL_hal_lowC_scaffolds$V1
dmrseq_output_CHG_hal_sig <- dmrseq_output_CHG_hal_sig[!matchCHG,]
matchCHH <- dmrseq_output_CHH_hal_sig$seqnames %in% LL_hal_lowC_scaffolds$V1
dmrseq_output_CHH_hal_sig <- dmrseq_output_CHH_hal_sig[!matchCHH,]

## Count regions showing increase or decrease in methylation

CG_hal <- ifelse(dmrseq_output_CG_hal_sig$stat>0, "decrease", "increase")
CHG_hal <- ifelse(dmrseq_output_CHG_hal_sig$stat>0, "decrease", "increase")
CHH_hal <- ifelse(dmrseq_output_CHH_hal_sig$stat>0, "decrease", "increase")
CG_lyr <- ifelse(dmrseq_output_CG_lyr_sig$stat>0, "decrease", "increase")
CHG_lyr <- ifelse(dmrseq_output_CHG_lyr_sig$stat>0, "decrease", "increase")
CHH_lyr <- ifelse(dmrseq_output_CHH_lyr_sig$stat>0, "decrease", "increase")

## To create the barplot, we first generate a dataframe with all
## the information we need about changes from both sides

context <- rep(c("CG", "CHG", "CHH"), 4)

species <- c(rep("halleri", 6), rep("lyrata", 6))

hypo_hyper <- c(rep("hyper", 3), rep("hypo", 3), rep("hyper", 3), rep("hypo", 3))

value <- c(table(CG_hal)[2], table(CHG_hal)[2], table(CHH_hal)[2],
           table(CG_hal)[1], table(CHG_hal)[1], table(CHH_hal)[1],
           table(CG_lyr)[2], table(CHG_lyr)[2], table(CHH_lyr)[2],
           table(CG_lyr)[1], table(CHG_lyr)[1], table(CHH_lyr)[1])

comparison <- c(rep("synthetic_4", 12))

dmrs_context_total_new <- data.frame(context, species, hypo_hyper, value, comparison)

dmrs_context_total <- rbind(dmrs_context_total, dmrs_context_total_new)

## To check the dataframe used with ggplot
## View(dmrs_context_total)
## Barplot code

## We create a vector to change the labels on top of the plot
## i.e. from halleri to halleri-side

species.labs <- c(paste0("halleri-side"),
                  paste0("lyrata-side"))

names(species.labs) <- c("halleri", "lyrata")

## Second plot showing the proportion of hyper and hypo methylated regions

cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

ggplot(dmrs_context_total, aes(x=factor(interaction(context, comparison), levels = c("CG.progenitors",
                                                                                     "CG.synthetic_1",
                                                                                     "CG.synthetic_4",
                                                                                     "CHG.progenitors",
                                                                                     "CHG.synthetic_1",
                                                                                     "CHG.synthetic_4",
                                                                                     "CHH.progenitors",
                                                                                     "CHH.synthetic_1",
                                                                                     "CHH.synthetic_4")), y=value, fill = comparison, col = factor(hypo_hyper, levels = c("hyper", "hypo")))) +
  geom_col(width = 0.5, show.legend = TRUE, position = "stack") +
  scale_fill_manual(values=cbPalette) +
  scale_color_manual(labels = c("hypo", "hyper"), values=c("white", "black")) +
  facet_grid(~species, labeller = labeller(species = species.labs)) +
  ylab("Number of DMRs") +
  theme_bw() +
  theme(legend.title = element_blank(), text = element_text(size=25),
        legend.text = element_text(size = 15),
        strip.background = element_rect(color="white", fill="white", size=1.5, linetype="solid"),
        axis.ticks = element_blank()) +
  scale_x_discrete(name = "DMR context", breaks = ,
                   labels = c("", "CG", "", "", "CHG", "", "", "CHH", "")) +
  ylim(0, 35000) +
  geom_vline(xintercept = c(3.5, 6.5), linetype = "longdash") +
  ggtitle("TKS (Hot Conditions)")

# Simple line plot for trend over comparisons

# We modify the original dataframe so that comparisons are integers instead of strings

dmrs_context_total_new <- data.frame(context = rep(c("CG", "CHG", "CHH"), 6),
                                     side = c(rep(c(rep("halleri", 3),
                                                    rep("lyrata", 3)), 3)),
                                     value = c(dmrs_context_total$value[1] + dmrs_context_total$value[4],
                                               dmrs_context_total$value[2] + dmrs_context_total$value[5],
                                               dmrs_context_total$value[3] + dmrs_context_total$value[6],
                                               dmrs_context_total$value[7] + dmrs_context_total$value[10],
                                               dmrs_context_total$value[8] + dmrs_context_total$value[11],
                                               dmrs_context_total$value[9] + dmrs_context_total$value[12],
                                               dmrs_context_total$value[13] + dmrs_context_total$value[16],
                                               dmrs_context_total$value[14] + dmrs_context_total$value[17],
                                               dmrs_context_total$value[15] + dmrs_context_total$value[18],
                                               dmrs_context_total$value[19] + dmrs_context_total$value[22],
                                               dmrs_context_total$value[20] + dmrs_context_total$value[23],
                                               dmrs_context_total$value[21] + dmrs_context_total$value[24],
                                               dmrs_context_total$value[25] + dmrs_context_total$value[28],
                                               dmrs_context_total$value[26] + dmrs_context_total$value[29],
                                               dmrs_context_total$value[27] + dmrs_context_total$value[30],
                                               dmrs_context_total$value[31] + dmrs_context_total$value[34],
                                               dmrs_context_total$value[32] + dmrs_context_total$value[35],
                                               dmrs_context_total$value[33] + dmrs_context_total$value[36]),
                                     comparison = c(rep("NvsP1", 6), rep("NvsS1", 6), rep("NvsS4",6)))


LL_TKS <- ggplot(dmrs_context_total_new, aes(x=comparison, y=value, group = context, colour = context)) +
  geom_line(show.legend = TRUE, size = 2) +
  #scale_x_continuous(labels = c("pro", "G1", "G2", "G3", "G4")) +
  #scale_fill_manual(values=cbPalette) +
  #scale_color_manual(labels = c("hypo", "hyper"), values=c("white", "black")) +
  facet_grid(~side) +
  ylab("") +
  xlab("Generation") +
  theme_bw() +
  theme(legend.title = element_blank(), text = element_text(size=15),
        legend.text = element_text(size = 10),
        #strip.background = element_rect(color="white", fill="white", size=1.5, linetype="solid"),
        axis.ticks = element_blank()) +
  ggtitle("TKS (Hot Conditions)") +
  ylim(0, 35000)


(HM_ALK | LL_ALK) /
  (HM_TKS | LL_TKS)
