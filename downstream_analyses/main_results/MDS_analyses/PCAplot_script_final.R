### This script produces a PCA plot with
### genome-wide methylation information

### import libraries

library(data.table)
library(pcaMethods)
library(tidyverse)
library(SummarizedExperiment)
library(devtools)
library(patchwork)

### import data

HM_halC <- fread("HM_allhal_filtered.cov.gz")
LL_halC <- fread("LL_allhal_filtered.cov.gz")
HM_lyrC <- fread("HM_alllyr_filtered.cov.gz")
LL_lyrC <- fread("LL_alllyr_filtered.cov.gz")

names(HM_halC) <- c("coords",
                    "hal_G1_1", "counts_11",
                 "hal_G1_2", "counts_12",
                 "hal_G1_3", "counts_13",
                 "hal_G4_1", "counts_21",
                 "hal_G4_2", "counts_22",
                 "hal_G4_3", "counts_23",
                 "RS7_G1_1", "counts_31",
                 "RS7_G1_2", "counts_32",
                 "RS7_G1_3", "counts_33",
                 "RS7_G4_1", "counts_41",
                 "RS7_G4_2", "counts_42",
                 "RS7_G4_3", "counts_43",
                 "ALK_G1_1", "counts_51",
                 "ALK_G1_2", "counts_52",
                 "ALK_G1_3", "counts_53",
                 "ALK_G4_1", "counts_61",
                 "ALK_G4_2", "counts_62",
                 "ALK_G4_3", "counts_63",
                 "TKS_G1_1", "counts_71",
                 "TKS_G1_2", "counts_72",
                 "TKS_G1_3", "counts_73",
                 "TKS_G5_1", "counts_81",
                 "TKS_G5_2", "counts_82",
                 "TKS_G5_3", "counts_83")

names(LL_halC) <- c("coords",
                    "hal_G1_1", "counts_11",
                    "hal_G1_2", "counts_12",
                    "hal_G4_1", "counts_21",
                    "hal_G4_2", "counts_22",
                    "RS7_G1_1", "counts_31",
                    "RS7_G1_2", "counts_32",
                    "RS7_G1_3", "counts_33",
                    "RS7_G4_1", "counts_41",
                    "RS7_G4_2", "counts_42",
                    "RS7_G4_3", "counts_43",
                    "ALK_G1_1", "counts_51",
                    "ALK_G1_2", "counts_52",
                    "ALK_G1_3", "counts_53",
                    "ALK_G4_1", "counts_61",
                    "ALK_G4_2", "counts_62",
                    "ALK_G4_3", "counts_63",
                    "TKS_G1_1", "counts_71",
                    "TKS_G1_2", "counts_72",
                    "TKS_G1_3", "counts_73",
                    "TKS_G5_1", "counts_81",
                    "TKS_G5_2", "counts_82",
                    "TKS_G5_3", "counts_83")

names(HM_lyrC) <- c("coords",
                    "hal_G1_1", "counts_11",
                    "hal_G1_2", "counts_12",
                    "hal_G1_3", "counts_13",
                    "hal_G4_1", "counts_21",
                    "hal_G4_2", "counts_22",
                    "hal_G4_3", "counts_23",
                    "RS7_G1_1", "counts_31",
                    "RS7_G1_2", "counts_32",
                    "RS7_G1_3", "counts_33",
                    "RS7_G4_1", "counts_41",
                    "RS7_G4_2", "counts_42",
                    "RS7_G4_3", "counts_43",
                    "ALK_G1_1", "counts_51",
                    "ALK_G1_2", "counts_52",
                    "ALK_G1_3", "counts_53",
                    "ALK_G4_1", "counts_61",
                    "ALK_G4_2", "counts_62",
                    "ALK_G4_3", "counts_63",
                    "TKS_G1_1", "counts_71",
                    "TKS_G1_2", "counts_72",
                    "TKS_G1_3", "counts_73",
                    "TKS_G5_1", "counts_81",
                    "TKS_G5_2", "counts_82",
                    "TKS_G5_3", "counts_83")

names(LL_lyrC) <- c("coords",
                    "hal_G1_1", "counts_11",
                    "hal_G1_2", "counts_12",
                    "hal_G4_1", "counts_21",
                    "hal_G4_2", "counts_22",
                    "RS7_G1_1", "counts_31",
                    "RS7_G1_2", "counts_32",
                    "RS7_G1_3", "counts_33",
                    "RS7_G4_1", "counts_41",
                    "RS7_G4_2", "counts_42",
                    "RS7_G4_3", "counts_43",
                    "ALK_G1_1", "counts_51",
                    "ALK_G1_2", "counts_52",
                    "ALK_G1_3", "counts_53",
                    "ALK_G4_1", "counts_61",
                    "ALK_G4_2", "counts_62",
                    "ALK_G4_3", "counts_63",
                    "TKS_G1_1", "counts_71",
                    "TKS_G1_2", "counts_72",
                    "TKS_G1_3", "counts_73",
                    "TKS_G5_1", "counts_81",
                    "TKS_G5_2", "counts_82",
                    "TKS_G5_3", "counts_83")

LL_halC_only <- dplyr::select(LL_halC, 
                              hal_G1_1, 
                              hal_G1_2,
                              hal_G4_1, 
                              hal_G4_2,
                              RS7_G1_1, 
                              RS7_G1_2, 
                              RS7_G1_3,
                              RS7_G4_1, 
                              RS7_G4_2, 
                              RS7_G4_3,
                              ALK_G1_1, 
                              ALK_G1_2, 
                              ALK_G1_3,
                              ALK_G4_1, 
                              ALK_G4_2, 
                              ALK_G4_3,
                              TKS_G1_1, 
                              TKS_G1_2, 
                              TKS_G1_3,
                              TKS_G5_1, 
                              TKS_G5_2, 
                              TKS_G5_3)

HM_lyrC_only <- dplyr::select(HM_lyrC, 
                              hal_G1_1, 
                              hal_G1_2,
                              hal_G1_3,
                              hal_G4_1, 
                              hal_G4_2,
                              hal_G4_3,
                              RS7_G1_1, 
                              RS7_G1_2, 
                              RS7_G1_3,
                              RS7_G4_1, 
                              RS7_G4_2, 
                              RS7_G4_3,
                              ALK_G1_1, 
                              ALK_G1_2, 
                              ALK_G1_3,
                              ALK_G4_1, 
                              ALK_G4_2, 
                              ALK_G4_3,
                              TKS_G1_1, 
                              TKS_G1_2, 
                              TKS_G1_3,
                              TKS_G5_1, 
                              TKS_G5_2, 
                              TKS_G5_3)

LL_lyrC_only <- dplyr::select(LL_lyrC, 
                              hal_G1_1, 
                              hal_G1_2,
                              hal_G4_1, 
                              hal_G4_2,
                              RS7_G1_1, 
                              RS7_G1_2, 
                              RS7_G1_3,
                              RS7_G4_1, 
                              RS7_G4_2, 
                              RS7_G4_3,
                              ALK_G1_1, 
                              ALK_G1_2, 
                              ALK_G1_3,
                              ALK_G4_1, 
                              ALK_G4_2, 
                              ALK_G4_3,
                              TKS_G1_1, 
                              TKS_G1_2, 
                              TKS_G1_3,
                              TKS_G5_1, 
                              TKS_G5_2, 
                              TKS_G5_3)

HM_halC_only <- dplyr::select(HM_halC, 
                              hal_G1_1, 
                              hal_G1_2, 
                              hal_G1_3,
                              hal_G4_1, 
                              hal_G4_2, 
                              hal_G4_3,
                              RS7_G1_1, 
                              RS7_G1_2, 
                              RS7_G1_3,
                              RS7_G4_1, 
                              RS7_G4_2, 
                              RS7_G4_3,
                              ALK_G1_1, 
                              ALK_G1_2, 
                              ALK_G1_3,
                              ALK_G4_1, 
                              ALK_G4_2, 
                              ALK_G4_3,
                              TKS_G1_1, 
                              TKS_G1_2, 
                              TKS_G1_3,
                              TKS_G5_1, 
                              TKS_G5_2, 
                              TKS_G5_3)



### mds plot
## halleri

HM_hal_methsTR <- asin(2 * HM_halC_only/100 - 1)

# top value can be changes

HM_hal_mds_meth <- limma::plotMDS(HM_hal_methsTR, top = 100000, 
                           plot = FALSE)$cmdscale.out

HM_df_hal <- data.frame(dim1 = HM_hal_mds_meth[,1], 
                 dim2 = HM_hal_mds_meth[,2], 
                 names = colnames(HM_halC_only), 
                 treat = c(rep("hal_G1", 3),
                           rep("hal_G4", 3),
                           rep("RS7_G1", 3),
                           rep("RS7_G4", 3),
                           rep("ALK_G1", 3),
                           rep("ALK_G4", 3),
                           rep("TKS_G1", 3),
                           rep("TKS_G5", 3)))

HM_mds_hal <- ggplot() + geom_point(data = HM_df_hal, mapping = aes_(x = ~dim1, 
                                                y = ~dim2, 
                                                color = ~treat), 
                      size = 5) +
  theme_bw() + 
  ggtitle("Cold conditions - halleri side") +
  theme(legend.text = element_text(size = 15),
        legend.title = element_blank(),
        axis.text = element_text(size = 16),
        axis.title = element_text(size = 16))

LL_hal_methsTR <- asin(2 * LL_halC_only/100 - 1)

LL_hal_mds_meth <- limma::plotMDS(LL_hal_methsTR, top = 10000, 
                                  plot = FALSE)$cmdscale.out

LL_df_hal <- data.frame(dim1 = LL_hal_mds_meth[,1], 
                        dim2 = LL_hal_mds_meth[,2], 
                        names = colnames(LL_halC_only), 
                        treat = c(rep("hal_G1", 2),
                                  rep("hal_G4", 2),
                                  rep("RS7_G1", 3),
                                  rep("RS7_G4", 3),
                                  rep("ALK_G1", 3),
                                  rep("ALK_G4", 3),
                                  rep("TKS_G1", 3),
                                  rep("TKS_G5", 3)))

LL_mds_hal <- ggplot() + geom_point(data = LL_df_hal, mapping = aes_(x = ~dim1, 
                                                                     y = ~dim2, 
                                                                     color = ~treat), 
                                    size = 5) +
  theme_bw() + 
  ggtitle("Hot conditions - halleri side") +
  theme(legend.text = element_text(size = 15),
        legend.title = element_blank(),
        axis.text = element_text(size = 16),
        axis.title = element_text(size = 16))


all_hal_methsTR <- asin(2 * all_halC_only/100 - 1)

all_hal_mds_meth <- limma::plotMDS(all_hal_methsTR, top = 10000, 
                                  plot = FALSE)$cmdscale.out

all_df_hal <- data.frame(dim1 = all_hal_mds_meth[,1], 
                        dim2 = all_hal_mds_meth[,2], 
                        names = colnames(all_halC_only), 
                        treat = c(rep("HM_hal_G1", 3),
                                  rep("HM_hal_G4", 3),
                                  rep("HM_RS7_G1", 3),
                                  rep("HM_RS7_G4", 3),
                                  rep("HM_ALK_G1", 3),
                                  rep("HM_ALK_G4", 3),
                                  rep("HM_TKS_G1", 3),
                                  rep("HM_TKS_G5", 3),
                                  rep("LL_hal_G1", 2),
                                  rep("LL_hal_G4", 2),
                                  rep("LL_RS7_G1", 3),
                                  rep("LL_RS7_G4", 3),
                                  rep("LL_ALK_G1", 3),
                                  rep("LL_ALK_G4", 3),
                                  rep("LL_TKS_G1", 3),
                                  rep("LL_TKS_G5", 3)))

all_mds_hal <- ggplot() + geom_point(data = all_df_hal, mapping = aes_(x = ~dim1, 
                                                                     y = ~dim2, 
                                                                     color = ~treat), 
                                    size = 5) +
  theme_bw() + 
  ggtitle("All conditions - halleri side") +
  theme(legend.text = element_text(size = 15),
        legend.title = element_blank(),
        axis.text = element_text(size = 16),
        axis.title = element_text(size = 16))

## lyrata

HM_lyr_methsTR <- asin(2 * HM_lyrC_only/100 - 1)

HM_lyr_mds_meth <- limma::plotMDS(HM_lyr_methsTR, top = 10000, 
                                  plot = FALSE)$cmdscale.out

HM_df_lyr <- data.frame(dim1 = HM_lyr_mds_meth[,1], 
                        dim2 = HM_lyr_mds_meth[,2], 
                        names = colnames(HM_lyrC_only), 
                        treat = c(rep("lyr_G1", 3),
                                  rep("lyr_G4", 3),
                                  rep("RS7_G1", 3),
                                  rep("RS7_G4", 3),
                                  rep("ALK_G1", 3),
                                  rep("ALK_G4", 3),
                                  rep("TKS_G1", 3),
                                  rep("TKS_G5", 3)))

HM_mds_lyr <- ggplot() + geom_point(data = HM_df_lyr, mapping = aes_(x = ~dim1, 
                                                                     y = ~dim2, 
                                                                     color = ~treat), 
                                    size = 5) +
  theme_bw() + 
  ggtitle("Cold conditions - lyrata side") +
  theme(legend.text = element_text(size = 15),
        legend.title = element_blank(),
        axis.text = element_text(size = 16),
        axis.title = element_text(size = 16))

LL_lyr_methsTR <- asin(2 * LL_lyrC_only/100 - 1)

LL_lyr_mds_meth <- limma::plotMDS(LL_lyr_methsTR, top = 10000, 
                                  plot = FALSE)$cmdscale.out

LL_df_lyr <- data.frame(dim1 = LL_lyr_mds_meth[,1], 
                        dim2 = LL_lyr_mds_meth[,2], 
                        names = colnames(LL_lyrC_only), 
                        treat = c(rep("lyr_G1", 2),
                                  rep("lyr_G4", 2),
                                  rep("RS7_G1", 3),
                                  rep("RS7_G4", 3),
                                  rep("ALK_G1", 3),
                                  rep("ALK_G4", 3),
                                  rep("TKS_G1", 3),
                                  rep("TKS_G5", 3)))

LL_mds_lyr <- ggplot() + geom_point(data = LL_df_lyr, mapping = aes_(x = ~dim1, 
                                                                     y = ~dim2, 
                                                                     color = ~treat), 
                                    size = 5) +
  theme_bw() + 
  ggtitle("Hot conditions - lyrata side") +
  theme(legend.text = element_text(size = 15),
        legend.title = element_blank(),
        axis.text = element_text(size = 16),
        axis.title = element_text(size = 16))


### MDS plot main Figure

#halleri cold conditions

HM_hal_methsTR <- asin(2 * HM_halC_only/100 - 1)

HM_hal_mds_meth <- limma::plotMDS(HM_hal_methsTR, top = 100000, 
                                  plot = FALSE)$cmdscale.out

HM_df_hal <- data.frame(dim1 = HM_hal_mds_meth[,1], 
                        dim2 = HM_hal_mds_meth[,2], 
                        names = colnames(HM_halC_only), 
                        treat = c(rep("halleri", 6),
                                  rep("synthetic", 6),
                                  rep("natural-ALK", 6),
                                  rep("natural-TKS_G1", 6)))

HM_mds_hal <- ggplot() + geom_point(data = HM_df_hal, mapping = aes_(x = ~dim1, 
                                                                     y = ~dim2, 
                                                                     color = ~treat), 
                                    size = 7) +
  theme_bw() + 
  ggtitle("Cold conditions - halleri side") +
  theme(legend.text = element_text(size = 15),
        legend.title = element_blank(),
        axis.text = element_text(size = 16),
        axis.title = element_text(size = 16))

HM_mds_hal

# halleri hot conditions

LL_hal_methsTR <- asin(2 * LL_halC_only/100 - 1)

LL_hal_mds_meth <- limma::plotMDS(LL_hal_methsTR, top = 100000, 
                                  plot = FALSE)$cmdscale.out

LL_df_hal <- data.frame(dim1 = LL_hal_mds_meth[,1], 
                        dim2 = LL_hal_mds_meth[,2], 
                        names = colnames(LL_halC_only), 
                        treat = c(rep("halleri", 4),
                                  rep("synthetic", 6),
                                  rep("natural-ALK", 6),
                                  rep("natural-TKS_G1", 6)))

LL_mds_hal <- ggplot() + geom_point(data = LL_df_hal, mapping = aes_(x = ~dim1, 
                                                                     y = ~dim2, 
                                                                     color = ~treat), 
                                    size = 7) +
  theme_bw() + 
  ggtitle("Hot conditions - halleri side") +
  theme(legend.text = element_text(size = 15),
        legend.title = element_blank(),
        axis.text = element_text(size = 16),
        axis.title = element_text(size = 16))

LL_mds_hal


#lyrata cold conditions

HM_lyr_methsTR <- asin(2 * HM_lyrC_only/100 - 1)

HM_lyr_mds_meth <- limma::plotMDS(HM_lyr_methsTR, top = 100000, 
                                  plot = FALSE)$cmdscale.out

HM_df_lyr <- data.frame(dim1 = HM_lyr_mds_meth[,1], 
                        dim2 = HM_lyr_mds_meth[,2], 
                        names = colnames(HM_lyrC_only), 
                        treat = c(rep("lyrata", 6),
                                  rep("synthetic", 6),
                                  rep("natural-ALK", 6),
                                  rep("natural-TKS_G1", 6)))

HM_mds_lyr <- ggplot() + geom_point(data = HM_df_lyr, mapping = aes_(x = ~dim1, 
                                                                     y = ~dim2, 
                                                                     color = ~treat), 
                                    size = 7) +
  theme_bw() + 
  ggtitle("Cold conditions - lyrata side") +
  theme(legend.text = element_text(size = 15),
        legend.title = element_blank(),
        axis.text = element_text(size = 16),
        axis.title = element_text(size = 16))

HM_mds_lyr

# lyrata hot conditions

LL_lyr_methsTR <- asin(2 * LL_lyrC_only/100 - 1)

LL_lyr_mds_meth <- limma::plotMDS(LL_lyr_methsTR, top = 100000, 
                                  plot = FALSE)$cmdscale.out

LL_df_lyr <- data.frame(dim1 = LL_lyr_mds_meth[,1], 
                        dim2 = LL_lyr_mds_meth[,2], 
                        names = colnames(LL_lyrC_only), 
                        treat = c(rep("lyrata", 4),
                                  rep("synthetic", 6),
                                  rep("natural-ALK", 6),
                                  rep("natural-TKS_G1", 6)))

LL_mds_lyr <- ggplot() + geom_point(data = LL_df_lyr, mapping = aes_(x = ~dim1, 
                                                                     y = ~dim2, 
                                                                     color = ~treat), 
                                    size = 7) +
  theme_bw() + 
  ggtitle("Hot conditions - lyrata side") +
  theme(legend.text = element_text(size = 15),
        legend.title = element_blank(),
        axis.text = element_text(size = 16),
        axis.title = element_text(size = 16))

LL_mds_lyr
