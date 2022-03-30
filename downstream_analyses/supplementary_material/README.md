# Supplementary results

To find which script is related to which result:

## Chapter 2

 - `Supplementary Figure 1` was obtained with `HE_CoverageC_chromosomes.R` and `coverage_data`
 - `Supplementary Figure 2` was obtained with `allDMRsOverlapsFunctionalGeneRegions.R` and annotation files for _A. halleri_ ([here](https://datadryad.org/stash/dataset/doi:10.5061/dryad.gn4hh)) and _A. lyrata_ (here)
 - `Supplementary Figures 3 & 4` were obtained with `DMRsSpatialDistributionKaryploteR.R`, ARPEGGIO output data and `agp` files from `files`
 - `Supplementary Figure 5` was obtained with `VisualizingRemapping.R` and `agp` files from `files`
 - `Supplementary Figure 6` was obtained with [D-GENIES](http://dgenies.toulouse.inra.fr/)
 - `Supplementary Figures 7, 8, 9 & 10` were obtained with `EnrichmentPlot.R` and recalibrated methylation values from `rc_methylation` in `files`
 - `Supplementary Figure 11 & 12` were obtained with `PCAplot_script_final.R` (first part, see `main_results` folder)
 - `Supplementary Figure 13` was obtained with `DMRs_Barplot_synthetics.R`, ARPEGGIO output data from `files` and low coverage regions from `coverage_data` in `files`
 - `Supplementary Figure 14` was obtained with `DMRs_barPlot_overGenerations.R`, ARPEGGIO outputs to download from a link in `files` and low coverage regions from `fles/coverage_data`
 - `Supplementary Figure 15` was obtained with `allDMRsOverlapsFunctionalRegions.R`, annotation files for _A. halleri_ ([here](https://datadryad.org/stash/dataset/doi:10.5061/dryad.gn4hh)) and _A. lyrata_ (here) and ARPEGGIO outputs both available and to download from a link in `files`

 ## Chapter 3

 - `Supplementary Figure 1` was obtained with `RNAseq_edgeR_fullReport.Rmd`, count files (download here link TBA) and files with genes falling in low coverage regions in `main_results/low_coverageFiles`
 - `Supplementary Figure 2 & 3` were obtained with `DEG_overTime.R`, lists of DEGs in `main_results/DEGs` and file with genes falling in low coverage regions in `main_results/low_coverageFiles`
 - `Supplementary Figures 4 & 5` were obtained with `SpatialExpression_v2.R`, low coverage regions from `fles/coverage_data`, annotation files for _A. halleri_ ([here](https://datadryad.org/stash/dataset/doi:10.5061/dryad.gn4hh)) and _A. lyrata_ (here), lists of DEGs in `main_results/DEGs` and `agp` files from `files/agp_files`
 - `Supplementary Figure 6` was obtained with `HeatmapSpatialDistributionGenes_v2.R`, count files (download here link TBA) normalized to TPM counts with `TPM_RPKM_converter.R`, annotation files for _A. halleri_ ([here](https://datadryad.org/stash/dataset/doi:10.5061/dryad.gn4hh)) and _A. lyrata_ (here) and `agp` files from `files/agp_files`
 - `Supplementary Figures 7 & 8` were obtained with `RawDE_vs_RawDM_correlation.R`, ARPEGGIO output data from `files` and lists of DEGs in `main_results/DEGs`
 - `Supplementary Figure 9` was obtained with `chisquare.R`
 - `Supplementary Figures 10-16` were obtained with `fromATIDtoMethPlot.R`, a command line script using ARPEGGIO outputs and Bismark cov files to be generated with ARPEGGIO
 - `Supplementary Table 2` was obtained with `main_results/GenesOfInterest.R` with ARPEGGIO outputs available in `supplementary_material/files` and lists of DEGs in `DEGs`
