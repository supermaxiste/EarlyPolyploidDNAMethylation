## MDS analyses

The starting point of MDS analyses are `cov` files from Bismark which can be obtained from ARPEGGIO runs. If these files are not available, we provide intermediate files (from step 3 onward, see below) together with the final files used for MDS plotting (after performing all of the steps outlined below). To perform the analyses for all progenitors' sides and conditions, there are a total of 96 files. Each condition has 48 files and includes the following samples:

 - _A. halleri_ generation 1 (3 replicates)
 - _A. halleri_ generation 4 (3 replicates)
 - _A. lyrata_ generation 1 (3 replicates)
 - _A. lyrata_ generation 4 (3 replicates)
 - _A. kamchatica_ synthetic generation 1 (3 replicates for each side: 6 replicates total)
 - _A. kamchatica_ synthetic generation 4 (3 replicates for each side: 6 replicates total)
 - _A. kamchatica_ natural generation 1, Alaska line (3 replicates for each side: 6 replicates total)
 - _A. kamchatica_ natural generation 4, Alaska line (3 replicates for each side: 6 replicates total)
 - _A. kamchatica_ natural generation 1, Takashima line (3 replicates for each side: 6 replicates total)
 - _A. kamchatica_ natural generation 4, Takashima line (3 replicates for each side: 6 replicates total)

 The following steps provide a guide to analyse the data:

  1) After downloading all `cov` files, in order to merge them, all need to be sorted. Assuming that we're working with a file called `filename_bismark.cov.gz`, the command to sort it is:

  ```
  zcat filename_bismark.cov.gz | sort -k1,1 -k2,2n | pigz > filename_bismark.cov.gz
  ```

  This command should be run for all samples and for _A. kamchatica_ samples make sure to name the files differently for each progenitors' side, for example `filename_hal_bismark.cov.gz` and `filename_lyr_bismark.cov.gz`

  2) Next we reformat the columns from the `cov` files by merging the first there columns with the spatial coordinates (providing unique names to each position) and summing the last two columns to obtain the total count for both methylated and unmethylated cytosines (coverage). We do this as follows:

  ```
  data=path/to/*cov.gz
  output=output/path/

  for filename in $data; do
          base=$(basename $filename _bismark.cov.gz)
          zcat $filename | awk '{print $1"_"$2"_"$3"\t"$4"\t"$5+$6}' | sort > output/path/${base}_minimal.cov
          pigz output/path/${base}_minimal.cov
  done
  ```
  Make sure to modify `data` and `output` paths.

  3) Files can now be merged based on their overlapping cytosines. For each sample, two values will be kept: the methylation proportion and the total number of cytosines. Assuming we have files from a samples named `filename_1_lyr_minimal.cov`, `filename_2_lyr_minimal.cov` and `filename_3_lyr_minimal.cov`, we merge them as follows:

  ```
  join -j 1 -o 1.1,1.2,1.3,2.2,2.3 filename_1_lyr_minimal.cov filename_2_lyr_minimal.cov > filename_12_lyr_minimal.cov
  join -j 1 -o 1.1,1.2,1.3,1.4,1.5,2.2,2.3 filename_12_lyr_minimal.cov filename_3_lyr_minimal.cov > filename_123_lyr_minimal.cov
  ```
  The structure of the output file `filename_123_lyr_minimal.cov` will be as follows:

   - Column 1: cytosine coordinates
   - Column 2-3: methylation proportion and total cytosines for replicate 1
   - Column 4-5: methylation proportion and total cytosines for replicate 2
   - Column 6-7: methylation proportion and total cytosines for replicate 3

  Output files from this step for each sample in each condition can be obtained at the following link: TBA

   4) Merged files can now be merged further in four final files (2 conditions x 2 progenitors' sides). We provide an extract from our own scripts for one progenitor side and both conditions (HM = cold, LL = hot) since the number of replicates is different:

   ```
   # halleri side, HM conditions

   join -j 1 -o 1.1,1.2,1.3,1.4,1.5,1.6,1.7,2.2,2.3,2.4,2.5,2.6,2.7 HM_RS7_G1_123hal_minimal.cov HM_RS7_G4_123hal_minimal.cov > HM_synthe_allgen_hal_minimal.cov

   join -j 1 -o 1.1,1.2,1.3,1.4,1.5,1.6,1.7,2.2,2.3,2.4,2.5,2.6,2.7 HM_ALK_G1_123hal_minimal.cov HM_ALK_G4_123hal_minimal.cov > HM_ALK_allgen_hal_minimal.cov

   join -j 1 -o 1.1,1.2,1.3,1.4,1.5,1.6,1.7,2.2,2.3,2.4,2.5,2.6,2.7 HM_TKS_G1_123hal_minimal.cov HM_TKS_G5_123hal_minimal.cov > HM_TKS_allgen_hal_minimal.cov

   join -j 1 -o 1.1,1.2,1.3,1.4,1.5,1.6,1.7,2.2,2.3,2.4,2.5,2.6,2.7 HM_hal_G1_123_minimal.cov HM_hal_G4_123_minimal.cov > HM_pro_allgen_hal_minimal.cov

   join -j 1 -o 1.1,1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9,1.10,1.11,1.12,1.13,2.2,2.3,2.4,2.5,2.6,2.7,2.8,2.9,2.10,2.11,2.12,2.13 HM_pro_allgen_hal_minimal.cov HM_synthe_allgen_hal_minimal.cov > HM_prosyn_hal.cov

   join -j 1 -o 1.1,1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9,1.10,1.11,1.12,1.13,2.2,2.3,2.4,2.5,2.6,2.7,2.8,2.9,2.10,2.11,2.12,2.13 HM_ALK_allgen_hal_minimal.cov HM_TKS_allgen_hal_minimal.cov > HM_nat_hal.cov

   join -j 1 -o 1.1,1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9,1.10,1.11,1.12,1.13,1.14,1.15,1.16,1.17,1.18,1.19,1.20,1.21,1.22,1.23,1.24,1.25,2.2,2.3,2.4,2.5,2.6,2.7,2.8,2.9,2.10,2.11,2.12,2.13,2.14,2.15,2.16,2.17,2.18,2.19,2.20,2.21,2.22,2.23,2.24,2.25 HM_prosyn_hal.cov HM_nat_hal.cov > HM_allhal.cov

   #Final ordering in the file:
   #halleri_G1+G4  HM_RS7K_G1+G4   HM_ALK_G1+G4    HM_TKS_G1+5

   # halleri side, LL conditions

   join -j 1 -o 1.1,1.2,1.3,1.4,1.5,1.6,1.7,2.2,2.3,2.4,2.5,2.6,2.7 LL_RS7_G1_123hal_minimal.cov LL_RS7K_G4_123hal_minimal.cov > LL_synthe_allgen_hal_minimal.cov

   join -j 1 -o 1.1,1.2,1.3,1.4,1.5,1.6,1.7,2.2,2.3,2.4,2.5,2.6,2.7 LL_ALK_G1_123hal_minimal.cov LL_ALK_G4_123hal_minimal.cov > LL_ALK_allgen_hal_minimal.cov

   join -j 1 -o 1.1,1.2,1.3,1.4,1.5,1.6,1.7,2.2,2.3,2.4,2.5,2.6,2.7 LL_TKS_G1_123hal_minimal.cov LL_TKS_G5_123hal_minimal.cov > LL_TKS_allgen_hal_minimal.cov

   join -j 1 -o 1.1,1.2,1.3,1.4,1.5,2.2,2.3,2.4,2.5 LL_hal_G1_12_minimal.cov LL_hal_G4_12_minimal.cov > LL_pro_allgen_hal_minimal.cov

   join -j 1 -o 1.1,1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9,2.2,2.3,2.4,2.5,2.6,2.7,2.8,2.9,2.10,2.11,2.12,2.13 LL_pro_allgen_hal_minimal.cov LL_synthe_allgen_hal_minimal.cov > LL_prosyn_hal.cov

   join -j 1 -o 1.1,1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9,1.10,1.11,1.12,1.13,2.2,2.3,2.4,2.5,2.6,2.7,2.8,2.9,2.10,2.11,2.12,2.13 LL_ALK_allgen_hal_minimal.cov LL_TKS_allgen_hal_minimal.cov > LL_nat_hal.cov

   join -j 1 -o 1.1,1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9,1.10,1.11,1.12,1.13,1.14,1.15,1.16,1.17,1.18,1.19,1.20,1.21,1.22,1.23,1.24,1.25,2.2,2.3,2.4,2.5,2.6,2.7,2.8,2.9,2.10,2.11,2.12,2.13,2.14,2.15,2.16,2.17,2.18,2.19,2.20,2.21,2.22,2.23,2.24,2.25 LL_prosyn_hal.cov LL_nat_hal.cov > LL_allhal.cov

   #Final ordering in the file:
   #halleri_G1+G4  HM_RS7K_G1+G4   HM_ALK_G1+G4    HM_TKS_G1+5
   ```  

   It's important to keep in mind the way files are merged to know which column represents which sample. In all of our files the order is the following. For cold conditions:

    - Column 1: cytosine coordinates
    - Columns 2-7: Progenitor generation 1 (3 replicates)
    - Columns 8-13: Progenitor generation 4 (3 replicates)
    - Columns 14-19: _A. kamchatica_ synthetic generation 1 (3 replicates)
    - Columns 20-25: _A. kamchatica_ synthetic generation 4 (3 replicates)
    - Columns 26-31: _A. kamchatica_ natural generation 1, Alaska line (3 replicates)
    - Columns 32-37: _A. kamchatica_ natural generation 4, Alaska line (3 replicates)
    - Columns 38-43: _A. kamchatica_ natural generation 1, Takashima line (3 replicates)
    - Columns 44-49: _A. kamchatica_ natural generation 4, Takashima line (3 replicates)

   For hot conditions:

   - Column 1: cytosine coordinates
   - Columns 2-5: Progenitor generation 1 (2 replicates)
   - Columns 6-9: Progenitor generation 4 (2 replicates)
   - Columns 10-15: _A. kamchatica_ synthetic generation 1 (3 replicates)
   - Columns 16-21: _A. kamchatica_ synthetic generation 4 (3 replicates)
   - Columns 22-27: _A. kamchatica_ natural generation 1, Alaska line (3 replicates)
   - Columns 28-33: _A. kamchatica_ natural generation 4, Alaska line (3 replicates)
   - Columns 34-39: _A. kamchatica_ natural generation 1, Takashima line (3 replicates)
   - Columns 40-45: _A. kamchatica_ natural generation 4, Takashima line (3 replicates)


  5) The last step is to filter all cytosines that have >=3 coverage. We use the same output files from the previous step:

   ```
   #HM

   awk -F " " '{ if(($3 >= 3) && ($5 >= 3) && ($7 >= 3) && ($9 >= 3) && ($11 >= 3) && ($13 >= 3) && ($15 >= 3) && ($17 >= 3) && ($19 >= 3) && ($21 >= 3) && ($23 >= 3) && ($25 >= 3) && ($27 >= 3) && ($29 >= 3) && ($31 >= 3) && ($33 >= 3) && ($35 >= 3) && ($37 >= 3) && ($39 >= 3) && ($41 >= 3) && ($43 >= 3) && ($45 >= 3) && ($47 >= 3) && ($49 >= 3)) { print } }' ${data}HM_allhal.cov > HM_allhal_filtered.cov

   #LL

   awk -F " " '{ if(($3 >= 3) && ($5 >= 3) && ($7 >= 3) && ($9 >= 3) && ($11 >= 3) && ($13 >= 3) && ($15 >= 3) && ($17 >= 3) && ($19 >= 3) && ($21 >= 3) && ($23 >= 3) && ($25 >= 3) && ($27 >= 3) && ($29 >= 3) && ($31 >= 3) && ($33 >= 3) && ($35 >= 3) && ($37 >= 3) && ($39 >= 3) && ($41 >= 3) && ($43 >= 3) && ($45 >= 3))  { print } }' ${data}LL_allhal.cov > LL_allhal_filtered.cov
   ```

  The outputs from this step can be found at this link: TBA
