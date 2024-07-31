######################################################################################

This is the electronic appendix (R code) to the article "Does combining numerous data types in multi-omics data improve or hinder performance in survival prediction? Insights from a large-scale benchmark study" (2024) by Yingxia Li1, Tobias Herold2, Ulrich Mansmann1, Roman Hornung1,3

1 Institute for Medical Information Processing, Biometry and Epidemiology, University of Munich, Marchioninistr. 15, 81377 Munich, Germany;

2 Laboratory for Leukemia Diagnostics, Department of Medicine III, LMU University Hospital, LMU Munich, Munich, Germany;

3 Munich Center for Machine Learning (MCML), Munich, Germany;

######################################################################################

Rrogram and Platfprm
- Program: R, version 4.1.2 (2021-11-01)
- Used Platform: Linux (x86-64)  (for the conduction of the analyses)
                 Windows10 64-bit (for the evaluation of the results)
-  The following output from sessionInfo() describes which R packages and versions were used:
  
  sessionInfo()
R version 4.1.2 (2021-11-01)
Platform: x86_64-pc-linux-gnu (64-bit)
Running under: SUSE Linux Enterprise Server 15 SP1

Matrix products: default

locale:
 [1] LC_CTYPE=en_US.UTF-8       LC_NUMERIC=C
 [3] LC_TIME=en_US.UTF-8        LC_COLLATE=en_US.UTF-8
 [5] LC_MONETARY=en_US.UTF-8    LC_MESSAGES=en_US.UTF-8
 [7] LC_PAPER=en_US.UTF-8       LC_NAME=C
 [9] LC_ADDRESS=C               LC_TELEPHONE=C
[11] LC_MEASUREMENT=en_US.UTF-8 LC_IDENTIFICATION=C

attached base packages:                                                                                                                                
[1] grid      stats     graphics  grDevices utils     datasets  methods
[8] base

other attached packages:
 [1] mlr_2.19.2          tidyr_1.3.1         survcomp_1.44.1
 [4] stringr_1.5.1       snow_0.4-4          Rmisc_1.5.1
 [7] lattice_0.20-45     RColorBrewer_1.1-3  ranger_0.13.1
[10] prioritylasso_0.2.5 plyr_1.8.7          pec_2022.03.06
[13] prodlim_2019.11.13  patchwork_1.2.0     ParamHelpers_1.14.1
[16] OpenML_1.12         ipflasso_1.1        survival_3.3-1
[19] gridExtra_2.3       glmnet_4.1-3        Matrix_1.4-1
[22] ggplot2_3.5.1       farff_1.1.1         dplyr_1.1.4
[25] bootstrap_2019.6    blockForest_0.2.4

loaded via a namespace (and not attached):                                                                 
 [1] httr_1.4.7          jsonlite_1.8.8      splines_4.1.2
 [4] foreach_1.5.2       rmeta_3.0           globals_0.14.0
 [7] survivalROC_1.0.3   timereg_2.0.2       numDeriv_2016.8-1.1
[10] pillar_1.9.0        backports_1.4.1     glue_1.6.2
[13] digest_0.6.36       checkmate_2.0.0     colorspace_2.0-3
[16] XML_3.99-0.17       pkgconfig_2.0.3     listenv_0.8.0
[19] purrr_1.0.2         scales_1.3.0        parallelMap_1.5.1
[22] lava_1.6.10         tzdb_0.4.0          tibble_3.2.1
[25] generics_0.1.2      cachem_1.1.0        withr_2.5.0
[28] cli_3.6.3           magrittr_2.0.3      memoise_2.0.1
[31] future_1.24.0       fansi_1.0.3         parallelly_1.31.0
[34] SuppDists_1.1-9.7   tools_4.1.2         data.table_1.14.2
[37] hms_1.1.3           lifecycle_1.0.4     BBmisc_1.13
[40] munsell_0.5.0       compiler_4.1.2      rlang_1.1.4
[43] iterators_1.0.14    gtable_0.3.0        codetools_0.2-18
[46] curl_5.2.1          R6_2.5.1            fastmap_1.2.0
[49] future.apply_1.8.1  utf8_1.2.2          fastmatch_1.1-4
[52] KernSmooth_2.23-20  shape_1.4.6         readr_2.1.5
[55] stringi_1.7.6       parallel_4.1.2      Rcpp_1.0.8.3
[58] vctrs_0.6.5         tidyselect_1.2.1

######################################################################################

Repository Structure

1. Data Subfolder

##########################
- **Purpose:** Contains scripts for downloading the OpenML data needed for reproducing the analyses.
- **Contents:**
  - `down_data.R`: Downloads data from OpenML.
  - `dataset_ids.RData`: Includes OpenML IDs for the datasets to be downloaded.

2. JobScripts Subfolder

##########################
- **Purpose:** Contains scripts for reproducing the benchmark study.
- **Contents:**
  - `AnalysisCluster_1_4_5.R`: For single blocks, combinations of 4 blocks, and combinations of 5 blocks.
  - `AnalysisCluster_2.R`: For combinations of 2 blocks.
  - `AnalysisCluster_3.R`: For combinations of 3 blocks.

3. Functions Subfolder
   
##########################
- **Purpose:** Contains scripts with functions used in the benchmark study.
- **Contents:**
  - Scripts whose labels contain `Functions_AnalysisCluster`: Functions for applying the different prediction methods to the different combinations.

4. Evaluations Subfolder
   
##########################
- **Purpose:** Contains scripts for evaluating results and reproducing figures and tables.
- **Contents:**
  - `Evaluation_AnalysisCluster_fivemethods.R`: Evaluates raw results.
  - `bootstrap analysis_ibrier.R` and `bootstrap analysis_cindex.R`: Performs the bootstrap analysis.
  - Scripts labeled `figures`: Reproduces figures shown in the paper and supplement, where `figures_2_S6.R` in addition produces the results of the analysis presented at the end of Section "Best-performing combinations of prediction methods and blocks per dataset".
  - `test_for_figure_2.R` and `tests_for_table_3.R`: Performs the statistical tests for Figure 2 and Table 3.

5. Results Subfolder
   
##########################
- **Purpose:** Contains results and figures from the benchmark study.
- **Contents:**
  - **rda_files Subfolder:**
    - `scenariogrid1.Rda`, `scenariogrid2.Rda`, `scenariogrid3.Rda`: Generated by the R scripts in `JobScripts`.
    - `resultsumsum.RData`, `resultsum.RData`, `CI_cindex.xlsx`, `CI_ibrier.xlsx`: Generated by R scripts in `Evaluations`.
  - **Figures Subfolder:**
    - Contains all figures from the paper and supplement. Most figures have two versions, one with the suffix "_raw" and one without. The versions with the suffix "_raw" were generated by the R code, and the versions without the suffix were subsequently edited for visual reasons.

 6. Features_1000
    
##########################
-   Purpose:  as a sensitivity analysis, just change the selected featrues from 2500 into 1000.

Further details on the benchmark study are available in the paper.

7. Full reproduction of the results:
   
#################################

- All R code needed to fully reproduce the analyses is available in 
  this electronic appendix.

- An MPI environment is required.

- The R scripts named "AnalysisCluster.R" in the Jobscripts subfolders require the 
  RMPISNOW shell script from the R package "snow".
  Therefore, before executing these scripts you need to install the RMPISNOW shell script 
  from the installed 'snow' R package or 'inst' directory of the package sources
  of the 'snow' R package in an appropriate location, preferably
  on your path. 
  See http://homepage.divms.uiowa.edu/~luke/R/cluster/cluster.html for more details.
  Subsequently, you need to create sh files, each for a different of the
  above R scripts. The following is the content of an example sh file "simulation_clustdata.sh":

  #!/bin/bash
  #SBATCH -o /myoutfiledirectory/myjob.%j.%N.out
  #SBATCH -D /myhomedirectory
  #SBATCH -J LargeStudy
  #SBATCH --get-user-env 
  #SBATCH --clusters=myclustername
  #SBATCH --partition=mypartitionname
  #SBATCH --qos=mypartitionname
  #SBATCH --nodes=??
  #SBATCH --ntasks-per-node=??
  #SBATCH --mail-type=end
  #SBATCH --mail-user=my@mail.de
  #SBATCH --time=??:??:??

  mpirun RMPISNOW < ./multi-omics-data/Jobscripts/AnalysisCluster.R

  The above sh file of course has to be adjusted to be useable (e.g., the "?"s have
  to replaced by actual numbers, the directories have to be adjusted and
  you need to specify your e-mail address; an e-mail will be sent to this address
  once the job is finished).
