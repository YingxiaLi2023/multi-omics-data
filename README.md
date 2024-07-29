######################################################################################

This is the electronic appendix (R code) to the article "Does combining numerous data types in multi-omics data improve or hinder performance in survival prediction? Insights from a large-scale benchmark study" by Yingxia Li1, Tobias Herold2, Ulrich Mansmann1, Roman Hornung1,3

1 Institute for Medical Information Processing, Biometry and Epidemiology, University of Munich, Marchioninistr. 15, 81377 Munich, Germany;

2 Laboratory for Leukemia Diagnostics, Department of Medicine III, LMU University Hospital, LMU Munich, Munich, Germany;

3 Munich Center for Machine Learning (MCML), Munich, Germany;

######################################################################################
## Rrogram and Platfprm



## Repository Structure

### 1. Data Subfolder
- **Purpose:** Contains scripts for downloading the OpenML data needed for reproducing the analyses.
- **Contents:**
  - `down_data.R`: Downloads data from OpenML.
  - `dataset_ids.RData`: Includes OpenML IDs for the datasets to be downloaded.

### 2. JobScripts Subfolder
- **Purpose:** Contains scripts for reproducing the benchmark study.
- **Contents:**
  - `AnalysisCluster_1_4_5.R`: For single blocks, combinations of 4 blocks, and combinations of 5 blocks.
  - `AnalysisCluster_2.R`: For combinations of 2 blocks.
  - `AnalysisCluster_3.R`: For combinations of 3 blocks.

### 3. Functions Subfolder
- **Purpose:** Contains scripts with functions used in the benchmark study.
- **Contents:**
  - Scripts whose labels contain `Functions_AnalysisCluster`: Functions for applying the different prediction methods to the different combinations.

### 4. Evaluations Subfolder
- **Purpose:** Contains scripts for evaluating results and reproducing figures and tables.
- **Contents:**
  - `Evaluation_AnalysisCluster_fivemethods.R`: Evaluates raw results.
  - `bootstrap analysis_ibrier.R` and `bootstrap analysis_cindex.R`: Performs the bootstrap analysis.
  - Scripts labeled `figures`: Reproduces figures shown in the paper and supplement, where `figures_2_S6.R` in addition produces the results of the analysis presented at the end of Section "Best-performing combinations of prediction methods and blocks per dataset".
  - `test_for_figure_2.R` and `tests_for_table_3.R`: Performs the statistical tests for Figure 2 and Table 3.

### 5. Results Subfolder
- **Purpose:** Contains results and figures from the benchmark study.
- **Contents:**
  - **rda_files Subfolder:**
    - `scenariogrid1.Rda`, `scenariogrid2.Rda`, `scenariogrid3.Rda`: Generated by the R scripts in `JobScripts`.
    - `resultsumsum.RData`, `resultsum.RData`, `CI_cindex.xlsx`, `CI_ibrier.xlsx`: Generated by R scripts in `Evaluations`.
  - **Figures Subfolder:**
    - Contains all figures from the paper and supplement. Most figures have two versions, one with the suffix "_raw" and one without. The versions with the suffix "_raw" were generated by the R code, and the versions without the suffix were subsequently edited for visual reasons.

### 6. Features_1000
-   Purpose:  as a sensitivity analysis, just change the selected featrues from 2500 into 1000.

Further details on the benchmark study are available in the paper.
