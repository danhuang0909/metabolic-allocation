Multidimensional Analysis of Single-Cell Data Reveals Metabolic Resource Reallocation Patterns in the Tumor Microenvironment of Colorectal Cancer
---------------------------------------------------------------------------------------
This resource provides the R code and processed data to reproduce key results described in Dan H, et al. Multidimensional Analysis of Single-Cell Data Reveals Metabolic Resource Reallocation Patterns in the Tumor Microenvironment of Colorectal Cancer

### Getting started
**1.** Clone Github repository. 
**1.** Clone Github repository and get the data from google drive. 
```
https://github.com/danhuang0909/metabolic-allocation.git
```
download the data from the  following link and store them in the ./data folder
https://www.synapse.org/Synapse:syn64302025/files/


![GitHub图像](/metabolic-allocation/figures/figure1a.tif)


**2.** Run the following script to get the figures and supplementary figures
```
cd ./Rscript
## integration of the scRNA-seq data
step1_integration_of_scRNA_data.r
## load the data and calculate the pathway allocation values
Rscript step0_load_data2.R
## Figure 1 and supplementary Figure1
Rscript step1_figure1_final.r
# Figure 2 and supplementary Figure2-4
Rscript step1_figure2_final.r
# Figure 3 and supplementary Figure5
Rscript step1_figure3_final.r
# Figure 4 and supplementary Figure6-7
Rscript step1_figure4_final.r
# Figure 5 and supplementary Figure8
Rscript step1_figure5_final.r
# Figure 6 and supplementary Figure9
Rscript step1_figure6_final.r
# Figure 7 and supplementary Figure10-11
Rscript step1_figure7_final.r
```
### Contact
danhuang2018dana@gmail.com
