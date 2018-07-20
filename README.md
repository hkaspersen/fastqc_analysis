# Information
R script for analysis of fastqc reports. The script relies on the .zip 
files created by fastqc, and use the fastqcr package 
(https://CRAN.R-project.org/package=fastqcr) for importing the data (no 
need to unpack the .zip files).

The script creates a few informative plots based on the data from fastqc 
and saves it in a new folder called "fastqc_results_TODAYS_DATE" in the 
output_dir location.

# Dependencies
The script depends on the following packages:

fastqcr

dplyr

ggplot2

tidyr

viridis

ggsci

scales

svglite

# Usage
Rscript fastqc_analysis.R zipfiles_location output_dir
