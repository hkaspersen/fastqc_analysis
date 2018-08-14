#!/usr/bin/env Rscript
args <- commandArgs(trailingOnly = TRUE)

report_loc <- args[1]
output_dir <- args[2]

# Libraries

library(fastqcr)
library(dplyr)
library(ggplot2)
library(tidyr)
library(viridis)
library(ggsci)
library(scales)
library(svglite)
library(R.utils)

# Import functions
sourceDirectory("src/")

# Check output directory
check_dir(output_dir)
output_dir <- paste0(output_dir, "/fastqc_results_", Sys.Date())

# Data and plotting
df_list <- get_fastqc_data(report_loc)
create_plots(df_list)
