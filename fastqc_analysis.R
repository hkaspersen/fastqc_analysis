#!/usr/bin/env Rscript
args <- commandArgs(trailingOnly = TRUE)

report_loc <- args[1]
output_dir <- args[2]
script.dir <- dirname(sys.frame(1)$ofile)

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
sourceDirectory(paste0(script.dir, "/src"))

# Create grouping data frame
group_df <- create_group_df(report_loc)

# Check output directory
check_dir(output_dir)
output_dir <- paste0(output_dir, "/fastqc_results_", Sys.Date())

# Data and plotting
df_list <- get_fastqc_data(report_loc)
create_plots(df_list)