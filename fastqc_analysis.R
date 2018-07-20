#!/usr/bin/env Rscript
args <- commandArgs(trailingOnly = TRUE)

report_loc <- args[1]
output_dir <- args[2]
libpath <- "/work/projects/nn9305k/lib/R"

.libPaths(c(libpath,.libPaths()))

# Libraries

library(fastqcr, lib.loc = libpath)
library(dplyr, lib.loc = libpath)
library(ggplot2, lib.loc = libpath)
library(tidyr, lib.loc = libpath)
library(viridis, lib.loc = libpath)
library(ggsci, lib.loc = libpath)
library(scales, lib.loc = libpath)
library(svglite, lib.loc = libpath)

# Functions

check_dir <- function(output_dir) {
  folder <- paste0("fastqc_results_", Sys.Date())
  dir.create(file.path(output_dir, folder), showWarnings = FALSE)
}

file_names <- function(filepath) {
  files <- list.files(path = filepath, pattern = "_fastqc")
  return(files)
}

create_groups <- function(df) {
  df <- df %>%
    mutate(group = gsub("(.*?)_R[1-2]_001", "\\1", ref))
  return(df)
}

get_fastqc_data <- function(filepath) {
  get_files <- file_names(filepath)
  
  data_list <- lapply(get_files,
                      FUN = function(file) {
                        qc_read(paste0(filepath, "/", file),
                                modules = "all",
                                verbose = FALSE)
                      })
  names(data_list) <- gsub("(.*?)_fastqc.zip", "\\1", get_files)
  data_list <- purrr::transpose(data_list)
  
  list_names <- names(data_list)
  list_numbers <- 1:length(list_names)
  for (i in list_numbers) {
    assign(list_names[i], bind_rows(data_list[[i]], .id = "ref"))
  }
  df_list <- list(summary,
                  basic_statistics,
                  per_base_sequence_quality,
                  per_tile_sequence_quality,
                  per_sequence_quality_scores,
                  per_base_sequence_content,
                  per_sequence_gc_content,
                  per_base_n_content,
                  sequence_length_distribution,
                  sequence_duplication_levels,
                  overrepresented_sequences,
                  adapter_content,
                  kmer_content,
                  total_deduplicated_percentage)
  names(df_list) <- list_names
  df_list <- lapply(df_list, create_groups)
  return(df_list)
}

filter_counts <- function(df) {
  df <- df %>%
    mutate(dupl = duplicated(Count)) %>%
    filter(dupl == FALSE)
  return(df)
}

prepare_seq_len_data <- function(list) {
  x <- split(list$sequence_length_distribution, list$sequence_length_distribution$group)
  x <- lapply(x, filter_counts)
  x <- bind_rows(x)
  return(x)
}

save_plots <- function(plot, title) {
  ggsave(paste0(output_dir, "/", title, ".svg"),
         plot,
         dpi = 100,
         device = "svg",
         units = "cm",
         height = 25,
         width = 25)
}

create_plots <- function(df_list) {
  p1 <- df_list$adapter_content %>%
    gather(key, value,-c(ref, Position, group)) %>%
    ggplot(aes(factor(
      Position, levels = unique(Position), ordered = TRUE
    ), value, color = key)) +
    geom_boxplot() +
    labs(
      x = "Position in Read",
      y = "Percent (%) Adapter Content",
      color = NULL,
      title = "Adapter content"
    ) +
    scale_colour_jama() +
    theme_classic() +
    theme(axis.text.x = element_text(
      angle = 90,
      hjust = 1,
      vjust = 0.4
    ),
    legend.position = "bottom")
  
  p2 <-
    ggplot(df_list$per_tile_sequence_quality, aes(factor(Tile), Mean)) +
    geom_boxplot() +
    labs(x = "Tile",
         y = "Mean Quality Score",
         title = "Per tile sequence quality") +
    theme_classic() +
    theme(axis.text.x = element_text(
      angle = 90,
      hjust = 1,
      vjust = 0.4
    ))
  
  p3 <- df_list$per_base_sequence_content %>%
    gather(key, value,-c(ref, Base, group)) %>%
    ggplot(aes(factor(
      Base, levels = unique(Base), ordered = TRUE
    ), value, color = key)) +
    geom_boxplot() +
    labs(
      x = "Position in Read",
      y = "Percent (%)",
      color = NULL,
      title = "Per base sequence content"
    ) +
    theme_classic() +
    scale_color_jama() +
    theme(axis.text.x = element_text(
      angle = 90,
      hjust = 1,
      vjust = 0.4
    ),
    legend.position = "bottom")
  
  p4 <- df_list$sequence_duplication_levels %>%
    gather(key, value,-c(ref, `Duplication Level`, group)) %>%
    ggplot(aes(
      factor(
        `Duplication Level`,
        levels = unique(`Duplication Level`),
        ordered = TRUE
      ),
      value,
      fill = key
    )) +
    geom_boxplot() +
    scale_fill_manual(values = c("#ef8a62",
                                 "#67a9cf")) +
    theme_classic() +
    labs(x = "Duplication Level",
         y = "Percent (%) of Sequences",
         fill = NULL,
         title = "Sequence duplication levels") +
    theme(legend.position = "bottom")
  
  p5 <- df_list %>%
    prepare_seq_len_data() %>%
    ggplot(aes(ref, Count, fill = group)) +
    geom_col(color = "black") +
    scale_y_continuous(labels = comma) +
    scale_fill_viridis(discrete = TRUE) +
    labs(x = NULL,
         y = "Total sequence length",
         title = "Sequence length per read size") +
    guides(fill = FALSE) +
    theme_classic() +
    theme(axis.text.x = element_text(
      angle = 90,
      hjust = 1,
      vjust = 0.4
    ),
    legend.position = "bottom") +
    facet_wrap( ~ Length)
  
  p6 <- df_list$per_sequence_quality_scores %>%
    ggplot(aes(factor(Quality), Count, fill = factor(Quality))) +
    geom_boxplot() +
    scale_y_continuous(labels = comma) +
    scale_fill_viridis(discrete = TRUE) +
    labs(x = "Quality",
         y = "Number of reads",
         title = "Per sequence quality scores") +
    guides(fill = FALSE) +
    theme_classic()
  
  p7 <- df_list$per_sequence_gc_content %>%
    ggplot(aes(factor(`GC Content`), Count)) +
    geom_boxplot() +
    labs(x = "GC content (%)",
         y = "Number of reads",
         title = "Per sequence GC content") +
    scale_y_continuous(labels = comma) +
    scale_x_discrete(breaks = as.character(seq(
      from = 0, to = 100, by = 10
    ))) +
    theme_classic()
  
  p8 <- df_list$per_base_n_content %>%
    ggplot(aes(factor(
      Base, levels = unique(Base), ordered = TRUE
    ), `N-Count`)) +
    geom_boxplot() +
    labs(x = "Position in read",
         title = "Per base N content") +
    theme_classic() +
    theme(axis.text.x = element_text(
      angle = 90,
      hjust = 1,
      vjust = 0.4
    ))
  save_plots(p1, "adapter_content")
  save_plots(p3, "per_base_sequence_content")
  save_plots(p4, "sequence_duplication_levels")
  save_plots(p5, "sequence_length_per_read_size")
  save_plots(p6, "per_sequence_quality_scores")
  save_plots(p7, "per_sequence_gc_content")
  save_plots(p8, "per_base_n_content")
}

# Check data
check_dir(output_dir)
output_dir <- paste0(output_dir, "/fastqc_results_", Sys.Date())

# Data and plotting
df_list <- get_fastqc_data(report_loc)
create_plots(df_list)
