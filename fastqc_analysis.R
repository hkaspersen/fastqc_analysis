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
  
  data_list$kmer_content <- NULL
  
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
                  total_deduplicated_percentage)
  names(df_list) <- list_names
  df_list <- lapply(df_list, create_groups)
  df_list$basic_statistics <- df_list$basic_statistics %>%
    spread(Measure,Value)
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

save_plots <- function(plot, title, height, width) {
  ggsave(paste0(output_dir, "/", title, ".svg"),
         plot,
         dpi = 100,
         device = "svg",
         units = "cm",
         height = height,
         width = width)
}

create_plots <- function(df_list) {
  p1 <- df_list$adapter_content %>%
    left_join(., df_list$basic_statistics[, c("ref", "Sequence length")], by = "ref") %>%
    rename(seqlen = "Sequence length") %>%
    gather(key,
           value, -c(ref,
                     Position,
                     group,
                     seqlen)) %>%
    ggplot(aes(factor(
      Position,
      levels = unique(Position),
      ordered = TRUE
    ), value, color = key)) +
    geom_boxplot(outlier.size = 0.5) +
    labs(
      x = "Position in Read",
      y = "Percent (%) Adapter Content",
      color = NULL,
      title = "Adapter content"
    ) +
    scale_colour_jama() +
    scale_y_continuous(limits = c(0, 20)) +
    theme_classic() +
    theme(
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      legend.position = "bottom"
    ) +
    facet_wrap(~ seqlen, scales = "free")
  
  p2 <- df_list$per_base_sequence_content %>%
    left_join(., df_list$basic_statistics[, c("ref", "Sequence length")], by = "ref") %>%
    rename(seqlen = "Sequence length") %>%
    gather(key, value, -c(ref, Base, group, seqlen)) %>%
    ggplot(aes(factor(
      Base, levels = unique(Base), ordered = TRUE
    ), value, color = key)) +
    geom_boxplot(outlier.size = 0.5) +
    labs(
      x = "Position in Read",
      y = "Percent (%)",
      color = NULL,
      title = "Per base sequence content"
    ) +
    theme_classic() +
    scale_color_jama() +
    theme(
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      legend.position = "bottom"
    ) +
    facet_wrap( ~ seqlen, scales = "free")
  
  p3 <- df_list$sequence_duplication_levels %>%
    left_join(., df_list$basic_statistics[, c("ref", "Sequence length")], by = "ref") %>%
    rename(seqlen = "Sequence length") %>%
    gather(key, value, -c(ref, `Duplication Level`, group, seqlen)) %>%
    ggplot(aes(
      factor(
        `Duplication Level`,
        levels = unique(`Duplication Level`),
        ordered = TRUE
      ),
      value,
      fill = key
    )) +
    geom_boxplot(outlier.size = 0.5) +
    scale_fill_manual(values = c("#ef8a62",
                                 "#67a9cf")) +
    theme_classic() +
    labs(x = "Duplication Level",
         y = "Percent (%) of Sequences",
         fill = NULL,
         title = "Sequence duplication levels") +
    theme(legend.position = "bottom",
          axis.text.x = element_text(
            angle = 90,
            hjust = 1,
            vjust = 0.4
          )) +
    facet_wrap( ~ seqlen, scales = "free")
  
  p4 <- df_list %>%
    prepare_seq_len_data() %>%
    ggplot(aes(factor(Length), Count)) +
    geom_boxplot(outlier.size = 0.5) +
    scale_y_continuous(labels = comma) +
    scale_fill_viridis(discrete = TRUE) +
    labs(x = "Read Size",
         y = "Total sequence length",
         title = "Sequence length per read size") +
    guides(fill = FALSE) +
    theme_classic() +
    theme(axis.text.x = element_text(
      angle = 90,
      hjust = 1,
      vjust = 0.4
    ))
  
  p5 <- df_list$per_sequence_quality_scores %>%
    left_join(., df_list$basic_statistics[, c("ref", "Sequence length")], by = "ref") %>%
    rename(seqlen = "Sequence length") %>%
    ggplot(aes(factor(Quality), Count, fill = factor(Quality))) +
    geom_boxplot(outlier.size = 0.5) +
    scale_y_continuous(labels = comma) +
    scale_fill_viridis(discrete = TRUE) +
    labs(x = "Quality",
         y = "Number of reads",
         title = "Per sequence quality scores") +
    guides(fill = FALSE) +
    theme_classic() +
    theme(axis.text.x = element_text(size = 7)) +
    facet_wrap( ~ seqlen)
  
  p6 <- df_list$per_sequence_gc_content %>%
    left_join(., df_list$basic_statistics[, c("ref", "Sequence length")], by = "ref") %>%
    rename(seqlen = "Sequence length") %>%
    ggplot(aes(factor(`GC Content`), Count)) +
    geom_boxplot(outlier.size = 0.5) +
    labs(x = "GC content (%)",
         y = "Number of reads",
         title = "Per sequence GC content") +
    scale_y_continuous(labels = comma) +
    scale_x_discrete(breaks = as.character(seq(
      from = 0, to = 100, by = 10
    ))) +
    theme_classic() +
    facet_wrap( ~ seqlen, scales = "free")
  
  p7 <- df_list$per_base_n_content %>%
    left_join(., df_list$basic_statistics[, c("ref", "Sequence length")], by = "ref") %>%
    rename(seqlen = "Sequence length") %>%
    ggplot(aes(factor(
      Base, levels = unique(Base), ordered = TRUE
    ), `N-Count`)) +
    geom_boxplot(fill = "#e6e6e6",
                 outlier.size = 0.5) +
    labs(x = "Position in read",
         title = "Per base N content") +
    guides(fill = FALSE) +
    theme_classic() +
    theme(axis.text.x = element_blank(),
          axis.ticks.x = element_blank()) +
    facet_wrap( ~ seqlen, scales = "free")
  
  p8 <- df_list$total_deduplicated_percentage %>%
    select(-group) %>%
    gather(key = "ref", value = `1`) %>%
    left_join(., df_list$basic_statistics[, c("ref", "Sequence length")], by = "ref") %>%
    rename(seqlen = "Sequence length",
           perc = `1`) %>%
    mutate(perc = as.numeric(perc)) %>%
    ggplot(aes(factor(seqlen), perc)) +
    geom_boxplot(fill = "#e6e6e6",
                 outlier.size = 0.5) + 
    scale_y_continuous(limits = c(0, 100)) +
    labs(x = "Read size",
         y = "Total percentage of deduplicated reads",
         title = "Total deduplicated percentage") +
    theme_classic()
  
  # This produces a huge figure - excluded for now
  # p9 <- df_list$per_tile_sequence_quality %>%
  #   left_join(., df_list$basic_statistics[, c("ref", "Sequence length")], by = "ref") %>%
  #   rename(seqlen = "Sequence length") %>%
  #   ggplot(aes(factor(Base,
  #                     levels = unique(Base),
  #                     ordered = TRUE),
  #              factor(Tile),
  #              fill = Mean)) +
  #   geom_tile() +
  #   scale_fill_viridis(limits = c(-10, 5),
  #                      breaks = c(-10, -8, -6, -4, -2, 0, 2, 4),
  #                      guide=guide_colourbar(ticks=T,nbin=50,barheight=.5,label=T,barwidth=10)) +
  #   labs(y = "Tiles",
  #        x = "Position in read",
  #        title = "Per tile sequence quality (mean)",
  #        fill = NULL) +
  #   theme_classic() +
  #   theme(legend.position="bottom",
  #         legend.justification="center",
  #         legend.direction="horizontal",
  #         legend.text=element_text(color="grey20"),
  #         axis.text.x = element_blank(),
  #         axis.ticks.x = element_blank(),
  #         axis.text.y = element_text(size = 6)) +
  #   facet_wrap(~seqlen, scales = "free")
    
  save_plots(p1, "adapter_content", 25, 35)
  save_plots(p2, "per_base_sequence_content", 25, 35)
  save_plots(p3, "sequence_duplication_levels", 25, 30)
  save_plots(p4, "sequence_length_per_read_size", 20, 20)
  save_plots(p5, "per_sequence_quality_scores", 25, 35)
  save_plots(p6, "per_sequence_gc_content", 25, 30)
  save_plots(p7, "per_base_n_content", 25, 35)
  save_plots(p8, "total_deduplicated_percentage", 20, 20)
}

# Check output directory
check_dir(output_dir)
output_dir <- paste0(output_dir, "/fastqc_results_", Sys.Date())

# Data and plotting
df_list <- get_fastqc_data(report_loc)
create_plots(df_list)
