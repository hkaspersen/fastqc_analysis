# Libraries

library(fastqcr)
library(tidyverse)
library(viridis)
library(ggsci)

# Functions

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

# Data and plotting
df_list <- get_fastqc_data("D:/R-Projects/fastqc_analysis/reports/")

df_list$adapter_content %>%
  gather(key, value, -c(ref, Position, group)) %>%
  ggplot(aes(factor(Position, levels = unique(Position), ordered = TRUE),value, color = key))+
  geom_boxplot()+
  labs(x = "Position in Read",
       y = "Percent (%) Adapter Content")+
  scale_colour_jama()+
  theme_classic()+
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.4),
        legend.position = "bottom")

ggplot(df_list$per_tile_sequence_quality, aes(factor(Tile), Mean))+
  geom_boxplot()+
  labs(x = "Tile",
       y = "Mean Quality Score")+
  theme_classic()+
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.4))

df_list$per_base_sequence_content %>%
  gather(key, value, -c(ref, Base, group)) %>%
  ggplot(aes(factor(Base, levels = unique(Base), ordered = TRUE), value, color = key, group = key))+
  geom_line()+
  labs(x = "Position in Read",
       y = "Percent (%)")+
  theme_classic()+
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.4))
