# Function that imports and wrangles fastqc data
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
  
  data_list$sequence_length_distribution <- NULL
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
