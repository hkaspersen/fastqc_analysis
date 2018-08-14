# Function that prepares the sequence length data for plotting
prepare_seq_len_data <- function(list) {
  x <- split(list$sequence_length_distribution, list$sequence_length_distribution$group)
  x <- lapply(x, filter_counts)
  x <- bind_rows(x)
  return(x)
}