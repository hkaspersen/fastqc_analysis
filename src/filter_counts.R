# Function that filter out duplicated counts
filter_counts <- function(df) {
  df <- df %>%
    mutate(dupl = duplicated(Count)) %>%
    filter(dupl == FALSE)
  return(df)
}