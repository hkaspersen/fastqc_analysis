# Function that creates groups from the isolate names
create_groups <- function(df) {
  df <- df %>%
    mutate(group = gsub("(.*?)_R[1-2]_001", "\\1", ref))
  return(df)
}
