# Function that creates groups from the isolate names
create_groups <- function(df, groupdf) {
  df <- df %>%
    left_join(groupdf, )
  return(df)
}
