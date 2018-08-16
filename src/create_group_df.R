# Creates a data frame with grouping information of the reports
create_group_df <- function(path) {
  folders <- folder_names(path)
  df <- lapply(folders, function(folder) get_grouping_variable(path, folder))
  df <- bind_rows(df)
  return(df)
}