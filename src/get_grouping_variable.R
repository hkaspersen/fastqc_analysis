# Identifies the names of the files and groups them accoring to their respective folders
get_grouping_variable <- function(path, folder) {
  files <- file_names(path, folder)
  df <- data.frame(files = files, group = folder)
  df$files <- gsub("(.*?)_fastqc.zip", "\\1", df$files)
  return(df)
}