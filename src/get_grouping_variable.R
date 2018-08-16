# Identifies the names of the files and groups them accoring to their respective folders
get_grouping_variable <- function(path, folder) {
  files <- file_names(path, folder)
  df <- data.frame(files = files, group = folder, stringsAsFactors = FALSE)
  df$files <- gsub("(.*?)_(.*?)_fastqc.zip", "\\1", df$files)
  colnames(df) <- c("ref","group")
  return(df)
}
