# Function that searches recursively for all filenames with fastqc.zip
file_names_recursive <- function(filepath) {
  files <- list.files(path = filepath, pattern = "_fastqc.zip", recursive = TRUE)
  return(files)
}
