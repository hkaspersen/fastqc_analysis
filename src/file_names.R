# Function that lists all the file names in the input folder
file_names <- function(filepath, folder) {
  files <- list.files(path = paste0(filepath, "/", folder), pattern = "_fastqc.zip")
  return(files)
}