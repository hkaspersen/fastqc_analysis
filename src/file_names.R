# Function that lists all the file names in the input folder
file_names <- function(filepath) {
  files <- list.files(path = filepath, pattern = "_fastqc.zip")
  return(files)
}