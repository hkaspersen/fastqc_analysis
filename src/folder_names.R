# Lists folders in the input folder
folder_names <- function(filepath) {
  folders <- list.files(path = filepath)
  return(folders)
}