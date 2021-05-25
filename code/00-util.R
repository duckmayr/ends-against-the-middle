#' Ensure the working directory contains desired subdirectories
#' 
#' For each desired subdirectory, the function checks if the subdirectory
#' already exists, and if it does not, creates it.
#' 
#' @param needed_dirs A character vector giving the desired subdirectories
#' 
#' @return Invisibly returns a character vector of the directories the function
#'     created because they were not present
prepare_directories <- function(needed_dirs = c("output", "plots")) {
    created_dirs <- character()
    for ( directory in needed_dirs ) {
        if ( dir.exists(directory) ) {
            next()
        } else {
            dir.create(directory)
            created_dirs <- c(created_dirs, directory)
        }
    }
    return(invisible(created_dirs))
}
