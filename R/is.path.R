is.path <- function(file_name) {
	# Checks if provided filename contains a valid path.
	is_path <- grepl("/",file_name)
	if (is_path) {
		path_exists <- dir.exists(dirname(file_name))
		if (!path_exists) { 
			stop("Path to file does not exist.") 
		}
	}
	return(is_path)
}

