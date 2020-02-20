cytoscape_ping <- function() {
	suppressPackageStartupMessages({
		library(RCy3)
		devtools::load_all()
	})
	tryCatch(
	 expr = {
		 cytoscapePing()
	 },
	 error = function(e) {
		 prompt <- paste("Please open Cytoscape.",
				 "When ready to proceed,",
				 "press [enter] to continue...\n")
		 pause(prompt)
	 },
	 finally = {
		 cytoscapePing()
	 }
	 )
}
