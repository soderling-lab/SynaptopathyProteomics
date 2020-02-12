# Rprofile.

# Startup.
# Note: startup messages are disable by bash alias: R --silent.
# This just prints the current version of R when launching an 
# interactive session.
.First <- function() {
	if(interactive()) {
		cat(R.version$version.string,"\n") 
	}
}

# Shutdown.
#.Last <- function() { }

# Use tar.
Sys.setenv(TAR = "/bin/tar")

# Don't prompt me to save when quiting.
utils::assignInNamespace("q", function(save = "no", status = 0, runLast = TRUE) {
				 .Internal(quit(save, status, runLast)) }, "base")    
utils::assignInNamespace("quit", function(save = "no", status = 0, runLast = TRUE) {
				 .Internal(quit(save, status, runLast)) }, "base")

# Set default CRAN mirror.
local({
	r <- getOption("repos")
	r["CRAN"] <- "http://cran.r-project.org"
	options(repos = r)
})

# Some hidden functions.
.env <- new.env()

# Create functions that can be called without parentheses.
print.command <- function (cmd) {
	  default.args <- attr(cmd, "default.args")
  if (length(default.args) == 0L) default.args <- list()
    res <- do.call(cmd, default.args, envir = parent.frame(2))
    if (attr(cmd, "print_result")) print(res)
      invisible(NULL)
}
make_command <- function(x, ..., print = TRUE) {
	class(x) <- c("command", class(x))
	attr(x, "default.args") <- list(...)
	attr(x, "print_result") <- print
	x
}

## Quit aliases.
# exit
exit <- function() { quit() }
exit <- make_command(exit,print = FALSE)
.env$exit <- exit

# q
q <- function() { quit() }
q <- make_command(q,print = FALSE)
.env$q <- q

# Add any functions.
attach(.env,warn.conflicts=FALSE)
