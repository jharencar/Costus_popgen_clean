# Log every R warning to a text file (and still echo to stderr for Slurm logs).
# Usage at top of a batch script (after library() calls):
#   source(file.path("QTL", "log_all_warnings.R"))  # or source("log_all_warnings.R")
#   start_warning_log("my_script_warnings.txt")
#
# Override path: WARNING_LOG=/path/to/warnings.txt Rscript ...

#' Begin logging warnings. Call once per script; registers on.exit() cleanup.
#' @param path Output file (default: env WARNING_LOG, else argument default_basename)
#' @param default_basename Used when path and WARNING_LOG are unset
start_warning_log <- function(path = NULL, default_basename = "script_warnings.txt") {
  if (!is.null(path) && !nzchar(path)) path <- NULL
  if (is.null(path)) {
    env_path <- Sys.getenv("WARNING_LOG", unset = "")
    path <- if (nzchar(env_path)) env_path else default_basename
  }
  options(warn = 1)
  con <- file(path, open = "wt")
  sink(con, type = "message", split = TRUE)
  msg <- paste0(
    "Warning log started: ", path,
    " at ", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"),
    "\n"
  )
  cat(msg)
  writeLines(msg, con)
  assign(".warning_log_con", con, envir = .GlobalEnv)
  assign(".warning_log_path", path, envir = .GlobalEnv)
  on.exit(stop_warning_log(), add = TRUE)
  message("Warning log: ", path, " (options(warn=1); all warnings recorded)")
  invisible(path)
}

#' Close warning log sink (called automatically via on.exit).
stop_warning_log <- function() {
  if (!exists(".warning_log_con", envir = .GlobalEnv, inherits = FALSE)) return(invisible(NULL))
  con <- get(".warning_log_con", envir = .GlobalEnv)
  path <- get(".warning_log_path", envir = .GlobalEnv)
  if (isOpen(con)) {
    writeLines(paste0("\nWarning log closed: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")), con)
  }
  sink(type = "message")
  if (isOpen(con)) close(con)
  rm(".warning_log_con", ".warning_log_path", envir = .GlobalEnv)
  message("Warnings saved to: ", path)
  invisible(path)
}

source_warning_log_util <- function() {
  candidates <- c("log_all_warnings.R", file.path("QTL", "log_all_warnings.R"))
  hit <- candidates[file.exists(candidates)]
  if (length(hit) == 0) {
    stop("log_all_warnings.R not found in working directory or QTL/", call. = FALSE)
  }
  source(hit[1], local = FALSE)
}
