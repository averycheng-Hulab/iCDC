#!/usr/bin/env Rscript

start_time <- Sys.time()

args <- commandArgs(trailingOnly = TRUE)

get_arg <- function(flag, default = NULL) {
  hit <- grep(paste0("^", flag, "="), args, value = TRUE)
  if (length(hit) == 0) return(default)
  sub(paste0("^", flag, "="), "", hit[[1]])
}

mode <- tolower(get_arg("--mode", "all"))
if (!mode %in% c("bulk", "scrna", "all")) {
  stop("Invalid --mode. Use --mode=bulk, --mode=scrna, or --mode=all")
}

cat("====================================\n")
cat(" iCDC demo runner\n")
cat(" Mode: ", mode, "\n", sep = "")
cat(" Working directory: ", getwd(), "\n", sep = "")
cat("====================================\n")


if (!dir.exists("bulk") || !dir.exists("scRNA") || !dir.exists("demo")) {
  stop("Please run this script from the repository root directory (must contain 'bulk/', 'scRNA/', and 'demo/').")
}


dir.create("demo_output", showWarnings = FALSE)
run_demo_script <- function(script_path, label) {
  cat("\n--- Running ", label, " demo ---\n", sep = "")
  if (!file.exists(script_path)) {
    stop("Cannot find demo script: ", script_path)
  }

  demo_env <- new.env(parent = globalenv())

  # Helpful info for scripts (optional)
  demo_env$DEMO_ROOT <- file.path(getwd(), "demo")
  demo_env$REPO_ROOT <- getwd()

  source(script_path, local = demo_env)

  cat("--- ", label, " demo finished ---\n", sep = "")
}


bulk_demo_path  <- file.path("demo", "demo_bulk.R")
scrna_demo_path <- file.path("demo", "demo_scRNA.R")

if (mode %in% c("bulk", "all")) {
  run_demo_script(bulk_demo_path, "bulk")
}

if (mode %in% c("scrna", "all")) {
  run_demo_script(scrna_demo_path, "scRNA")
}

end_time <- Sys.time()
cat("\n====================================\n")
cat("All requested demos completed.\n")
cat("Total time (sec): ", round(as.numeric(difftime(end_time, start_time, units = "secs")), 1), "\n", sep = "")
cat("Check outputs under 'demo/' (e.g., demo/bulk_output, demo/scRNA_output) or 'demo_output/'.\n")
cat("====================================\n")
