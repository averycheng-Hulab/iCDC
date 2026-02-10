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


if (!dir.exists("demo")) stop("Missing demo/; please run from repository root.")

dir.create("demo_output", showWarnings = FALSE)
run_demo_script <- function(script_path, label) {
  cat("
--- Running ", label, " demo ---
", sep = "")
  if (!file.exists(script_path)) {
    stop("Cannot find demo script: ", script_path)
  }

  demo_env <- new.env(parent = globalenv())

  demo_env$DEMO_ROOT <- file.path(getwd(), "demo")
  demo_env$REPO_ROOT <- getwd()

  demo_start <- Sys.time()
  source(script_path, local = demo_env)
  demo_end <- Sys.time()

  elapsed_sec <- as.numeric(difftime(demo_end, demo_start, units = "secs"))
  cat("--- ", label, " demo finished ---
", sep = "")
  cat("    Elapsed time: ", sprintf("%.2f", elapsed_sec), " sec (", 
      sprintf("%.2f", elapsed_sec/60), " min)
", sep = "")
}


bulk_demo_path  <- file.path("demo", "demo_bulk.R")
scrna_demo_path <- file.path("demo", "demo_scRNA.R")

if (mode %in% c("bulk", "all")) {
  run_demo_script(bulk_demo_path, "bulk")
}

if (mode %in% c("scrna", "all")) {
  run_demo_script(scrna_demo_path, "scRNA")
}

cat("====================================\n")
cat("All requested demos completed.\n")
cat("Check outputs under 'demo/bulk_output' or 'demo/scRNA_output/'.\n")
cat("====================================\n")

end_time <- Sys.time()
total_sec <- as.numeric(difftime(end_time, start_time, units = "secs"))
cat("Total elapsed time: ", sprintf("%.2f", total_sec), " sec (", sprintf("%.2f", total_sec/60), " min)\n", sep = "")
cat("Finished at: ", format(end_time, "%Y-%m-%d %H:%M:%S"), "\n", sep = "")
