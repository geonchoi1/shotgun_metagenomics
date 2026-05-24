#!/usr/bin/env Rscript
# Step 3: CooccurrenceAffinity for ARG vs MGE per source.
# Filters: alpha (alpha_mle) > 2, together >= 10, individual >= 20.
# Input:  $PROJECT/cross/mobile_arg/step2/{plasmid,mag,unbinned}_arg_mge_matrix.tsv
# Output: $PROJECT/cross/mobile_arg/step3/{src}_affinity.tsv  (filtered + full)

suppressPackageStartupMessages({
  if (!requireNamespace("CooccurrenceAffinity", quietly = TRUE)) {
    stop("R package 'CooccurrenceAffinity' not installed in env $ENV_R")
  }
  library(CooccurrenceAffinity)
})

PROJECT <- Sys.getenv("PROJECT")
if (PROJECT == "") stop("ERROR: export PROJECT=...")

in_dir  <- file.path(PROJECT, "cross/mobile_arg/step2")
out_dir <- file.path(PROJECT, "cross/mobile_arg/step3")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

sources <- c("plasmid", "mag", "unbinned")

run_one <- function(src) {
  f <- file.path(in_dir, paste0(src, "_arg_mge_matrix.tsv"))
  if (!file.exists(f) || file.info(f)$size == 0) {
    message("  skip ", src, ": missing ", f); return(invisible())
  }
  m <- read.table(f, header = TRUE, sep = "\t", check.names = FALSE,
                  comment.char = "", quote = "", stringsAsFactors = FALSE)
  rownames(m) <- m$contig; m$contig <- NULL
  if (ncol(m) < 2 || nrow(m) < 2) {
    message("  skip ", src, ": matrix too small (", nrow(m), "x", ncol(m), ")")
    return(invisible())
  }

  # CooccurrenceAffinity expects sites (rows) x species (cols). Already in that shape.
  # affinity() returns long-form data.frame of pairwise associations.
  af <- tryCatch(
    affinity(data = as.matrix(m), row.or.col = "col"),
    error = function(e) { message("  affinity() failed for ", src, ": ", e$message); NULL }
  )
  if (is.null(af)) return(invisible())

  # Keep ARG-MGE pairs only (one side ARG:, other side IS/INT/ICE)
  is_arg <- function(x) grepl("^ARG:", x)
  is_mge <- function(x) grepl("^(IS|INT|ICE):", x)
  pair_ok <- (is_arg(af$entity_1) & is_mge(af$entity_2)) |
             (is_arg(af$entity_2) & is_mge(af$entity_1))
  af <- af[pair_ok, , drop = FALSE]

  full_out <- file.path(out_dir, paste0(src, "_affinity_full.tsv"))
  write.table(af, full_out, sep = "\t", quote = FALSE, row.names = FALSE)

  # Filter: alpha_mle > 2, X (together) >= 10, marg_1 / marg_2 >= 20
  alpha_col <- if ("alpha_mle" %in% names(af)) "alpha_mle" else
               if ("alpha" %in% names(af)) "alpha" else NA
  tog_col   <- if ("X" %in% names(af)) "X" else
               if ("together" %in% names(af)) "together" else NA
  m1_col    <- if ("mA" %in% names(af)) "mA" else
               if ("marg_1" %in% names(af)) "marg_1" else NA
  m2_col    <- if ("mB" %in% names(af)) "mB" else
               if ("marg_2" %in% names(af)) "marg_2" else NA

  if (any(is.na(c(alpha_col, tog_col, m1_col, m2_col)))) {
    message("  ", src, ": unexpected affinity() columns: ",
            paste(names(af), collapse = ","))
    return(invisible())
  }

  keep <- af[[alpha_col]] > 2 &
          af[[tog_col]]   >= 10 &
          af[[m1_col]]    >= 20 &
          af[[m2_col]]    >= 20
  af_f <- af[which(keep), , drop = FALSE]

  filt_out <- file.path(out_dir, paste0(src, "_affinity.tsv"))
  write.table(af_f, filt_out, sep = "\t", quote = FALSE, row.names = FALSE)
  message("  ", src, ": pairs full=", nrow(af), " filtered=", nrow(af_f),
          " -> ", filt_out)
}

for (s in sources) run_one(s)
message("[step3] DONE")
