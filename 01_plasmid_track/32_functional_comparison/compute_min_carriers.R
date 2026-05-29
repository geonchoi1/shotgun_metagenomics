#!/usr/bin/env Rscript
# Minimum carrier threshold for Fisher's exact test (R implementation).
# Uses R pwr::pwr.2p2n.test() — single-function, two-sided, supports unequal n.
# Asymptotically equivalent to Python statsmodels.stats.power.NormalIndPower
# (both Cohen's h + z-test).
#
# Usage:
#   Rscript compute_min_carriers.R --n1 188 --n2 300 [--OR 4] [--alpha 0.001] [--power 0.80]
#
# Defaults: OR=4 (strong env-specific effect), alpha=0.001 (BH-FDR proxy),
#           power=0.80 (Cohen 1988 convention).

# ---- CLI parsing ----
args <- commandArgs(trailingOnly = TRUE)
get_arg <- function(name, default = NULL) {
  i <- which(args == paste0("--", name))
  if (length(i) > 0 && i < length(args)) return(args[i + 1])
  if (!is.null(default)) return(default)
  stop(sprintf("Missing required argument: --%s", name))
}

n1     <- as.integer(get_arg("n1"))
n2     <- as.integer(get_arg("n2"))
OR_tgt <- as.numeric(get_arg("OR",    "4.0"))
alpha  <- as.numeric(get_arg("alpha", "0.001"))
target <- as.numeric(get_arg("power", "0.80"))

# ---- Given (k, OR, n1, n2), solve (p1, p2) such that:
#      n1*p1 + n2*p2 = k    (expected total carriers)
#      OR = [p1/(1-p1)] / [p2/(1-p2)]
make_pair <- function(k, OR, n1, n2) {
  gap <- function(p2) {
    p1 <- if (OR == 1) p2 else OR * p2 / (1 + (OR - 1) * p2)
    n1 * p1 + n2 * p2 - k
  }
  p2 <- tryCatch(uniroot(gap, c(1e-6, 0.99))$root, error = function(e) NA)
  if (is.na(p2)) return(c(NA_real_, NA_real_))
  p1 <- OR * p2 / (1 + (OR - 1) * p2)
  c(p1, p2)
}

suppressPackageStartupMessages(library(pwr))

cat(sprintf("\n=== MIN_CARRIERS auto-compute (R pwr::pwr.2p2n.test) ===\n"))
cat(sprintf("  n1            = %d\n", n1))
cat(sprintf("  n2            = %d\n", n2))
cat(sprintf("  Target OR     = %.1f\n", OR_tgt))
cat(sprintf("  Alpha         = %.0e\n", alpha))
cat(sprintf("  Target power  = %.2f\n", target))

# ---- Search minimum k ----
k_min <- NA_integer_
for (k in 2:400) {
  p <- make_pair(k, OR_tgt, n1, n2)
  if (any(is.na(p)) || p[1] >= 0.999 || p[2] >= 0.999) next
  h <- ES.h(p[1], p[2])  # Cohen's h
  pw <- tryCatch(
    pwr.2p2n.test(h = abs(h), n1 = n1, n2 = n2,
                  sig.level = alpha, alternative = "two.sided")$power,
    error = function(e) NA_real_
  )
  if (!is.na(pw) && pw >= target) {
    k_min <- k
    break
  }
}

if (is.na(k_min)) {
  cat("\n  WARN: no k<=400 achieves the target power. Relax alpha or accept lower power.\n\n")
  quit(status = 1)
}

p <- make_pair(k_min, OR_tgt, n1, n2)
pw_check <- pwr.2p2n.test(h = abs(ES.h(p[1], p[2])), n1 = n1, n2 = n2,
                          sig.level = alpha, alternative = "two.sided")$power

cat(sprintf("\n  Recommended MIN_CARRIERS = %d  (power = %.1f%%)\n", k_min, pw_check * 100))
cat(sprintf("\n  -> export MIN_CARRIERS=%d\n", k_min))
cat(sprintf("     bash run.sh\n\n"))

# Sensitivity at adjacent k
cat("  [Sensitivity at adjacent k]:\n")
for (kk in unique(c(max(2, k_min - 20), max(2, k_min - 10), k_min,
                    k_min + 10, k_min + 25, k_min + 50))) {
  pk <- make_pair(kk, OR_tgt, n1, n2)
  if (any(is.na(pk)) || pk[1] >= 0.999) next
  pkw <- pwr.2p2n.test(h = abs(ES.h(pk[1], pk[2])), n1 = n1, n2 = n2,
                       sig.level = alpha, alternative = "two.sided")$power
  mark <- if (kk == k_min) " <- rec" else ""
  cat(sprintf("    k=%3d -> power = %5.1f%%%s\n", kk, pkw * 100, mark))
}

# Cross-check options (manual)
cat("\n  Cross-check equivalents:\n")
cat(sprintf("    R base : power.prop.test(n=%d, p1=%.4f, p2=%.4f, sig.level=%.4f, alternative='two.sided')\n",
            min(n1, n2), p[1], p[2], alpha))
cat(sprintf("    R pwr  : pwr.2p2n.test(h=ES.h(%.4f,%.4f), n1=%d, n2=%d, sig.level=%.4f)\n",
            p[1], p[2], n1, n2, alpha))
cat(sprintf("    Python : statsmodels.stats.power.NormalIndPower (see compute_min_carriers.py)\n"))
cat("\n")
