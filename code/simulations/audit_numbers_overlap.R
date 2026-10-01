## audit_numbers_overlap.R -- are the "independent" re-runs of the same cell really independent?
## Share of replicates whose G-independent p-values (GiViTI, Stk.LR, u.sym, HL10) coincide exactly, by rep index.
files <- setdiff(list.files(pattern = "^runL_.*_pvalues\\.csv$"), "runL_validate_base_pvalues.csv")
PL <- lapply(setNames(files, sub("^runL_(.*)_pvalues\\.csv$", "\\1", files)), function(f) read.csv(f, check.names = FALSE))
pairs <- list(c("base_asym", "base_asym_m10"), c("base_asym", "base_asym_m5"), c("base_asym_m10", "base_asym_m5"),
              c("base_cauchit", "base_cauchit_m10"), c("low_probit", "low_probit_m10"),
              c("s2_probit", "s2_probit_m25sc"), c("s2_probit", "s2_probit_m10"), c("s2_probit_m25sc", "s2_probit_m10"),
              c("s2_t4", "s2_t4_m10"))
same <- function(a, b) mean(abs(a - b) < 1e-12, na.rm = TRUE)
for (pr in pairs) { A <- PL[[pr[1]]]; B <- PL[[pr[2]]]
  for (lk in unique(B$link)) for (n in intersect(unique(A$n), unique(B$n))) {
    a <- A[A$link == lk & A$n == n, ]; b <- B[B$link == lk & B$n == n, ]; a <- a[order(a$rep), ]; b <- b[order(b$rep), ]
    if (nrow(a) != nrow(b)) next
    cat(sprintf("  %-15s vs %-16s %-7s n=%-5d identical: GiViTI %.3f  Stk.LR %.3f  u.sym %.3f  HL10 %.3f | rows %d\n", pr[1], pr[2], lk, n,
      same(a$GiViTI, b$GiViTI), same(a$Stk.LR, b$Stk.LR), same(a$u.sym, b$u.sym), same(a$HL10, b$HL10), nrow(a))) } }
## where in the sequence do identical replicates sit? (first rep of each worker chunk would be every 100th for B=2000, 20 workers)
A <- PL[["base_asym_m10"]]; B <- PL[["base_asym_m5"]]
a <- A[A$link == "cloglog" & A$n == 1000, ]; b <- B[B$link == "cloglog" & B$n == 1000, ]
idx <- which(abs(a$GiViTI - b$GiViTI) < 1e-12)
cat("\n  m10 vs m5 cloglog n=1000: identical reps", length(idx), "; first 30 indices:", head(idx, 30), "\n")
A <- PL[["base_asym"]]; a <- A[A$link == "cloglog" & A$n == 1000, ]; b <- PL[["base_asym_m10"]]; b <- b[b$link == "cloglog" & b$n == 1000, ]
idx <- which(abs(a$GiViTI - b$GiViTI) < 1e-12)
cat("  base_asym vs m10 cloglog n=1000: identical reps", length(idx), "; first 30 indices:", head(idx, 30), "\n")
