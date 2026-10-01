## audit_numbers_peek.R -- are run L and run L2 files the same simulated data (same seeds)?
a <- read.csv("runL_base_asym_pvalues.csv", check.names = FALSE)
b <- read.csv("runL_base_asym_m10_pvalues.csv", check.names = FALSE)
c5 <- read.csv("runL_base_asym_m5_pvalues.csv", check.names = FALSE)
for (col in c("GiViTI", "Stk.LR", "u.sym", "HL10")) {
  k <- a$link == "cloglog" & a$n == 1000
  kb <- b$link == "cloglog" & b$n == 1000
  cat(col, " identical a vs m10:", isTRUE(all.equal(a[[col]][k], b[[col]][kb])),
      " a vs m5:", isTRUE(all.equal(a[[col]][k], c5[[col]][c5$link == "cloglog" & c5$n == 1000])),
      " cor:", cor(a[[col]][k], b[[col]][kb], use = "complete"), "\n")
}
cat("rows per cell:\n"); print(table(a$link, a$n)); print(table(b$link, b$n))
p <- read.csv("runL_s2_probit_pvalues.csv", check.names = FALSE); q <- read.csv("runL_s2_probit_m25sc_pvalues.csv", check.names = FALSE)
cat("s2 probit g.sym identical base vs m25sc:", isTRUE(all.equal(p$g.sym, q$g.sym)), " GiViTI:", isTRUE(all.equal(p$GiViTI, q$GiViTI)), "\n")
cat("names m25sc:", paste(names(q), collapse = " "), "\n")
for (f in list.files(pattern = "^runL_.*_pvalues\\.csv$")) { z <- read.csv(f, check.names = FALSE)
  cat(sprintf("%-32s rows %5d  links %-30s ns %s  B/cell %s\n", f, nrow(z), paste(unique(z$link), collapse = ","),
    paste(unique(z$n), collapse = ","), paste(unique(as.vector(table(z$link, z$n))), collapse = ","))) }
