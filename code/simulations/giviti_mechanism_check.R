## giviti_mechanism_check.R -- why does GiViTI collapse on a heavy-tailed symmetric link (cauchit)?
## Records, per replicate, the degree GiViTI's forward selection chose, the LR increments it tested
## (degree 3 given 2, degree 4 given 3), its statistic and p-value, under the logistic null and a
## cauchit truth, on the base design at n = 2500 (where run L measured GiViTI 0.237 vs EDGE-sym 0.586).
suppressPackageStartupMessages(library(givitiR))
set.seed(20260923)
n <- 2500; B <- as.integer(Sys.getenv("GM_B", "300"))
one <- function(link) {
  x <- runif(n, -3, 3); d <- rbinom(n, 1, 0.5); eta <- 0.6 * x + 0.5 * d
  y <- rbinom(n, 1, if (link == "cauchit") pcauchy(eta) else plogis(eta))
  fit <- suppressWarnings(glm(y ~ x + d, family = binomial()))
  e <- fitted(fit); le <- qlogis(e)
  dev <- sapply(1:4, function(k) suppressWarnings(glm(y ~ poly(le, k, raw = TRUE), family = binomial()))$deviance)
  tst <- tryCatch(suppressWarnings(givitiCalibrationTest(y, e, devel = "internal")), error = function(er) NULL)
  c(m = if (is.null(tst)) NA else length(tst$estimate) - 1, stat = if (is.null(tst)) NA else unname(tst$statistic),
    p = if (is.null(tst)) NA else tst$p.value, lr3 = dev[2] - dev[3], lr4 = dev[3] - dev[4], lr2 = dev[1] - dev[2],
    rangele = diff(range(le)))
}
for (lk in c("logit", "cauchit")) {
  M <- t(replicate(B, one(lk)))
  cat(sprintf("\n== %s, n=%d, B=%d ==\n", lk, n, B))
  cat("  selected degree m:        ", paste(names(table(M[, "m"])), table(M[, "m"]), sep = ":", collapse = "  "), "\n")
  cat(sprintf("  LR increment deg2|1 mean %.2f | deg3|2 mean %.2f, P(>3.84) %.3f | deg4|3 mean %.2f, P(>3.84) %.3f\n",
    mean(M[, "lr2"]), mean(M[, "lr3"]), mean(M[, "lr3"] > 3.84), mean(M[, "lr4"]), mean(M[, "lr4"] > 3.84)))
  cat(sprintf("  GiViTI stat mean %.2f | P(p<=0.05) %.3f | by m: %s\n", mean(M[, "stat"], na.rm = TRUE),
    mean(M[, "p"] <= 0.05, na.rm = TRUE),
    paste(sapply(sort(unique(M[, "m"])), function(k) sprintf("m=%d rej %.3f (stat %.1f)", k,
      mean(M[M[, "m"] == k, "p"] <= .05, na.rm = TRUE), mean(M[M[, "m"] == k, "stat"], na.rm = TRUE))), collapse = "; ")))
  cat(sprintf("  fixed-degree LR tests vs degree 1 (no selection): deg2 %.3f  deg3 %.3f  deg4 %.3f\n",
    mean(M[, "lr2"] > qchisq(.95, 1)), mean(M[, "lr2"] + M[, "lr3"] > qchisq(.95, 2)),
    mean(M[, "lr2"] + M[, "lr3"] + M[, "lr4"] > qchisq(.95, 3))))
  write.csv(M, sprintf("giviti_mechanism_%s.csv", lk), row.names = FALSE)
}
