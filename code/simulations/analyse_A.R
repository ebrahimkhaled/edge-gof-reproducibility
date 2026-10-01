## analyse_A.R -- read runA and answer the two referee points in plain numbers.
S <- read.csv("runA_pstar_giviti_summary.csv", stringsAsFactors = FALSE)
P <- read.csv("runA_pstar_giviti_pvalues.csv", stringsAsFactors = FALSE)
B <- unique(S$B)[1]
mc <- function(r) sqrt(r * (1 - r) / B)

cat("=========================================================================\n")
cat(" POINT 3 -- SIZE as p* grows.  Omega subtracts a rank-p* matrix, so this\n")
cat(" is where the construction is supposed to break.  nominal 0.05, MCSE",
    sprintf("%.4f", mc(0.05)), "\n")
cat("=========================================================================\n")
for (G in sort(unique(S$G))) {
  cat("\n-- G =", G, ifelse(G == 10, "(G < p* at p*=20)", "(G = p* at p*=20)"), "--\n")
  z <- S[S$G == G & S$scenario == "cloglog", c("pstar", "test", "size_raw")]
  w <- reshape(z, idvar = "test", timevar = "pstar", direction = "wide")
  names(w) <- sub("size_raw.", "p*=", names(w), fixed = TRUE)
  print(w, row.names = FALSE, digits = 3)
}

cat("\n=========================================================================\n")
cat(" POINT 3 -- SIZE-ADJUSTED POWER as p* grows (signal held fixed;\n")
cat(" the added covariates are pure noise, so any change IS the p* effect)\n")
cat("=========================================================================\n")
for (sc in c("cloglog", "quad", "binint")) {
  for (G in sort(unique(S$G))) {
    cat("\n--", sc, " G =", G, "--\n")
    z <- S[S$G == G & S$scenario == sc, c("pstar", "test", "power_adj")]
    w <- reshape(z, idvar = "test", timevar = "pstar", direction = "wide")
    names(w) <- sub("power_adj.", "p*=", names(w), fixed = TRUE)
    print(w, row.names = FALSE, digits = 3)
  }
}

cat("\n=========================================================================\n")
cat(" POINT 2 -- GiViTI, on identical samples.  Does the rival the paper cites\n")
cat(" 10x as its motivation actually beat EDGE?\n")
cat("=========================================================================\n")
for (sc in c("cloglog", "quad", "binint")) {
  z <- S[S$scenario == sc & S$G == 10 & S$test %in% c("EDGE.poly3", "GiViTI", "HL", "HL_F", "Stukel"),
         c("pstar", "test", "size_raw", "power_adj", "na_alt")]
  cat("\n--", sc, "(G=10) --\n"); print(z, row.names = FALSE, digits = 3)
}

cat("\n=========================================================================\n")
cat(" DECLINE RATES -- a test that does not return a value is not a comparator\n")
cat("=========================================================================\n")
z <- aggregate(cbind(na_null, na_alt) ~ test + pstar + G, data = S, FUN = mean)
print(z[z$na_null > 0 | z$na_alt > 0, ], row.names = FALSE, digits = 3)
if (!any(z$na_null > 0 | z$na_alt > 0)) cat("  none: every test returned a value in every cell\n")

cat("\n=========================================================================\n")
cat(" HEADLINE for the paper\n")
cat("=========================================================================\n")
e <- S[S$test == "EDGE.poly3" & S$scenario == "cloglog", c("G", "pstar", "size_raw", "power_adj")]
cat("EDGE-poly3 on cloglog:\n"); print(e, row.names = FALSE, digits = 3)
worst <- e[which.max(abs(e$size_raw - 0.05)), ]
cat(sprintf("\nworst size distortion for EDGE-poly3: %.3f at G=%d p*=%d (%.1f MCSE from nominal)\n",
            worst$size_raw, worst$G, worst$pstar, abs(worst$size_raw - 0.05) / mc(0.05)))
