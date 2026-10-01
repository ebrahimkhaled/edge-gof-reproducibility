## pkg280_verify_fix_battery.R -- run.all.gof(G = "auto") (regression F8) and the ensemble-row Note (regression F5),
## for one copy of ebrahim.gof loaded with pkgload::load_all(). Run it on the pre-fix clone and on the fixed clone.
## With a third argument "slow", also compares the slow battery at G = "auto" and at the resolved number.
## Run: Rscript pkg280_verify_fix_battery.R <package dir> <label> [slow]

args <- commandArgs(trailingOnly = TRUE)
PKG <- args[1]; LAB <- args[2]; SLOW <- length(args) >= 3 && args[3] == "slow"
options(width = 220, warn = 1)
suppressPackageStartupMessages(pkgload::load_all(PKG, export_all = FALSE, quiet = TRUE))
ns  <- asNamespace("ebrahim.gof")
sha <- system2("git", c("-C", shQuote(PKG), "rev-parse", "--short", "HEAD"), stdout = TRUE)
cat(LAB, "| ebrahim.gof", as.character(utils::packageVersion("ebrahim.gof")), "via load_all | HEAD", sha, "|",
    format(Sys.time()), "\n")
mk <- function(seed, n, c0 = 0, s = 1) {
  set.seed(seed)
  dat <- data.frame(xa = runif(n, -3, 3), db = rbinom(n, 1, 0.5))
  dat$out <- rbinom(n, 1, plogis(c0 + s * (0.6 * dat$xa + 0.5 * dat$db)))
  suppressWarnings(glm(out ~ xa + db, family = binomial(), data = dat))
}
bat <- function(...) tryCatch(as.data.frame(suppressMessages(run.all.gof(..., install = "no"))),
                              error = function(e) paste("ERROR:", conditionMessage(e)))
side <- function(a, b, la, lb) {
  if (is.character(a) || is.character(b)) { cat("  ", la, ":", if (is.character(a)) a else "ok", "|", lb, ":",
                                                if (is.character(b)) b else "ok", "\n"); return(invisible()) }
  stopifnot(identical(a$Test, b$Test))
  num <- function(x) c(x$Statistic, x$df, x$p_value)
  cat(sprintf("  rows %d | Test identical %s | NA pattern of (Statistic, df, p) identical %s | max abs diff %.3g | rows whose Note differs: %s\n",
              nrow(a), identical(a$Test, b$Test), identical(is.na(num(a)), is.na(num(b))),
              max(abs(num(a) - num(b)), na.rm = TRUE), paste(a$Test[a$Note != b$Note], collapse = ", ")))
  cat("  rows with no p-value under", la, ":", paste(a$Test[is.na(a$p_value)], collapse = ", "), "\n")
  m <- merge(a[, c("Test", "Statistic", "df", "p_value", "Note")], b[, c("Test", "Statistic", "df", "p_value", "Note")],
             by = "Test", sort = FALSE, suffixes = c(paste0(".", la), paste0(".", lb)))
  print(m, row.names = FALSE)
}

## ---- F8: G = "auto" at the top level ----
cat("\n==== F8. run.all.gof(G = 'auto') on the base design, n = 1000 (auto = 40) ====\n")
fb <- mk(101, 1000)
side(bat(fb, G = "auto", include_slow = FALSE), bat(fb, G = 40, include_slow = FALSE), "auto", "G40")
cat("\n  (y, predicted_probs, X) form:\n")
side(bat(as.numeric(fb$y), fitted(fb), X = model.matrix(fb), G = "auto", include_slow = FALSE),
     bat(as.numeric(fb$y), fitted(fb), X = model.matrix(fb), G = 40, include_slow = FALSE), "auto", "G40")
cat("\n  bad G values:\n")
for (G in list("many", "AUTO", c(10, 20), NA_character_))
  cat(sprintf("    G = %-14s -> %s\n", deparse(G), {
    r <- bat(fb, G = G, tests = c("EF", "DEF.poly3"), include_slow = FALSE)
    if (is.character(r)) r else paste(r$Test, format(r$p_value, digits = 4), r$Note, collapse = " || ") }))

## ---- F5: ensemble-row Note ----
cat("\n==== F5. ensemble rows and control, gof_demo wrong model (n = 800, auto = 32) ====\n")
gof_demo <- ebrahim.gof::gof_demo
wrong <- glm(outcome ~ age + bmi + sex + treatment, data = gof_demo, family = binomial())
plain <- bat(wrong, include_slow = FALSE)
ens <- grepl("^Ensemble", plain$Test)
showrows <- function(b, label) {
  if (is.character(b)) { cat("  ", label, ":", b, "\n"); return(invisible()) }
  e <- grepl("^Ensemble", b$Test); dd <- b$Test %in% c("DEF.poly2", "DEF.poly3", "DEF.stukel", "DEF.sym")
  pe <- b$p_value[e]; pp <- plain$p_value[ens]
  cat(sprintf("  %-58s ensemble p %s | same as plain battery: %s | ensemble Note: '%s'\n", label,
              paste(format(pe, digits = 8), collapse = " / "),
              if (identical(b$Test[e], plain$Test[ens])) isTRUE(all(abs(pe - pp) == 0)) else "rows differ",
              paste(unique(b$Note[e]), collapse = "' / '")))
  if (any(b$Note[dd] != "")) cat(sprintf("      DEF Notes: %s\n", paste(b$Test[dd], b$Note[dd], sep = ": ", collapse = " | ")))
}
showrows(plain, "no control")
ctl3 <- list(DEF.poly2 = list(weights = "score"), DEF.poly3 = list(weights = "score"), DEF.stukel = list(weights = "score"))
b3 <- bat(wrong, include_slow = FALSE, control = ctl3)
showrows(b3, "score form for the three DEF rows")
if (!is.character(b3)) {
  shown <- setNames(b3$p_value, b3$Test)[c("DEF.poly2", "DEF.poly3", "DEF.stukel")]
  cat(sprintf("      printed DEF p %s | CCT of those %.8g | def.ensemble.gof(fit) %.8g | def.ensemble.gof(fit, score) %s\n",
              paste(format(shown, digits = 7), collapse = " / "), ns$.combine_pvalues(shown, "cct"),
              def.ensemble.gof(wrong)$p_value,
              tryCatch(format(def.ensemble.gof(wrong, weights = "score")$p_value, digits = 8), error = function(e) "n/a")))
}
showrows(bat(wrong, include_slow = FALSE, control = list(DEF.poly3 = list(G = 10))), "control DEF.poly3 G = 10 (the battery G)")
showrows(bat(wrong, include_slow = FALSE, control = list(DEF.poly3 = list(weights = "unit"))), "control DEF.poly3 weights = 'unit'")
showrows(bat(wrong, include_slow = FALSE, control = list(DEF.sym = list(weights = "score"))), "control DEF.sym weights = 'score' (not a member)")
showrows(bat(wrong, include_slow = FALSE, control = list(DEF.stukel = list(G = "auto"))), "control DEF.stukel G = 'auto' (32)")
showrows(bat(wrong, include_slow = FALSE, control = list(DEF.poly2 = list(G = 20))), "control DEF.poly2 G = 20")
cat("  battery G = 'auto' (32):\n")
plain <- bat(wrong, include_slow = FALSE, G = "auto"); ens <- if (is.character(plain)) FALSE else grepl("^Ensemble", plain$Test)
showrows(plain, "G = 'auto', no control")
showrows(bat(wrong, include_slow = FALSE, G = "auto", control = list(DEF.poly3 = list(G = "auto"))), "G = 'auto', control DEF.poly3 G = 'auto'")
showrows(bat(wrong, include_slow = FALSE, G = "auto", control = list(DEF.poly3 = list(G = 32))), "G = 'auto', control DEF.poly3 G = 32")
showrows(bat(wrong, include_slow = FALSE, G = "auto", control = list(DEF.poly3 = list(G = 10))), "G = 'auto', control DEF.poly3 G = 10")

## ---- slow battery at G = "auto" ----
if (SLOW) {
  cat("\n==== slow battery, base design n = 300 (auto = 12) ====\n")
  f3 <- mk(301, 300)
  t0 <- proc.time()[["elapsed"]]
  a  <- bat(f3, G = "auto", include_slow = TRUE)
  b  <- bat(f3, G = 12, include_slow = TRUE)
  b2 <- bat(f3, G = 12, include_slow = TRUE)
  cat(sprintf("  three slow batteries in %.0f s\n", proc.time()[["elapsed"]] - t0))
  if (!is.character(a) && !is.character(b)) {
    cat("  G = 12 run twice, rows whose numbers differ between the two runs (stochastic rows):",
        paste(b$Test[!mapply(function(i) identical(c(b$Statistic[i], b$p_value[i]), c(b2$Statistic[i], b2$p_value[i])),
                             seq_len(nrow(b)))], collapse = ", "), "\n")
    side(a, b, "auto", "G12")
  }
}
cat("\ndone", format(Sys.time()), "\n")
