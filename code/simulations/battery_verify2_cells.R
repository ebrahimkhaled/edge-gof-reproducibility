## battery_verify2_cells.R -- independent pre-launch check of the cell table, launch plan and seeds (no cell is run).
## Output: battery/_review/v2/cells_check.log
## ---- archive paths (inserted by make_archive.py; the convention is in README.md, "How to run") -------------
EDGE_ARCHIVE_ROOT <- local({
  r <- Sys.getenv("EDGE_ARCHIVE_ROOT")
  if (!nzchar(r)) {
    f <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))
    r <- if (length(f)) file.path(dirname(normalizePath(f[1], winslash = "/")), "..", "..") else getwd()
  }
  r <- normalizePath(r, winslash = "/", mustWork = FALSE)
  Sys.setenv(EDGE_ARCHIVE_ROOT = r)        # so worker processes started from here resolve the same root
  r
})
edge_path <- function(...) file.path(EDGE_ARCHIVE_ROOT, ...)
edge_battery <- function(...) {           # the author's battery/ folder: analysis/ -> results/analysis, rest -> results/blocks
  if (...length() == 0L) return(edge_path("results", "blocks"))
  p <- file.path(...)
  ifelse(p == "analysis" | startsWith(p, "analysis/"), edge_path("results", p), edge_path("results", "blocks", p))
}
edge_out <- function(...) { d <- edge_path("output", ...); dir.create(d, showWarnings = FALSE, recursive = TRUE); d }
## ---------------------------------------------------------------------------------------------------------------

SIMDIR <- edge_path("code/simulations")
OUT <- edge_battery("_review", "v2")
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
sink(file.path(OUT, "cells_check.log"), split = TRUE)
suppressPackageStartupMessages(library(data.table))
source(file.path(SIMDIR, "_battery_cells.R"))
C <- battery_cells()

cat("== the five cells of check (1)\n")
want <- data.frame(block = c("1b", "1b", "2", "3", "4"),
                   cell = c("sparse49_n100", "sparse49_n200", "probit_auc_n3500", "stk_short_n1000", "sparse49_cloglog_n300"))
for (i in seq_len(nrow(want))) {
  z <- C[C$block == want$block[i] & C$cell == want$cell[i], ]
  cat(sprintf("%-3s %-24s gen %-7s link %-8s n %5d B %5d seed_family %-5s cell_id %4d seed_base %.0f G_list %s null %s/%s\n",
              z$block, z$cell, z$generator, z$link, z$n, z$B, z$seed_family, z$cell_id, z$seed_base, z$G_list,
              z$null_block, z$null_cell))
}

cat("\n== counts per block (cells, sum B)\n")
tb <- as.data.table(C)[, .(cells = .N, rows = sum(B)), by = block]
print(tb)
cat(sprintf("all blocks: %d cells, %d rows; outside block 0: %d cells, %d rows\n", nrow(C), sum(C$B),
            sum(C$block != "0"), sum(C$B[C$block != "0"])))
rl <- C[C$generator == "real", ]
cat(sprintf("real one-off cells: %d, orders %s, extra rows beyond one per fit: %d\n", nrow(rl),
            paste(rl$orders, collapse = ","), sum(rl$B) - nrow(rl)))
cat(sprintf("outside block 0 without the row orders: %d rows\n", sum(C$B[C$block != "0"]) - (sum(rl$B) - nrow(rl))))

cat("\n== cells_expected()\n")
ce <- cells_expected(C)
print(ce[!ce$ok, ])
cat(sprintf("%d checks, %d not ok\n", nrow(ce), sum(!ce$ok)))

cat("\n== seed overlaps (E8.7)\n")
so <- seed_overlaps(C)
cat(sprintf("pairs %d; kinds: %s; offsets: %s\n", nrow(so), paste(names(table(so$kind)), table(so$kind), collapse = "; "),
            paste(sort(unique(so$offset)), collapse = ", ")))
cat("pairs involving a new cell:", sum(so$family_a == "new" | so$family_b == "new"), "\n")
## the row-order seeds 20260915 ... 20261914 against every cell seed range
z <- C[!is.na(C$seed_base), ]
hit <- z[z$seed_base + 1 <= 20261914 & z$seed_base + z$B >= 20260915, c("block", "cell", "seed_base", "B")]
cat("cells whose seed range meets 20260915-20261914:", nrow(hit), "\n"); if (nrow(hit)) print(hit)

cat("\n== new seed bases are disjoint from each other and within integer range\n")
nw <- C[C$seed_family == "new", ]
cat(sprintf("new cells %d, max seed %.0f (< %.0f), duplicated bases %d\n", nrow(nw), max(nw$seed_base + nw$B),
            .Machine$integer.max, anyDuplicated(nw$seed_base)))

cat("\n== launch plan against E4\n")
P <- launch_plan(C)
print(P[, c("step", "action", "block", "note")], row.names = FALSE)
b2null <- sort(unique(C$null_cell[C$block == "2"]))
s2 <- sort(strsplit(P$cells[2], ",")[[1]])
cat("step 2 cells = the null cells of every block 2 cell:", identical(s2, b2null), "; all in 1b:",
    all(s2 %in% C$cell[C$block == "1b"]), "; count", length(s2), "\n")
s6 <- strsplit(P$cells[6], ",")[[1]]
cat("step 6 cells = every block 3 cell at n = 1000:", setequal(s6, C$cell[C$block == "3" & C$n %in% 1000L]),
    "; count", length(s6), "; any n != 1000:", any(C$n[match(paste("3", s6), paste(C$block, C$cell))] != 1000L), "\n")
e4 <- c("0", "1b", "2", "1a", "1b", "3", "3", "4", "5", "6", "7")
cat("run steps block order = 0; 1b (nulls of 2); 2; 1a; 1b (rest); 3 (n = 1000); 3 (rest); 4; 5; 6; 7:",
    identical(P$block[P$action == "run"], e4), "\n")
cat("summary steps:", paste(P$block[P$action == "summary"], collapse = ", "), "\n")
cat("block 8 (E9) in the plan:", any(P$block == "8"), "\n")

cat("\n== battery/cells.csv against battery_cells() (E8.8: 'the null named in battery/cells.csv')\n")
CC <- fread(edge_battery("cells.csv"), colClasses = "character")
cat(sprintf("cells.csv rows %d, battery_cells() rows %d; columns only in the table: %s; only in the file: %s\n", nrow(CC), nrow(C),
            paste(setdiff(names(C), names(CC)), collapse = ","), paste(setdiff(names(CC), names(C)), collapse = ",")))
k1 <- paste(CC$block, CC$cell); k2 <- paste(C$block, C$cell)
cat("same cell keys:", setequal(k1, k2), "\n")
m <- match(k2, k1)
nb <- ifelse(is.na(C$null_block), "", C$null_block); nc <- ifelse(is.na(C$null_cell), "", C$null_cell)
cat("null_block equal:", all(nb == CC$null_block[m]), "; null_cell equal:", all(nc == CC$null_cell[m]), "\n")
dB <- which(as.character(C$B) != CC$B[m])
cat("cells whose B differs between cells.csv and battery_cells():", length(dB), "\n")
if (length(dB)) print(data.frame(cell = k2[dB], file_B = CC$B[m][dB], table_B = C$B[dB]))
sb <- which(format(C$seed_base, scientific = FALSE, trim = TRUE) != CC$seed_base[m] & !(is.na(C$seed_base) & CC$seed_base[m] %in% c("", "NA")))
cat("cells whose seed_base text differs:", length(sb), "\n")
if (length(sb)) print(head(data.frame(cell = k2[sb], file = CC$seed_base[m][sb], table = C$seed_base[sb]), 10))
sink()
