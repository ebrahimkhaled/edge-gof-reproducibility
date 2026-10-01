## _cohort_p3.R -- the Diabetes 130-US cohort as paper 3 uses it, shared by run_I_bigdata_p3.R (Section 7) and
## run_M_blockR_realdata.R when RR_COHORT=p3 (the corrupted-record experiment on the cohort).
##
## Paper 2 split the 101,763 ADMISSIONS at random. A referee pointed out that patients recur (patient_nbr), so the
## same patients sat in both halves and the outcomes within a half were clustered, against the independence the
## tests assume. Paper 3 follows the convention of the data set's own authors (Strack et al. 2014): one admission
## per patient, the first by encounter_id, and no admission that ended in death or discharge to hospice, since
## those patients cannot be readmitted. The random split is then over patients. The model formula is unchanged.
cohort_p3 <- function(simdir) {
  D <- read.csv(file.path(simdir, "..", "data_large", "diabetic_data.csv"),
                stringsAsFactors = FALSE, na.strings = c("?", ""))
  n_raw <- nrow(D)
  D <- D[D$gender %in% c("Female", "Male"), ]
  D <- D[!(D$discharge_disposition_id %in% c(11, 13, 14, 19, 20, 21)), ]   # expired or hospice
  D <- D[order(D$patient_nbr, D$encounter_id), ]
  D <- D[!duplicated(D$patient_nbr), ]                                       # first admission per patient
  D$y <- as.integer(D$readmitted == "<30")
  age_mid <- c("[0-10)" = 5, "[10-20)" = 15, "[20-30)" = 25, "[30-40)" = 35, "[40-50)" = 45, "[50-60)" = 55,
               "[60-70)" = 65, "[70-80)" = 75, "[80-90)" = 85, "[90-100)" = 95)
  D$age_n <- unname(age_mid[D$age])
  D$female <- as.integer(D$gender == "Female")
  D$insulin <- factor(D$insulin, levels = c("No", "Steady", "Up", "Down"))
  D$change <- as.integer(D$change == "Ch"); D$diabetesMed <- as.integer(D$diabetesMed == "Yes")
  f <- y ~ age_n + female + time_in_hospital + num_lab_procedures + num_procedures + num_medications +
           number_outpatient + number_emergency + number_inpatient + number_diagnoses + insulin +
           change + diabetesMed
  D <- D[stats::complete.cases(D[, all.vars(f)]), ]
  set.seed(20260930)
  dev <- sample(nrow(D), floor(nrow(D) / 2))                 # one row per patient, so this splits patients
  stopifnot(!any(D$patient_nbr[dev] %in% D$patient_nbr[-dev]))
  list(Ddev = D[dev, ], Dval = D[-dev, ], f = f, n_raw = n_raw, n = nrow(D))
}
