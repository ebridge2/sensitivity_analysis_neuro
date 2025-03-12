# docker run -ti --entrypoint /bin/bash -v /cis/home/ebridge2/Documents/research/causal-repos/sensitivity_analysis_neuro/:/outputs ebridge2/sensitivity_ana:0.0.1
require(tidyverse)
require(dplyr)
require(lme4)
require(lmerTest)
require(broom.mixed)
require(boot)
require(parallel)
source("./mediation_analysis.R")

# test runs a subset of analyses for diagnostic purposes
settings <- list("test"=list(nroi=5, nadv=5, nrep=100), 
                 "full"=list(nroi=68, nadv=14, nrep=1000))
setting <- "full"
nroi <- settings[[setting]]$nroi; nrep=settings[[setting]]$nrep
nadv <- settings[[setting]]$nadv

## -------------- Data Pre-Processing --------------- ##

preproc.dat <- readRDS("../../data/abcd/preproc_abcd_data.rds")

sorted_data <- preproc.dat$Volume$All %>%
  select(participant_id, contains("smri_vol_cdk"))

id_col <- "participant_id"

old_vol_cols <- setdiff(names(sorted_data), id_col)
new_vol_cols <- old_vol_cols %>%
  str_remove("smri_vol_cdk_") %>%
  str_replace("lh$", ".LH") %>%
  str_replace("rh$", ".RH") %>%
  toupper() %>%
  paste0("Y.", .)

new_names <- c(id_col, new_vol_cols)
names(new_names) <- c(id_col, old_vol_cols)

names(sorted_data) <- new_names[names(sorted_data)]

lh_cols <- names(sorted_data)[str_detect(names(sorted_data), "\\.LH$")]
rh_cols <- names(sorted_data)[str_detect(names(sorted_data), "\\.RH$")]
other_cols <- names(sorted_data)[!str_detect(names(sorted_data), "\\.(LH|RH)$") & 
                                   names(sorted_data) != id_col]

lh_cols <- sort(lh_cols)
rh_cols <- sort(rh_cols)
other_cols <- sort(other_cols)

final_col_order <- c(id_col, lh_cols, rh_cols, other_cols)

vol.dat <- sorted_data[, final_col_order]

imputed_demographics <- readRDS("../../data/abcd/abcd_demo_dat_cleaned.rds")$Complete.Data

full.data <- imputed_demographics %>%
  filter(
    # Either Black=1 and White=0, or Black=0 and White=1
    (Black == 1 & White == 0 & Asian == 0 & Latinx == 0) |
      (Black == 0 & White == 1 & Asian == 0 & Latinx == 0)
  ) %>%
  rename_with(~ paste0("X.", .), -participant_id) %>%
  left_join(vol.dat, by="participant_id") %>%
  drop_na()

Ys <- full.data %>% select(contains("Y.", ignore.case=FALSE))

Xs <- full.data %>%
  select(contains("X.")) %>%
  mutate(Race = ifelse(X.Black == 1, "Black", "White"),
         Race = factor(Race, ordered=FALSE, levels=c("White", "Black")),
         `X.Right-handed`=factor(`X.Right-handed`, ordered=FALSE, levels=c(0, 1)),
         X.scanner_model = factor(X.scanner_model, ordered=FALSE),
         X.anesthesia_exposure = factor(X.anesthesia_exposure, ordered=FALSE, levels=c(0, 1)),
         X.male = factor(X.male, ordered=FALSE, levels=c(0, 1))) %>%
  select(-X.Black, -X.Asian, -X.Latinx, -X.White)

## --------------- Total Effects ---------------- ##

total_effect_results <- do.call(rbind, lapply(names(Ys)[1:nroi], function(roi) {
  time.st = Sys.time()
  print(roi)
  Y.vol <- data.frame(Ys[[roi]]); names(Y.vol) = roi
  data <- cbind(Xs, Y.vol)
  
  formula <- as.formula(paste0(roi, " ~ Race + X.age + X.male + X.scanner_model + X.site + (1|X.family.id)"))
  
  # Get bootstrap results
  boot_results <- robust_parameter_analysis(data, formula, "RaceBlack", n_boot = nrep, n_perm=nrep)
  time.end = Sys.time()
  
  data.frame(
    model = "Total Effect",
    roi = roi,
    estimate = boot_results$estimate,
    lower.ci = boot_results$lower.ci,
    upper.ci = boot_results$upper.ci,
    p.value = boot_results$p.value
  )
}))

## ---------------- Mediator Total Effects ------------- ##

adv_start = 13
adversity_vars <- c(
  "X.ADI", "X.fam.conflict", "X.mat.hardship", "X.trauma.hist", 
  "X.PM2.5", "X.NO2", "X.racism.state", "X.discrim.feeling", 
  "X.COI", "X.SVI", "X.School.Ach", "X.School.Impr", 
  "X.School.District.LogInc", "X.School.SpecEdu", "X.School.Gifted",
  "X.income", "X.parental_education"
)

if (setting == "test") {
  adversity_vars <- adversity_vars[1:nadv]
}

mediator_results <- do.call(rbind, lapply(adversity_vars, function(adv_var) {
  adv_results <- do.call(rbind, lapply(names(Ys)[1:nroi], function(roi) {
    Y.vol <- data.frame(Ys[[roi]]); names(Y.vol) = roi
    data <- cbind(Xs, Y.vol)
    
    formula <- as.formula(paste0(roi, " ~ ", adv_var, " + X.age + X.male + X.scanner_model + X.site + (1|X.family.id)"))
    
    boot_results <- robust_parameter_analysis(data, formula, adv_var, n_boot = 100, 100)
    
    data.frame(
      model = "Mediator Effect",
      roi = roi,
      predictor = adv_var,
      estimate = boot_results$estimate,
      lower.ci = boot_results$lower.ci,
      upper.ci = boot_results$upper.ci,
      p.value = boot_results$p.value,
      stringsAsFactors = FALSE
    )
  }))
}))