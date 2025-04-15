# docker run -ti --entrypoint /bin/bash -v /cis/home/ebridge2/Documents/research/graphstats_repos/causal_batch/:/causal_batch neurodata/cdcorr:0.0.1
# cd causal_batch/docs/causal_cdcorr_paper/Figure_reproduction/
require(tidyverse)
require(jsonlite)
require(np)

# Read cortical volume data and exclude total volume columns
dsk.vol.dat <- read_csv("../data/abcd/mri_y_smr_vol_dsk.csv", col_names = TRUE) %>%
  rename("participant_id"=src_subject_id) %>%
  mutate(participant_id = str_replace(participant_id, "^([A-Z]+)_([A-Z]+.*)$", "sub-\\1\\2")) %>%
  dplyr::select(-contains("total")) %>%
  left_join(read_csv("../data/abcd/mri_y_qc_incl.csv") %>%
              dplyr::select(src_subject_id, eventname, imgincl_t1w_include, imgincl_t2w_include) %>%
              rename("participant_id"=src_subject_id) %>%
              mutate(participant_id = str_replace(participant_id, "^([A-Z]+)_([A-Z]+.*)$", "sub-\\1\\2")),
            by=c("participant_id"="participant_id", "eventname"="eventname")) %>%
  ungroup() %>%
  mutate(passed.qc.anat = imgincl_t1w_include == 1 & imgincl_t2w_include == 1) %>%
  filter(eventname == "baseline_year_1_arm_1", passed.qc.anat) %>%
  dplyr::select(contains("vol"), "participant_id")

# Read DTI average FA data
dsk.fa.dat <- read_csv("../data/abcd/mri_y_dti_fa_fs_wm_dsk.csv", col_names = TRUE) %>%
  rename("participant_id"=src_subject_id) %>%
  mutate(participant_id = str_replace(participant_id, "^([A-Z]+)_([A-Z]+.*)$", "sub-\\1\\2")) %>%
  dplyr::select(-contains("total")) %>%
  left_join(read_csv("../data/abcd/mri_y_qc_incl.csv") %>%
              dplyr::select(src_subject_id, eventname, imgincl_dmri_include) %>%
              rename("participant_id"=src_subject_id) %>%
              mutate(participant_id = str_replace(participant_id, "^([A-Z]+)_([A-Z]+.*)$", "sub-\\1\\2")),
            by=c("participant_id"="participant_id", "eventname"="eventname")) %>%
  ungroup() %>%
  rename(passed.qc.dmri=imgincl_dmri_include) %>%
  mutate(passed.qc.dmri = as.logical(passed.qc.dmri)) %>%
  filter(eventname == "baseline_year_1_arm_1", passed.qc.dmri) %>%
  dplyr::select(contains("dmdtifp1"), "participant_id", -"dmdtifp1_399", -"dmdtifp1_400", -"dmdtifp1_401")

# Read T1w average white matter intensity data
dsk.t1w.wm.dat <- read_csv("../data/abcd/mri_y_smr_t1_white_dsk.csv", col_names = TRUE) %>%
  rename("participant_id"=src_subject_id) %>%
  mutate(participant_id = str_replace(participant_id, "^([A-Z]+)_([A-Z]+.*)$", "sub-\\1\\2")) %>%
  dplyr::select(-contains("total")) %>%
  left_join(read_csv("../data/abcd/mri_y_qc_incl.csv") %>%
              dplyr::select(src_subject_id, eventname, imgincl_t1w_include) %>%
              rename("participant_id"=src_subject_id) %>%
              mutate(participant_id = str_replace(participant_id, "^([A-Z]+)_([A-Z]+.*)$", "sub-\\1\\2")),
            by=c("participant_id"="participant_id", "eventname"="eventname")) %>%
  ungroup() %>%
  rename(passed.qc.t1w=imgincl_t1w_include) %>%
  mutate(passed.qc.t1w = as.logical(passed.qc.t1w)) %>%
  filter(eventname == "baseline_year_1_arm_1", passed.qc.t1w) %>%
  dplyr::select(contains("t1ww02"), "participant_id", -contains("cdk_mean"))

# Read T1w average white matter intensity data
dsk.t1w.gm.dat <- read_csv("../data/abcd/mri_y_smr_t1_gray_dsk.csv", col_names = TRUE) %>%
  rename("participant_id"=src_subject_id) %>%
  mutate(participant_id = str_replace(participant_id, "^([A-Z]+)_([A-Z]+.*)$", "sub-\\1\\2")) %>%
  dplyr::select(-contains("total")) %>%
  left_join(read_csv("../data/abcd/mri_y_qc_incl.csv") %>%
              dplyr::select(src_subject_id, eventname, imgincl_t1w_include) %>%
              rename("participant_id"=src_subject_id) %>%
              mutate(participant_id = str_replace(participant_id, "^([A-Z]+)_([A-Z]+.*)$", "sub-\\1\\2")),
            by=c("participant_id"="participant_id", "eventname"="eventname")) %>%
  ungroup() %>%
  rename(passed.qc.t1w=imgincl_t1w_include) %>%
  mutate(passed.qc.t1w = as.logical(passed.qc.t1w)) %>%
  filter(eventname == "baseline_year_1_arm_1", passed.qc.t1w) %>%
  dplyr::select(contains("t1wgray02"), "participant_id", -contains("cdk_mean"))

data <- list("dMRI"=dsk.fa.dat, "Volume"=dsk.vol.dat, "T1w.wm.intens"=dsk.t1w.wm.dat, 
             "T1w.gm.intens"=dsk.t1w.gm.dat)
y.contains.indicator <- list("dMRI"="dmdtifp1", "Volume"="vol", "T1w.wm.intens"="t1ww02", "T1w.gm.intens"="t1wgray02")
preproc.dat <- lapply(names(data), function(datname) {
  print(sprintf("Processing %s...", datname))
  dat.modality <- data[[datname]]
  
  all.dat <- demo_dat %>%
    left_join(dat.modality, by=c("participant_id"="participant_id")) %>%
    na.omit()
  
  Ys <- all.dat %>% dplyr::select(contains(y.contains.indicator[[datname]]))
  
  Xs.full <- all.dat %>% dplyr::select(male:parental_education) %>%
    dplyr::select(-c("gen_ability", "exec_func", "learn_mem"))
  
  Xs.scaled.full <- Xs.full %>%
    mutate(across(c(male, White, Black, Asian, Latinx, age, income, 
                    anesthesia_exposure, `Right_handed`, 
                    parental_education), ~as.vector(scale(.))))
  
  Xs.full <- Xs.full %>%
    mutate(male=factor(male, levels=c(0, 1), ordered=FALSE), White=factor(White, levels=c(0, 1), ordered=FALSE), 
           Black=factor(Black, levels=c(0, 1), ordered=FALSE), Asian=factor(Asian, levels=c(0, 1), ordered=FALSE),
           Latinx=factor(Latinx, levels=c(0, 1), ordered=FALSE), anesthesia_exposure=factor(anesthesia_exposure, levels=c(0, 1), ordered=FALSE),
           Right_handed=factor(Right_handed, levels=c(0, 1), ordered=FALSE), income=factor(income, levels=1:10, ordered=TRUE), 
           parental_education=factor(parental_education, levels=1:21, ordered=TRUE))
  
  bws.scott <- apply(Xs.full, 2, causalBatch:::scotts_rule)
  # bws.np <- npudensbw(dat=Xs.full, bwmethod="cv.ml", cores = parallel::detectCores() - 1)$bw
  
  Ts.full <- all.dat$site
  
  return(list(All=all.dat, Ys=Ys, Xs=Xs.full, Xs.scaled=Xs.scaled.full, Ts=Ts.full,  bws.scott=bws.scott)) #, bws.np=bws.np)
})

names(preproc.dat) <- names(data)

saveRDS(preproc.dat, "../data/abcd/preproc_abcd_data.rds")

write.csv(preproc.dat$dMRI$All, "../data/abcd/dMRI/all_data.csv", row.names = FALSE)
write.csv(preproc.dat$dMRI$Ys, "../data/abcd/dMRI/ys_data.csv", row.names = FALSE)
write.csv(preproc.dat$dMRI$Xs, "../data/abcd/dMRI/xs_data.csv", row.names = FALSE)
write.csv(preproc.dat$dMRI$Ts, "../data/abcd/dMRI/ts_data.csv", row.names = FALSE)