# Environment ------------------------------------------------------------------
library(tidyverse)
library(here)


# Data I/O ---------------------------------------------------------------------
# data on server
# prefix_path <- here("/project", "deid_bblrepo1")
# data on local drive
prefix_path <- here("data", "raw")

health <- here(prefix_path, "n1601_dataFreeze", "health",
               "n1601_health_20170421.csv") %>% 
    read_csv(col_types = cols())

psychopathology_bifactor <- here(
        prefix_path, "n1601_dataFreeze", "clinical",
        "n1601_goassess_itemwise_bifactor_scores_20161219.csv"
    ) %>% 
    read_csv(col_types = cols())

psychopathology_corrtraits <- here(
    prefix_path, "n1601_dataFreeze", "clinical",
    "n1601_goassess_itemwise_corrtraits_scores_20161219.csv"
) %>% 
    read_csv(col_types = cols())

cognition <- here(prefix_path, "n1601_dataFreeze", "cnb",
                  "n1601_cnb_factor_scores_tymoore_20151006.csv") %>% 
    read_csv(col_types = cols())

goassess <- here(prefix_path, "n1601_dataFreeze", "clinical",
                 "n1601_diagnosis_dxpmr_20170509.csv") %>% 
    read_csv(col_types = cols())
    
t1_qa <- here(prefix_path, "n1601_dataFreeze", "neuroimaging", "t1struct", 
              "n1601_t1QaData_20170306.csv") %>% 
    read_csv(col_types = cols())

ct <- here(prefix_path, "n1601_dataFreeze", "neuroimaging",
                     "t1struct", "n1601_glasserCTMultiLabel_20161220.csv") %>% 
    read_csv(col_types = cols())

gmd <- here(prefix_path, "n1601_dataFreeze", "neuroimaging",
            "t1struct", "n1601_glasserGMDMultiLabel_20161220.csv") %>% 
    read_csv(col_types = cols())

vol <- here(prefix_path, "n1601_dataFreeze", "neuroimaging",
            "t1struct", "n1601_glasserVolMultiLabel_20161220.csv") %>% 
    read_csv(col_types = cols())

subcortical <- here(prefix_path, "n1601_dataFreeze", "neuroimaging",
                    "t1struct", "n1601_freesurferAsegVol_20180828.csv") %>% 
    read_csv(col_types = cols())
brainsummary <- here(prefix_path, "n1601_dataFreeze", "neuroimaging",
                     "t1struct", "n1601_ctVol20170412.csv") %>% 
    read_csv(col_types = cols())
demog <- here(prefix_path, "n1601_dataFreeze", "demographics",
              "n1601_demographics_go1_20161212.csv") %>% 
    read_csv(col_types = cols())

tse <- here("data", "raw", "pnc_full_data.csv") %>% 
    read_csv(col_types = cols()) %>% 
    select(bblid, starts_with("ptd"), -ptd005) %>% 
    rowwise() %>% 
    mutate(tse = sum(ptd001, ptd002, ptd003, ptd004, ptd006, ptd007, ptd008, ptd009, na.rm = TRUE)) %>% 
    mutate(assault = sum(ptd003, ptd004, ptd006, ptd009, na.rm = TRUE)) %>% 
    mutate(nonassault = sum(ptd001, ptd002, ptd007, ptd008, ptd009, na.rm = TRUE)) %>% 
    ungroup() %>% 
    mutate(nonassault = if_else(nonassault >= 1, 1, 0)) %>% 
    mutate(assault = if_else(assault >= 1, 1, 0)) %>% 
    rowwise() %>%
    mutate(trauma_type = sum(nonassault + assault * 2)) %>% 
    ungroup() %>% 
    mutate(trauma_type = if_else(trauma_type >= 2, 2, trauma_type)) %>% 
    mutate(trauma_type = factor(trauma_type, levels = c(0, 1, 2), 
                                  labels = c("No TSE", "Non-Assaultive", "Assaultive"))) %>% 
    select(-c(nonassault, assault))

ses <- here("data", "raw", "pnc_longitudinal_h2o.csv") %>% 
    read_csv(show_col_types = FALSE) %>% 
    select(bblid, envSES)


# Preprocessing ----------------------------------------------------------------
preprocessed_df <- t1_qa %>% 
    left_join(health, by = c("bblid", "scanid")) %>% 
    left_join(ses, by = c("bblid")) %>% 
    left_join(tse, by = c("bblid")) %>% 
    left_join(demog, by = c("bblid", "scanid")) %>%
    left_join(goassess, by = c("bblid", "scanid")) %>%
    left_join(psychopathology_bifactor, by = c("bblid", "scanid")) %>% 
    left_join(psychopathology_corrtraits, by = c("bblid", "scanid")) %>% 
    left_join(cognition, by = c("bblid", "scanid")) %>% 
    left_join(brainsummary, c("bblid", "scanid")) %>% 
    left_join(vol, by = c("bblid", "scanid")) %>% 
    left_join(ct, by = c("bblid", "scanid")) %>% 
    left_join(gmd, by = c("bblid", "scanid")) %>%
    left_join(subcortical, by = c("bblid", "scanid")) %>%
    filter(healthExcludev2 == 0 & t1Exclude == 0) %>%  
    mutate(age = round(ageAtScan1 / 12, 2), .after = "ageAtScan1") %>% 
    mutate(sex = factor(sex, levels = c(1, 2), labels = c("Male", "Female"))) %>% 
    mutate(race2 = factor(race2, levels = 1:3, labels = c("White", "Black", "Others"))) %>% 
    mutate(handednessv2 = factor(handednessv2)) %>% 
    mutate(is_td = if_else(goassessDxpmr4 == "1TD", 1, 0)) %>%
    select(bblid, scanid, averageManualRating,
           age, sex, race2, medu1, tse, race2, trauma_type, envSES, is_td, 
           overall_psychopathology_4factorv2,
           mood_corrtraitsv2:fear_corrtraitsv2,
           Overall_Efficiency_Ar,
           F1_Social_Cognition_Efficiency_Ar:F4_Executive_Efficiency_Ar,
           mprage_antsCT_vol_TBV,
           matches("glasser"),
           contains(c("Amygdala", "Caudate", "Thalamus.Proper", "Putamen", "Pallidum",
                      "Cerebellum.Cortex", "Hippocampus"))) %>%
    mutate(
        overall_functioning = (-overall_psychopathology_4factorv2 + Overall_Efficiency_Ar) / 2,
        .before = overall_psychopathology_4factorv2) %>% 
    drop_na()

preprocessed_df$tse_ar <- preprocessed_df %>% 
    #lm(formula = tse ~ sex * age + I(age^2)) %>% 
    lm(formula = tse ~ age * sex) %>% 
    resid()
preprocessed_df$srs <- preprocessed_df %>% 
    lm(formula = -overall_functioning ~ tse) %>% 
    resid()
#preprocessed_df$srs_ses <- preprocessed_df %>% 
#    lm(formula = -overall_functioning ~ tse + envSES + medu1) %>% 
#    resid()

write_rds(preprocessed_df, here("data", "processed", "preprocessed_data.rds"))
