# Environment ------------------------------------------------------------------
library(tidyverse)
library(here)
library(quantreg)
library(glue)


# Data I/O ---------------------------------------------------------------------
preprocessed_df <- here("data", "processed", "preprocessed_data.rds") %>% 
    read_rds()

normative_df <- filter(preprocessed_df, is_td == 1 & tse == 0)
atypical_df <- filter(preprocessed_df, !c(is_td == 1 & tse == 0))

# extract glasser labels
glasser_labels <- preprocessed_df %>% 
    select(contains("glasser")) %>% 
    colnames()

subcortical_labels <- preprocessed_df %>% 
    select(contains(c("Amygdala", "Caudate", "Thalamus.Proper", "Putamen", 
                      "Pallidum", "Hippocampus"))) %>% 
    colnames() 

all_labels <- c(glasser_labels , subcortical_labels)


# QRI Calculation --------------------------------------------------------------
qri_score <- select(atypical_df, bblid, scanid)

for (label in all_labels) {
    
    region_formula <- glue("`{label}` ~ sex * age + I(age^2)") %>% 
        as.formula()
    
    fit_005 <- normative_df %>% 
        rq(formula = region_formula, tau = 0.05)
    
    fit_095 <- normative_df %>% 
        rq(formula = region_formula, tau = 0.95)
    
    
    pred_005_ci <- stats::predict(
        fit_005, newdata = atypical_df, interval = "confidence"
    ) %>% as_tibble()
    pred_095_ci <- stats::predict(
        fit_095, newdata = atypical_df, interval = "confidence"
    ) %>% as_tibble()
    
    score <- atypical_df %>%
        select(all_of(label)) %>% 
        bind_cols(pred_005_ci %>% select(lower)) %>% 
        bind_cols(pred_095_ci %>% select(higher)) %>% 
        mutate(score = case_when(
            !!sym(label) < lower ~ 1,
            !!sym(label)  > higher ~ -1,
            TRUE ~ 0
            )
        ) %>% 
        select(score) %>% 
        rename(!!quo_name(label) := score)
    
    qri_score <- qri_score %>% 
        bind_cols(score)
}

atypical_qri_df <- qri_score %>% 
    rowwise() %>% 
    mutate(aQRI = mean(c_across(all_of(all_labels)), na.rm = TRUE)) %>% 
    ungroup() %>% 
    right_join(
        atypical_df %>% select(-all_of(all_labels)) , by = c("bblid", "scanid")
    )

# save data for next scripts
write_rds(
    atypical_qri_df,
    here("data", "processed", "preprocessed_data_atypical_withQRI.rds")
)
