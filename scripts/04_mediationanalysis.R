# Environment ------------------------------------------------------------------
library(tidyverse)
library(here)
library(tidymodels)
library(ggplot2)
library(ggthemes)
library(doParallel)
library(quantreg)
library(glue)
library(ggridges)
library(ggseg)
library(ggsegGlasser)
library(lavaan)
library(flextable)
library(officer)

sect_properties <- prop_section(
    page_size = page_size(orient = "landscape",
                          width = 8.3, height = 11.7),
    type = "continuous",
    page_margins = page_mar()
)


# Data I/O ---------------------------------------------------------------------
atypical_qri_df <- 
    here("data", "processed", "preprocessed_data_atypical_withQRI.rds") %>% 
    read_rds()

centiles_score <- here("data", "raw", "PNC_centiles.csv") %>% 
    read_csv(show_col_types = FALSE) %>% 
    rename("scanid" = "ID")

outcome_labels <- atypical_qri_df %>% 
    select(overall_psychopathology_4factorv2:F4_Executive_Efficiency_Ar) %>% 
    colnames()


# Correlations between aQRI and Centile Scores ---------------------------------
atypical_qri_df %>% 
    select(scanid, aQRI) %>% 
    left_join(centiles_score, by = "scanid") %>% 
    select(-scanid, -age, -sex) %>% 
    mutate(aQRI = -1 * aQRI) %>% # aQRI < 5% is +1 so reverse the scores
    correlation::correlation(method = "spearman", p_adjust = "fdr") %>% 
    as_tibble() %>% 
    filter(Parameter1 == "aQRI") %>% 
    mutate(p = if_else(p < 0.001, "p < 0.001", as.character(round(p, 3)))) %>% 
    mutate_if(is.numeric, round, digits = 3) %>% 
    unite("95% CI", CI_low:CI_high, sep = ", ", remove = TRUE) %>% 
    select(-c(S, CI, Method, n_Obs)) %>% 
    flextable() %>% 
    bold(part = "header") %>% 
    save_as_docx(path = here("outputs", "tables", "aQRI_centiles_cor_table.docx"), 
                 pr_section = sect_properties)


# Correlation between Stress Load and Others -----------------------------------
atypical_qri_df %>% 
    select(aQRI, tse,
           overall_functioning:F4_Executive_Efficiency_Ar) %>% 
    correlation::correlation(method = "spearman", p_adjust = "fdr")

corr_res <- atypical_qri_df %>% 
    select(aQRI, tse,
           overall_functioning:F4_Executive_Efficiency_Ar) %>% 
    rename(
        "Traumatic Stress Load" = tse,
        "Overall Functioning" = overall_functioning,
        "Psychopathology (g)" = overall_psychopathology_4factorv2,
        "Mood" = mood_corrtraitsv2, 
        "Psychosis" = psychosis_corrtraitsv2,
        "Externalizing" = externalizing_corrtraitsv2,
        "Fear" = fear_corrtraitsv2,
        "Cognitive Efficiency (g)" = Overall_Efficiency_Ar,
        "Social Cognition Efficiency" = F1_Social_Cognition_Efficiency_Ar,
        "Complex Reasoning Efficiency" = F2_Complex_Reasoning_Efficiency_Ar,
        "Memory Efficiency" = F3_Memory_Efficiency_Ar,
        "Executive Efficiency" = F4_Executive_Efficiency_Ar
    ) %>% 
    psych::corr.test(method = "spearman", adjust = "fdr")
pdf(here("outputs", "figs", "aqri_spearman_corr.pdf"), height = 8, width = 8)
corrplot::corrplot(
    corr_res$r,
    col = colorRampPalette(c("#0C6291", "#FBFEF9", "#A63446"))(256),
    method = "square",
    mar = rep(0, 4),
    p.mat = corr_res$p,
    tl.col = "#1C1C1C",
    tl.srt = 45,
    type = c("upper"),
    diag = FALSE,
    rect.col = "lightgrey",
    number.digits = 2
)
dev.off()


# Mediation Analysis -----------------------------------------------------------
output_list <- list()
for (label in outcome_labels) {
    set.seed(1234)
    model <- glue(
        ' # direct effect
        {label} ~ c*tse_ar
        # mediator
        aQRI ~ a*tse_ar
        {label} ~ b*aQRI
    
        # indirect effect (a*b)
        ab := a*b
        # total effect
        total := c + (a*b)
        '
    )
    fit <- sem(model,
               data = atypical_qri_df,
               meanstructure = TRUE,
               se = "bootstrap",
               bootstrap = 1000)
    # summary(fit, fit.measure = TRUE, standardized = TRUE, ci = TRUE)
    output_list[[label]] <- standardizedSolution(fit)
}

outcome_items <- c(
    "Psychopathology (g)", "Mood", "Psychosis", "Externalizing", "Fear",
    "Cognitive Efficiency (g)", "Social Cognition Efficiency",
    "Complex Reasoning Efficiency", "Memory Efficiency",
    "Executive Efficiency"
)

lapply(output_list, slice, arg1 = c(1, 10)) %>% 
    bind_rows() %>% 
    as_tibble() %>% 
    mutate(outcome = rep(outcome_items, each = 2), .before = lhs) %>% 
    mutate(id = rep(1:2, 10), .before = lhs) %>% 
    mutate_if(is.numeric, ~round(., digits = 3)) %>% 
    unite("ci", ci.lower:ci.upper, remove = FALSE, sep = ", ") %>% 
    mutate(pvalue = if_else(pvalue < 0.001, "< 0.001", as.character(pvalue))) %>% 
    select(outcome, id, est.std, pvalue, ci) %>% 
    rename("est" = est.std) %>% 
    pivot_wider(names_from = id, values_from = est:ci, names_glue = "{id}_{.value}") %>% 
    select(outcome, `1_est`, `1_ci`, `1_pvalue`, `2_est`, `2_ci`, `2_pvalue`) %>% 
    flextable() %>%
    ftExtra::span_header() %>% 
    align(align = "center", part = "all") %>% 
    set_header_labels(
        outcome = "Outcome Variable",
        `1_est` = "Standardized Beta",
        `1_pvalue` = "p value",
        `1_ci` = "95% CI",
        `2_est` = "Standardized Beta",
        `2_pvalue` = "p value",
        `2_ci` = "95% CI",
        `1` = "Direct Effect"
    ) %>% 
    bold(part = "header") %>% 
    save_as_docx(path = here("outputs", "tables", "aQRI_meidations.docx"), 
                 pr_section = sect_properties)

# overall functioning
psych::mediate(overall_functioning ~ tse_ar + (aQRI),
               data = atypical_qri_df,
               std = TRUE,
               n.iter = 1000) %>%
    print(short = FALSE)
set.seed(1234)
model <- glue(
    ' # direct effect
    overall_functioning ~ c*tse_ar
    # mediator
    aQRI ~ a*tse_ar
    overall_functioning ~ b*aQRI
    
    # indirect effect (a*b)
    ab := a*b
    # total effect
    total := c + (a*b)
    '
)
fit <- sem(model,
           data = atypical_qri_df,
           meanstructure = TRUE,
           se = "bootstrap",
           bootstrap = 1000)
summary(fit, fit.measure = TRUE, standardized = TRUE, ci = TRUE)
standardizedSolution(fit)

# check if sex effect on aQRI exists -------------------------------------------
wilcox.test(formula = aQRI ~ sex, data = atypical_qri_df)
atypical_qri_df %>% 
    mutate(tse_ar_rank = factor(ntile(tse_ar, 3))) %>% 
    ggplot(aes(x = sex, y = aQRI)) +
        facet_wrap(~tse_ar_rank) +
        geom_boxplot() +
        stat_compare_means() +
        theme_pander() +
        theme(plot.margin = margin(2, 2, 2, 2, "mm"))


# QRI and Trauma ---------------------------------------------------------------
lm(formula = aQRI ~ trauma_type + sex + age + I(age^2), data = atypical_qri_df) %>% 
    report::report()

traumatype_qri_boxplot <- atypical_qri_df %>%
    ggplot(aes(x = trauma_type, y = aQRI)) +
        gg.layers::geom_boxplot2(width = 0.5, width.errorbar = 0.2) +
        stat_compare_means(
            aes(label=..p.adj..),
            comparisons = list(
                c("No TSE", "Non-Assaultive"),
                c("Non-Assaultive", "Assaultive"),
                c("No TSE", "Assaultive")
            ),
            label.y = c(0.07, 0.09, 0.11),
            label = "p.signif"
        ) +
        labs(x = "", y = "Averaged Quantile Regression Index (aQRI)") +
        theme_pander() +
        theme(plot.margin = margin(2, 2, 2, 2, "mm"))
traumatype_qri_boxplot 
ggsave(filename = here("outputs", "figs", "traumatype_qri_boxplot.pdf"),
       plot = traumatype_qri_boxplot, width = 6, height = 4)



fit <- atypical_qri_df %>% 
    lm(formula = overall_functioning ~ tse_ar + aQRI + sex * age + I(age^2) + envSES)

fit %>% 
    lm.beta::lm.beta() %>% 
    broom::tidy(conf.int = TRUE)
broom::glance(fit)
report::report(fit)
#sjPlot::plot_model(fit, type = "pred", terms = c("tse_ar", "aQRI [0, 0.1]"))

