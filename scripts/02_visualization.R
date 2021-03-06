# Environment ------------------------------------------------------------------
library(tidyverse)
library(here)
library(ggpubr)
library(ggthemes)
library(gt)
library(gtsummary)


# Data I/O ---------------------------------------------------------------------
preprocessed_df <- here("data", "processed", "preprocessed_data.rds") %>% 
    read_rds()

normativel_df <- filter(preprocessed_df, is_td == 1 & tse == 0)
atypical_df <- filter(preprocessed_df, !c(is_td == 1 & tse == 0))


# Distribution -----------------------------------------------------------------
tse_distplot <- atypical_df %>% 
    ggplot(aes(x = tse)) +
     geom_histogram(binwidth = 1) +
        labs(x = "Traumatic Stress Events") +
        theme_pander()

tsl_distplot <- atypical_df %>% 
    ggplot(aes(x = tse_ar)) +
        geom_histogram(binwidth = 0.5) +
        labs(x = "Traumatic Stress Load") +
        theme_pander()

ggarrange(tse_distplot, tsl_distplot, nrow = 2) %>% 
    ggexport(filename = here("outputs", "figs", "tse_distributions.pdf"),
             width = 6, height = 8)


# Table ------------------------------------------------------------------------
tse_comparison_table <- atypical_df %>% 
    mutate(tse_ar_rank = ntile(tse_ar, 3)) %>% 
    bind_rows(normativel_df %>% mutate(tse_ar_rank = 0)) %>% 
    mutate(tse_ar_rank = factor(tse_ar_rank, levels = c(0, 1, 2, 3),
                                labels = c("Normative", "Low", "Moderate", "High"))) %>% 
    mutate(mprage_antsCT_vol_TBV = as.numeric(scale(mprage_antsCT_vol_TBV))) %>% 
    select(tse_ar_rank, age:mprage_antsCT_vol_TBV) %>% 
    select(-c(is_td, tse, trauma_type)) %>% 
    tbl_summary(
        by = "tse_ar_rank",
        statistic = list(all_continuous() ~ "{mean} ({sd})"),
        label = list(
            age ~ "Age",
            sex ~ "Sex",
            medu1 ~ "Maternal Education",
            envSES ~ "Environmental SES",
            overall_functioning ~ "Overall Functioning",
            overall_psychopathology_4factorv2 ~ "Psychopathology (g)",
            mood_corrtraitsv2 ~ "Mood",
            psychosis_corrtraitsv2 ~ "Psychosis",
            externalizing_corrtraitsv2 ~ "Externalizing",
            fear_corrtraitsv2 ~ "Fear",
            Overall_Efficiency_Ar ~ "Cognitive Efficiency (g)",
            F1_Social_Cognition_Efficiency_Ar ~ "Social Cognition Efficiency",
            F2_Complex_Reasoning_Efficiency_Ar ~ "Complex Reasoning Efficiency",
            F3_Memory_Efficiency_Ar ~ "Memory Efficiency",
            F4_Executive_Efficiency_Ar ~ "Executive Efficiency",
            mprage_antsCT_vol_TBV ~ "Total Brain Volume (Z)"
        )
    ) %>% 
    add_p() %>% 
    add_q(method = "fdr") %>% 
    bold_p(q = TRUE) %>% 
    modify_spanning_header(
        c("stat_2", "stat_3", "stat_4") ~ "**Traumatic Stress Load by Rank**"
    )
tse_comparison_table

tse_comparison_table %>% 
    as_flex_table() %>% 
    flextable::bold(part = "header") %>% 
    flextable::save_as_docx(
        path = here("outputs", "tables", "tse_comparison.docx")
    )

tse_comparison_table %>% 
    as_gt() %>% 
    gtsave(
        filename = "tse_comparison.html", 
        path = here("outputs", "tables")
    )


# Visualization ----------------------------------------------------------------
# control for sex and age effect
atypical_df %>% 
    select(tse_ar, overall_psychopathology_4factorv2,
           Overall_Efficiency_Ar, mprage_antsCT_vol_TBV) %>% 
    correlation::correlation(method = "spearman", p_adjust = "fdr")

tse_psychopathology_scatterplot <- atypical_df %>% 
    ggplot(aes(x = tse_ar, y = overall_psychopathology_4factorv2)) +
        geom_point(color = "grey30", alpha = 0.9, size = 3) +
        geom_smooth(method = "lm", color = "tomato3", fill = "grey80") +
        annotate("text", x = 4, y = -1.9, label = "italic(r)==0.24*','~italic(p)<0.001", size = 4.5, parse = TRUE) +
        labs(x = "Traumatic Stress Load",
             y = "Psychopathology (g)") +
        theme_pander() +
        theme(legend.position = "none", plot.margin = margin(2, 2, 2, 2, "mm"))
tse_cognition_scatterplot <- atypical_df %>% 
    ggplot(aes(x = tse_ar, y = Overall_Efficiency_Ar)) +
        geom_point(color = "grey30", alpha = 0.9, size = 3) +
        geom_smooth(method = "lm", color = "tomato3", fill = "grey80") +
        annotate("text", x = 4, y = -7, label = "italic(r)==-0.09*','~italic(p)==0.004", size = 4.5, parse = TRUE) +
        labs(x = "Traumatic Stress Load",
             y = "Cognitive Efficiency (g)") +
        theme_pander() +
        theme(legend.position = "none", plot.margin = margin(2, 2, 2, 2, "mm"))
tse_tbv_scatterplot <- atypical_df %>% 
    mutate(mprage_antsCT_vol_TBV = as.numeric(scale(mprage_antsCT_vol_TBV))) %>% 
    ggplot(aes(x = tse_ar, y = mprage_antsCT_vol_TBV)) +
        geom_point(color = "grey30", alpha = 0.8, size = 3) +
        geom_smooth(method = "lm", color = "tomato3", fill = "grey80") +
        annotate("text", x = 4, y = -2.8, label = "italic(r)==-0.09*','~italic(p)==0.004", size = 4.5, parse = TRUE) +
        labs(x = "Traumatic Stress Load",
             y = "Total Brain Volume") +
        theme_pander() +
        theme(legend.position = "none", plot.margin = margin(2, 2, 2, 2, "mm"))

ggarrange(tse_psychopathology_scatterplot,
          tse_cognition_scatterplot,
          tse_tbv_scatterplot,
          nrow = 1,
          ncol = 3) %>% 
    ggexport(
        filename = here("outputs", "figs", "tse_ar_scatterplot.pdf"),
        height = 4, width = 14
    )



# overall functioning and SRS (for presentation)
overallfunction_srs_scatterplot <- atypical_df %>% 
    ggplot(aes(x = tse_ar, y = -overall_functioning)) +
        geom_point(size = 3, color = "grey30", alpha = 0.9) +
        geom_smooth(method = "lm", color = "tomato3") +
        labs(x = "Traumatic Stress Load", y = "-Overall Functioning") +
        geom_point(aes(x = 1.2, y = 2.5), size = 3, color = "tomato3") +
        geom_label(
            label="Over-reacted", 
            x = 1.9,
            y = 2.5,
            label.padding = unit(0.55, "lines"), # Rectangle size around label
            label.size = 1,
            color = "grey10",
            fill = "tomato3"
        ) +
        geom_point(aes(x = 2.2, y = -1), size = 3, color = "royalblue") +
        geom_label(
            label="Resilient", 
            x = 2.7,
            y = -1,
            label.padding = unit(0.55, "lines"), # Rectangle size around label
            label.size = 1,
            color = "grey10",
            fill = "royalblue"
        ) +
        theme_pander()

overallfunction_srs_scatterplot
ggsave(plot = overallfunction_srs_scatterplot,
       filename = here("outputs", "figs", "overallfunction_srs_scatterplot.pdf"))
