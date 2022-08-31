# Environment ------------------------------------------------------------------
library(tidyverse)
library(here)
library(ggpubr)
library(ggthemes)
library(gt)
library(gtsummary)
library(patchwork)
library(flextable)
library(officer)

# docx page setup
sect_properties <- prop_section(
    page_size = page_size(),
    type = "continuous",
    page_margins = page_mar(
        bottom = 0.5, top = 0.5, right = 0.5, left = 0.5, gutter = 0
    )
)


# Data I/O ---------------------------------------------------------------------
preprocessed_df <- here("data", "processed", "preprocessed_data.rds") %>% 
    read_rds()

normativel_df <- filter(preprocessed_df, is_td == 1 & tse == 0)
atypical_df <- filter(preprocessed_df, !c(is_td == 1 & tse == 0))


# TSE Distribution -------------------------------------------------------------
tse_distplot <- atypical_df %>% 
    ggplot(aes(x = tse)) +
     geom_histogram(binwidth = 1) +
        labs(x = "Traumatic Stress Load")

tsl_distplot <- atypical_df %>% 
    ggplot(aes(x = tse_ar)) +
        geom_histogram(binwidth = 0.5) +
        labs(x = "Age-and-sex-adjusted Traumatic Stress Load (a.u.)")

tse_distribution_plot <- tse_distplot / 
    tsl_distplot &
    theme_pander() +
    theme(plot.margin = margin(5, 5, 5, 5, "mm"))
tse_distribution_plot

tse_distribution_plot %>% 
    ggexport(
        filename = here("outputs", "figs", "tse_distributions.pdf"),
        width = 6, 
        height = 6
    )


# Table ------------------------------------------------------------------------
tse_comparison_table <- atypical_df %>% 
    mutate(tse_ar_rank = ntile(tse, 3)) %>% 
    bind_rows(normativel_df %>% mutate(tse_ar_rank = 0)) %>% 
    mutate(tse_ar_rank = factor(tse_ar_rank, levels = c(0, 1, 2, 3),
                                labels = c("Normative", "Low", "Moderate", "High"))) %>% 
    mutate(mprage_antsCT_vol_TBV = as.numeric(scale(mprage_antsCT_vol_TBV))) %>% 
    select(tse_ar_rank, age:Overall_Efficiency_Ar, 
           F4_Executive_Efficiency_Ar, F3_Memory_Efficiency_Ar,
           F2_Complex_Reasoning_Efficiency_Ar, F1_Social_Cognition_Efficiency_Ar,
           mprage_antsCT_vol_TBV) %>% 
    select(-c(is_td, tse, trauma_type)) %>% 
    tbl_summary(
        by = "tse_ar_rank",
        statistic = list(all_continuous() ~ "{mean} ({sd})"),
        label = list(
            age ~ "Age",
            sex ~ "Sex",
            race2 ~ "Race",
            medu1 ~ "Maternal Education",
            envSES ~ "Environmental SES",
            overall_functioning ~ "Functioning Composite Score",
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
    fontsize(size = 8.5, part = "all") %>% 
    padding(padding.top = 1, padding.bottom = 1, part = "all") %>% 
    bold(part = "header") %>% 
    set_table_properties(width = 1, layout = "autofit") %>% 
    save_as_docx(
        path = here("outputs", "tables", "table1.docx"),
        pr_section = sect_properties
    )

tse_comparison_table %>% 
    as_gt() %>% 
    gtsave(
        filename = "table1.html", 
        path = here("outputs", "tables")
    )


# Visualization ----------------------------------------------------------------
# control for sex and age effect
corr_tse_outcomes <- atypical_df %>% 
    select(tse_ar, overall_psychopathology_4factorv2,
           Overall_Efficiency_Ar, mprage_antsCT_vol_TBV) %>% 
    correlation::correlation(method = "spearman", p_adjust = "fdr")
corr_tse_outcomes

tse_psychopathology_scatterplot <- atypical_df %>% 
    ggplot(aes(x = tse_ar, y = overall_psychopathology_4factorv2)) +
        #geom_point(color = "grey30", alpha = 0.9, size = 3) +
        geom_hex(binwidth = c(0.25, 0.25)) +
        scale_fill_gradient(low = "grey70", high = "grey10") +
        geom_smooth(method = "lm", color = "tomato3", fill = "grey80") +
        annotate(
            geom = "text", 
            x = 4, 
            y = -3, 
            label = "italic(r)==0.25*','~italic(p)<0.001", 
            size = 4.5, 
            parse = TRUE
        ) +
        ylim(-3, 3) +
        labs(x = "Age-and-sex-adjusted Traumatic Stress Load",
             y = "Psychopathology (g)")
tse_cognition_scatterplot <- atypical_df %>% 
    ggplot(aes(x = tse_ar, y = Overall_Efficiency_Ar)) +
        #geom_point(color = "grey30", alpha = 0.9, size = 3) +
        geom_hex(binwidth = c(0.25, 0.25)) +
        scale_fill_gradient(low = "grey70", high = "grey10") +
        geom_smooth(method = "lm", color = "tomato3", fill = "grey80") +
        annotate(
            geom = "text", 
            x = 4, 
            y = -3, 
            label = "italic(r)==-0.09*','~italic(p)==0.005", 
            size = 4.5, 
            parse = TRUE
        ) +
        ylim(-3, 3) +
        labs(x = "Age-and-sex-adjusted Traumatic Stress Load",
             y = "Cognitive Efficiency (g)")
tse_tbv_scatterplot <- atypical_df %>% 
    mutate(mprage_antsCT_vol_TBV = as.numeric(scale(mprage_antsCT_vol_TBV))) %>% 
    ggplot(aes(x = tse_ar, y = mprage_antsCT_vol_TBV)) +
        #geom_point(color = "grey30", alpha = 0.8, size = 3) +
        geom_hex(binwidth = c(0.25, 0.25)) +
        scale_fill_gradient(low = "grey70", high = "grey10") +
        geom_smooth(method = "lm", color = "tomato3", fill = "grey80") +
        annotate(
            geom = "text", 
            x = 4, 
            y = -3, 
            label = "italic(r)==-0.09*','~italic(p)==0.005", 
            size = 4.5, 
            parse = TRUE
        ) +
        ylim(-3, 3) +
        labs(x = "Age-and-sex-adjusted Traumatic Stress Load",
             y = "Total Brain Volume")

tse_ar_scatterplot <- tse_psychopathology_scatterplot +
    tse_cognition_scatterplot +
    tse_tbv_scatterplot &
    theme_pander() +
    theme(legend.position = "none", plot.margin = margin(2, 2, 2, 2, "mm")) 

tse_ar_scatterplot
tse_ar_scatterplot %>% 
    ggexport(
        filename = here("outputs", "figs", "tse_ar_scatterplot.pdf"),
        height = 4, width = 14
    )


# Not for the Main Analysis ----------------------------------------------------
# overall functioning and SRS (for presentation)
overallfunction_srs_scatterplot <- atypical_df %>% 
    ggplot(aes(x = tse_ar, y = -overall_functioning)) +
        geom_point(size = 3, color = "grey30", alpha = 0.9) +
        geom_smooth(method = "lm", color = "tomato3") +
        labs(x = "Traumatic Stress Load", y = "-Overall Functioning") +
        geom_point(aes(x = 1.2, y = 2.5), size = 3, color = "tomato3") +
        geom_label(
            aes(size = 30),
            label="Over-reacted", 
            x = 2.45,
            y = 2.5,
            label.padding = unit(0.55, "lines"), # Rectangle size around label
            label.size = 0.8,
            color = "grey10",
            fill = "tomato3"
        ) +
        geom_point(aes(x = 2.2, y = -1), size = 3, color = "royalblue") +
        geom_label(
            aes(size = 30),
            label="Resilient", 
            x = 3.3,
            y = -1,
            label.padding = unit(0.55, "lines"), # Rectangle size around label
            label.size =  0.8,
            color = "grey10",
            fill = "royalblue"
        ) +
        theme_pander() +
        theme(
            legend.position = "none",
            plot.margin = margin(5, 5, 5, 5, "mm")
        )
overallfunction_srs_scatterplot
ggsave(plot = overallfunction_srs_scatterplot,
       filename = here("outputs", "figs", "overallfunction_srs_scatterplot.pdf"),
       width = 6,
       height = 4
)
