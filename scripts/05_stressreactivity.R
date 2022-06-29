# Environment ------------------------------------------------------------------
library(tidyverse)
library(here)
library(ggpubr)
library(ggthemes)
library(glue)
library(ggridges)
library(ggseg)
library(ggsegGlasser)
library(ggsegYeo2011)


# Data I/O ---------------------------------------------------------------------
atypical_qri_df <- here("data", "processed", "preprocessed_data_atypical_withQRI.rds") %>% 
    read_rds()

glasser2yeo_df <- here("data", "raw", "glasser2yeo.csv") %>% 
    read_csv(show_col_types = FALSE)

# extract glasser labels
glasser_labels <- atypical_qri_df %>% 
    select(contains("glasser")) %>% 
    colnames()

subcortical_labels <- atypical_qri_df %>% 
    select(contains(c("Amygdala", "Caudate", "Thalamus.Proper", "Putamen", 
                      "Pallidum", "Hippocampus"))) %>% 
    colnames()

all_labels <- c(glasser_labels , subcortical_labels)

# Correlations: QRI against SRS ------------------------------------------------
atypical_qri_df %>% 
    select(aQRI, srs) %>% 
    correlation::correlation(method = "spearman")

corr_res <- atypical_qri_df %>% 
    select(aQRI, srs) %>% 
    psych::corr.test(method = "spearman")

# create correlation label for the hexagon plot
corr_res_r_label <- glue("italic(r)=={round(corr_res$r[1, 2], 3)}")
corr_res_p_label <- case_when(
    corr_res$p.adj < 0.001 ~ "~italic(p)<0.001",
    TRUE ~ as.character(glue("~italic(p)=={round(corr_res$p.adj, 3)}"))
)
corr_res_label <- paste(corr_res_r_label, corr_res_p_label, sep = "*','")

srs_qri_scatterplot <- atypical_qri_df %>% 
    ggplot(aes(x = srs, y = aQRI)) +
        geom_hex(binwidth = c(0.15, 0.01)) +
        #geom_point(size = 3, alpha = 0.9, color = "grey30") +
        geom_smooth(method = "lm", color = "tomato3", fill = "grey80") +
        #scale_fill_viridis_c(alpha = 0.95) +
        scale_fill_gradient(low = "grey75", high = "grey10") +
        annotate("text", x = 3, y = 0.19, label = corr_res_label, size = 4.5, parse = TRUE) +
        labs(x = "Stress Reactivity Score (SRS)",
             y = "Averaged Quantile Regression Index (aQRI)") +
        ggthemes::theme_pander() +
        theme(plot.margin = margin(2, 2, 2, 2, "mm"),
              legend.position = "none")
srs_qri_scatterplot
ggsave(filename = here("outputs", "figs", "srs_qri_scatterplot.pdf"),
       plot = srs_qri_scatterplot, height = 4, width = 5.5)

atypical_qri_df %>% 
    select(age, srs) %>% 
    correlation::correlation(method = "spearman")
age_srs_scatterplot <- atypical_qri_df %>% 
    ggplot(aes(x = age, y = srs)) +
    geom_point(size = 3, alpha = 0.9, color = "grey30") +
    geom_smooth(method = "lm", color = "tomato3", fill = "grey80") +
    #scale_fill_viridis_c(alpha = 0.95) +
    scale_fill_gradient(low = "grey75", high = "grey10") +
    #annotate("text", x = 3, y = 0.19, label = corr_res_label, size = 4.5, parse = TRUE) +
    labs(x = "Age",
         y = "Stress Reactivity Score") +
    ggthemes::theme_pander() +
    theme(plot.margin = margin(2, 2, 2, 2, "mm"),
          legend.position = "none")
age_srs_scatterplot
ggsave(filename = here("outputs", "figs", "age_srs_scatterplot.pdf"),
       plot = age_srs_scatterplot, height = 4, width = 5.5)


wilcox.test(srs ~ sex, atypical_qri_df)

# Distributions: QRI against 3 SRS Ranks ---------------------------------------
atypical_qri_df %>% 
    mutate(srs_rank = factor(ntile(srs, 3))) %>% 
    ggplot(aes(x = aQRI, y = srs_rank, fill = srs_rank)) +
        geom_density_ridges(scale = 2,
                            size = 0.45,
                            alpha = 0.95,
                            quantile_lines = TRUE,
                            quantile_fun = function(x, ...) mean(x)) +
        scale_x_continuous(limits = c(-0.08, 0.15)) +
        scale_discrete_manual("vline_color",
                              values = c("tomato3", "tomato3"), 
                              name = NULL) +
        scale_fill_manual(values = c("grey30", "grey50", "grey70", "grey90")) +
        labs(x = "Averaged Quantile Regression Index",
             y = "Stress Reactivity Score (Lower Rank = More Resilient)") +
        theme_pander() +
        scale_y_discrete(expand = expansion(add = c(0.3, 1.9))) +
        stat_compare_means(comparisons = list(c("1", "2"), c("2", "3"), c("1", "3"))) +
        theme(legend.position = "none",
              plot.margin = margin(2, 2, 2, 2, "mm"))


# plot distribution of each quantile of SRS
cortical_srs_long_df <- atypical_qri_df %>% 
    mutate(srs_rank = factor(ntile(srs, 3))) %>% 
    pivot_longer(
        cols = contains("glasser"),
        names_to = c("variable"),
        values_to = "region_QRI"
    ) %>% 
    mutate(variable = str_remove(variable, "mprage_glasser_")) %>% 
    separate(variable,
             into = c("modality", "label"),
             sep = "_",
             extra = "merge")

cortical_qri_df <- cortical_srs_long_df %>% 
    group_by(label, srs_rank) %>% 
    summarize(aQRI = mean(region_QRI, na.rm = TRUE)) %>% 
    ungroup() %>% 
    mutate(label = str_replace(label, "Right", "rh_R")) %>% 
    mutate(label = str_replace(label, "Left", "lh_L")) 

srs_rank_labels <- c("SRS Rank = 1", "SRS Rank = 2", "SRS Rank = 3")
srs_rank_levels <- c(1, 2, 3)

top_10_cortical_value <- cortical_qri_df %>% 
    arrange(desc(abs(aQRI))) %>% 
    slice(round(nrow(.) * 0.1)) %>% 
    pull(aQRI)

cortical_srs_brainplot <- cortical_qri_df %>% 
    mutate(aQRI = if_else(abs(aQRI) <= top_10_cortical_value, 0, aQRI)) %>% 
    mutate(srs_rank = factor(
            srs_rank, 
            levels = srs_rank_levels, 
            labels = srs_rank_labels
        )
    ) %>% 
    add_row(label = c("lh_???"), srs_rank = srs_rank_labels) %>% 
    add_row(label = c("rh_???"), srs_rank = srs_rank_labels) %>% 
    ggplot() +   
        geom_brain(atlas = glasser, 
                   color = "grey95", 
                   mapping = aes(fill = aQRI)) +
        scale_fill_gradientn(colors = c("royalblue", "grey90", "tomato3"),
                             na.value = "grey90",
                             limit = c(-0.1, 0.1)) +
        labs(fill = "aQRI") +
        theme_void() +
        facet_wrap(~srs_rank, ncol = 1, strip.position = "left") +
        theme(legend.position = "right") +
        theme(strip.text.y = element_text(size = 12),
              plot.margin = margin(2, 2, 2, 2, "mm"))
cortical_srs_brainplot
ggsave(filename = here("outputs", "figs", "srs_cortical_brainplot.pdf"),
       plot = cortical_srs_brainplot, width = 10, height = 5)

subcortical_qri_df <- atypical_qri_df %>% 
    mutate(srs_rank = factor(ntile(srs, 3))) %>% 
    select(-contains("glasser")) %>% 
    pivot_longer(cols = contains(subcortical_labels),
                 names_to = c("label"),
                 values_to = "region_QRI") %>% 
    mutate(label = gsub(".", "-", label, fixed = TRUE)) %>% 
    group_by(label, srs_rank) %>% 
    summarize(aQRI = mean(region_QRI, na.rm = TRUE), .groups = "drop") %>% 
    ungroup()

top_10_subcortical_value <- subcortical_qri_df %>% 
    arrange(desc(abs(aQRI))) %>% 
    slice(round(nrow(.) * 0.1)) %>% 
    pull(aQRI)

aseg_mapping_df <- subcortical_qri_df %>% 
    mutate(aQRI = if_else(abs(aQRI) < top_10_subcortical_value, 0, aQRI)) %>% 
    mutate(srs_rank = factor(srs_rank, levels = srs_rank_levels, labels = srs_rank_labels)) %>% 
    # add these rows to avoid NA group when use "facet_wrap"
    add_row(label = "Right-Lateral-Ventricle", srs_rank = srs_rank_labels) %>%
    add_row(label = "Left-Lateral-Ventricle", srs_rank = srs_rank_labels) %>%
    add_row(label = "Right-VentralDC", srs_rank = srs_rank_labels) %>%
    add_row(label = "Left-VentralDC", srs_rank = srs_rank_labels) %>% 
    add_row(label = "x3rd-ventricle", srs_rank = srs_rank_labels) %>% 
    add_row(label = "x4th-ventricle", srs_rank = srs_rank_labels) %>% 
    add_row(label = "brain-stem", srs_rank = srs_rank_labels) %>% 
    add_row(label = "cc-anterior", srs_rank = srs_rank_labels) %>% 
    add_row(label = "cc-central", srs_rank = srs_rank_labels) %>% 
    add_row(label = "cc-mid-anterior", srs_rank = srs_rank_labels) %>% 
    add_row(label = "cc-mid-posterior", srs_rank = srs_rank_labels) %>% 
    add_row(label = "cc-posterior", srs_rank = srs_rank_labels) %>% 
    add_row(label = "right-cerebellum-white-matter", srs_rank = srs_rank_labels) %>% 
    add_row(label = "right-cerebellum-cortex", srs_rank = srs_rank_labels) %>% 
    add_row(label = NA, srs_rank = srs_rank_labels) %>% 
    right_join(as_tibble(aseg))

subcortial_srs_brainplot <- aseg_mapping_df %>% 
    ggplot() +   
        geom_brain(atlas = aseg, 
                   color = "grey95", 
                   side = "axial",
                   mapping = aes(fill = aQRI)) +
        scale_fill_gradientn(colors = c("royalblue", "grey90", "tomato3"),
                             na.value = "grey90",
                             limits = c(-0.1, 0.1)) +
        labs(fill = "aQRI") +
        theme_void() +
        facet_wrap(~srs_rank, ncol = 1, strip.position = "left") +
        theme(legend.position = "right") +
        theme(strip.text.y = element_text(size = 12),
              plot.margin = margin(2, 2, 2, 2, "mm"))
subcortial_srs_brainplot
ggsave(filename = here("outputs", "figs", glue("srs_subcortical_brainplot.pdf")),
       plot = subcortial_srs_brainplot, width = 10, height = 5)



# plot brain of each TSE load over each modality
modality_qri_df <- cortical_srs_long_df %>% 
    group_by(label, modality, srs_rank) %>% 
    summarize(aQRI = mean(region_QRI, na.rm = TRUE)) %>% 
    ungroup() %>% 
    mutate(label = str_replace(label, "Right", "rh_R")) %>% 
    mutate(label = str_replace(label, "Left", "lh_L")) 

for (mod in c("gmd", "vol", "ct")) {
    figure <- modality_qri_df %>% 
        mutate(aQRI = if_else(abs(aQRI) < 0.037, 0, aQRI)) %>% 
        filter(modality == mod) %>% 
        mutate(srs_rank = factor(srs_rank, levels = srs_rank_levels, labels = srs_rank_labels)) %>% 
        add_row(label = c("lh_???"), srs_rank = srs_rank_labels) %>% 
        add_row(label = c("rh_???"), srs_rank = srs_rank_labels) %>% 
        ggplot() +   
            geom_brain(atlas = glasser, 
                       color = "white", 
                       mapping = aes(fill = aQRI)) +
            scale_fill_gradientn(colors = c("royalblue", "grey95", "tomato3"),
                                 na.value = "grey95",
                                 limits = c(-0.1, 0.1)) +
            theme_void() +
            facet_wrap(~srs_rank, ncol = 1, strip.position = "left") +
            theme(legend.position = "right") +
            theme(strip.text.y = element_text(size = 12),
                  plot.margin = margin(2, 2, 2, 2, "mm"))
    
    ggsave(filename = here("outputs", "figs", glue("srs_cortical_{mod}_brain.pdf")),
           plot = figure, width = 10, height = 5)
}


# compare upper and lower tertile ----------------------------------------------
all_res <- NULL
all_res$label <- all_labels
all_res$statistic <- NULL
all_res$pvalue <- NULL

for (label in all_labels) {
    f <- as.formula(glue("`{label}` ~ srs_rank"))
    
    res <- atypical_qri_df  %>% 
        select(srs, all_of(label)) %>% 
        correlation::correlation(method = "spearman")
   
   
    
    all_res$statistic <- c(all_res$statistic, res$r)
    all_res$pvalue <- c(all_res$pvalue, res$p)
}

all_res$label[(all_res$pvalue %>% p.adjust(method = "fdr")) < 0.05]


# Comparison between High and Low SRS ------------------------------------------

network_qri_df <- as_tibble(
    list(
        network = rep(1:7),
        p = rep(NA, 7),
        effectsize = rep(NA, 7)
    )
)

for (network in 1:7) {
    
    idx <- glasser2yeo_df %>% 
        filter(yeo == network) %>% 
        select(glasser) %>% 
        pull()
    
    
    tmp_df <- atypical_qri_df %>% 
        mutate(srs_rank = factor(ntile(srs, 3))) %>% 
        select(contains("glasser"), bblid, srs_rank) %>% 
        select(idx, idx + 360, idx + 720, bblid, srs_rank) %>%
        pivot_longer(cols = contains("glasser"),
                     names_to = c("label"),
                     values_to = "QRI") %>% 
        filter(srs_rank %in% c(1, 3)) %>% 
        mutate(srs_rank = droplevels(srs_rank)) %>%  
        group_by(bblid, srs_rank) %>% 
        summarize(
            network_qri = mean(QRI),
            .groups = "drop"
        )
    
    tmp_df %>% 
        group_by(srs_rank) %>% 
        summarize(
            network_qri = mean(network_qri),
            .groups = "drop"
        ) %>% 
        print()
    
    
    network_qri_df[network, 2] <- 
        wilcox.test(formula = network_qri ~ srs_rank, data = tmp_df) %>% 
        broom::tidy() %>% 
        pull(p.value)
    
    network_qri_df[network, 3] <- 
        rstatix::wilcox_effsize(formula = network_qri ~ srs_rank, data = tmp_df) %>% 
        pull(effsize)
}

comparison_df <- network_qri_df %>% 
    mutate(adj.p = p.adjust(network_qri_df$p, method = "fdr"), .after = "p") %>% 
    mutate(effectsize = if_else(adj.p < 0.05, effectsize, 0)) %>% 
    filter(adj.p < 0.05)

# visualize Yeo 7 Networks
comparison_brainplot <- as_tibble(
        list(
            label = paste0("lh_7Networks_", comparison_df$network),
            es = comparison_df$effectsize
        )
    ) %>% 
    ggseg(
        atlas = "yeo7", 
        mapping = aes(fill = es), 
        hemisphere = "left",
        color = "white"
    ) +
    scale_fill_viridis_c(begin = 0.2, na.value = "grey85") +
    labs(x = "", fill = "Effect Size", title = "Yeo 7 Networks:  High (SRS Rank = 3) - Low (SRS Rank = 1)") +
    theme_brain(text.size = 12, text.family = "sans") 
comparison_brainplot
ggsave(
    filename = here("outputs", "figs", "comparison_high_low_srs_brainplot.pdf"),
    plot = comparison_brainplot,
    width = 8,
    height = 3
)


tmp_df <- atypical_qri_df %>% 
    mutate(srs_rank = factor(ntile(srs, 3))) %>% 
    select(contains(c("Left.Amygdala")), bblid, srs_rank) %>% 
    pivot_longer(cols = !c(bblid, srs_rank),
                 names_to = c("label"),
                 values_to = "QRI") %>% 
    filter(srs_rank %in% c(1, 3)) %>% 
    mutate(srs_rank = droplevels(srs_rank)) %>% 
    group_by(bblid, srs_rank) %>% 
    summarize(
        network_qri = mean(QRI),
        .groups = "drop"
    )

tmp_df %>% 
     group_by(srs_rank) %>% 
     summarize(
         network_qri_mean = mean(network_qri),
         network_qri_sd = sd(network_qri),
         .groups = "drop"
     ) %>% 
     print()
 
wilcox.test(formula = network_qri ~ srs_rank, data = tmp_df) %>% 
     broom::tidy() 

rstatix::wilcox_effsize(formula = network_qri ~ srs_rank, data = tmp_df)
 
 
atypical_qri_df %>% 
     select(contains(c("Left.Amygdala", "Putamen")), bblid, srs) %>% 
     pivot_longer(cols = !c(bblid, srs),
                  names_to = c("label"),
                  values_to = "QRI") %>% 
     group_by(bblid, srs) %>% 
     summarize(
         network_qri = mean(QRI),
         .groups = "drop"
     ) %>% 
     select(srs, network_qri) %>% 
     correlation::correlation(method = "spearman")
 
 
 
# Figure 4C --------------------------------------------------------------------

combined_qri_df <- cortical_qri_df %>% 
    select(label, srs_rank, aQRI) %>% 
    mutate(modality = "cortical", .after = "label") %>% 
    bind_rows(modality_qri_df) %>% 
    bind_rows(subcortical_qri_df %>% mutate(modality = "subcortical")) %>% 
    filter(srs_rank %in% c(1, 3)) %>% 
    mutate(srs_rank = factor(srs_rank, levels = c(1, 3), labels = c("SRS Rank = 1", "SRS Rank = 3"))) %>% 
    mutate(modality = factor(
            modality, 
            levels = c("cortical", "ct", "vol", "gmd", "subcortical"),
            labels = c("Cortical", "Thickness", "Volume", "Density", "Subcortical")
        )
    )
 
aqri_comparison_boxplot <- combined_qri_df %>% 
    ggplot(aes(x = srs_rank, y = aQRI, fill = srs_rank)) +
        gg.layers::geom_boxplot2(width = 0.6, width.errorbar = 0.3) +
        facet_wrap(~modality, strip.position = "bottom", nrow = 1) +
        scale_fill_manual(values = c("grey50", "grey90")) +
        labs(x = "", y = "Averaged Quantile Regression Index (aQRI)", fill = "") +
        theme_pander() +
        theme(axis.text.x = element_blank(),
              plot.margin = margin(2, 2, 2, 2, "mm"),
              legend.position = "top",
              legend.key.size = unit(1.2, "cm")) +
        stat_compare_means(aes(label = ..p.signif..), label.y = 0.1, label.x = 1.5)

aqri_comparison_boxplot
ggsave(filename = here("outputs", "figs", "aqri_comparison_boxplot.pdf"),
       plot = aqri_comparison_boxplot, width = 14, height = 4)
     
