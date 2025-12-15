# Environment ------------------------------------------------------------------
library(tidyverse)
library(here)
library(ggthemes)
library(ggtext)
library(patchwork)
library(lme4)
library(flextable)
library(officer)

source(here("src", "R", "plot_ggseg_brain.R"))
source(here("src", "R", "extract_mediation_summary.R"))

sect_properties <- prop_section(
    page_size = page_size(
        orient = "landscape",
        width = 8.3, height = 11.7
    ),
    type = "continuous",
    page_margins = page_mar(gutter = 0)
)

.save_effect2docx <- function(df, pathname) {
    df %>%
        mutate(
            p.value = case_when(
                p.value == 0 ~ "< 2e-16",
                p.value >= 0.001 ~ as.character(round(p.value, 3)),
                p.value < 0.001 ~ paste(formatC(p.value, format = "e", digits = 2))
            ),
            p.adj = case_when(
                p.adj == 0 ~ "< 2e-16",
                p.adj >= 0.001 ~ as.character(round(p.adj, 3)),
                p.adj < 0.001 ~ paste(formatC(p.adj, format = "e", digits = 2))
            ),
            estimate = round(estimate, 3),
            std.error = round(std.error, 3),
            statistic = round(statistic, 3)
        ) %>%
        arrange(desc(abs(estimate))) %>% 
        rename_all(
            ~c("Dataset", "Variable Name", "ggseg Label", "Term", "Std. β", 
               "SE", "T", "p", "p(FDR)")
        ) %>% 
        flextable() %>% 
        bold(part = "header") %>% 
        italic(j = c(3, 7:9), part = "header") %>%
        autofit() %>% 
        save_as_docx(
            path = here("outputs", "tables", glue::glue("{pathname}.docx")),
            pr_section = sect_properties
        )
    return(cat("your file is saved"))
}

write_rds(master_srs_brain_df, here("outputs", "caches", "master_srs_brain_df.rds"))

# Data IO ----------------------------------------------------------------------
gmv_variables_df <- here("data", "raw", "included_variables.xlsx") %>%
    readxl::read_excel(sheet = "GMV") %>%
    select(variable_name, ggseg_label)

# ABCD data folder path
# using data release 5.1
abcd_folder_path <- "path_to_abcd_folder" %>%
    here("abcd-data-release-5.1/core")

# general
demo_y_age <- here(abcd_folder_path, "abcd-general", "abcd_y_lt.csv") %>%
    read_csv(show_col_types = FALSE) %>%
    select(src_subject_id, eventname, interview_age)
mri_y_adm_info <- here(abcd_folder_path, "imaging", "mri_y_adm_info.csv") %>%
    read_csv(show_col_types = FALSE)
mri_y_qc_clfind <- here(abcd_folder_path, "imaging", "mri_y_qc_clfind.csv") %>%
    read_csv(show_col_types = FALSE)
mri_y_qc_motion <- here(abcd_folder_path, "imaging", "mri_y_qc_motion.csv") %>%
    read_csv(show_col_types = FALSE)
mri_y_qc_incl <- here(abcd_folder_path, "imaging", "mri_y_qc_incl.csv") %>%
    read_csv(show_col_types = FALSE)
mri_y_vol_dsk <- here(
    abcd_folder_path, "imaging", "mri_y_smr_vol_dsk.csv"
) %>%
    read_csv(show_col_types = FALSE)
mri_y_vol_aseg <- here(
    abcd_folder_path, "imaging", "mri_y_smr_vol_aseg.csv"
) %>%
    read_csv(show_col_types = FALSE)


master_srs_brain_df <- read_rds(
    file = here("outputs", "caches", "master_srs_df.rds")
) %>%
    select(
        src_subject_id, dataset, site_id_l, avg_srs, race_ethnicity_3, 
        demo_sex_v2, `social services`:`family SES`
    ) %>%
    group_by(src_subject_id) %>%
    slice(1) %>%
    ungroup() %>%
    right_join(
        mri_y_qc_motion,
        by = c("src_subject_id"),
        relationship = "many-to-many"
    ) %>%
    left_join(demo_y_age, by = c("src_subject_id", "eventname")) %>%
    right_join(mri_y_qc_clfind, by = c("src_subject_id", "eventname")) %>%
    right_join(mri_y_qc_incl, by = c("src_subject_id", "eventname")) %>%
    right_join(mri_y_adm_info, by = c("src_subject_id", "eventname")) %>%
    right_join(mri_y_vol_dsk, by = c("src_subject_id", "eventname")) %>%
    right_join(mri_y_vol_aseg, by = c("src_subject_id", "eventname")) %>%
    filter(
        # https://wiki.abcdstudy.org/release-notes/imaging/quality-control.html
        mrif_score %in% c(1, 2) &
            imgincl_t1w_include == 1
    ) %>%
    mutate(across(contains("smri_"), ~ as.numeric(scale(.)))) %>%
    drop_na(avg_srs)


# Linear Mixed Model -----------------------------------------------------------
cgmv_lmm_d_fit <- master_srs_brain_df %>%
    filter(dataset == "Discovery") %>% 
    lmerTest::lmer(
        formula = smri_vol_cdk_total ~ scale(interview_age) * avg_srs +
            demo_sex_v2 + race_ethnicity_3 + (1 | site_id_l) +
            (1 | src_subject_id)
    )
broom.mixed::tidy(cgmv_lmm_d_fit)

cgmv_lmm_d_res <- broom.mixed::tidy(cgmv_lmm_d_fit) %>%
    mutate(estimate = round(estimate, 3)) %>%
    mutate(p.value = case_when(
        p.value == 0 ~ "< 2e-16",
        p.value >= 0.001 ~ as.character(round(p.value, 3)),
        p.value < 0.001 ~ paste(formatC(p.value, format = "e", digits = 2)))
    ) %>% 
    slice(2, 3, 7) %>% 
    select(-c(1:2)) %>% 
    mutate(dataset = "Discovery", .before = term) %>% 
    mutate(modality = "Cortical GMV", .before = term)
cgmv_lmm_h_fit <- master_srs_brain_df %>%
    filter(dataset == "Holdout") %>% 
    lmerTest::lmer(
        formula = smri_vol_cdk_total ~ scale(interview_age) * avg_srs +
            demo_sex_v2 + race_ethnicity_3 + (1 | site_id_l) +
            (1 | src_subject_id)
    )
cgmv_lmm_h_res <- broom.mixed::tidy(cgmv_lmm_h_fit) %>%
    mutate(estimate = round(estimate, 3)) %>%
    mutate(p.value = case_when(
        p.value == 0 ~ "< 2e-16",
        p.value >= 0.001 ~ as.character(round(p.value, 3)),
        p.value < 0.001 ~ paste(formatC(p.value, format = "e", digits = 2)))
    ) %>% 
    slice(2, 3, 7) %>% 
    select(-c(1:2)) %>%
    mutate(dataset = "Holdout", .before = term) %>% 
    mutate(modality = "Cortical GMV", .before = term)

annotation_text <- glue::glue(
    paste0(
        "Effect of Age: *β* = {cgmv_lmm_d_res$estimate[1]}, ",
        "*p* {cgmv_lmm_d_res$p.value[1]}<br>",
        "Effect of SRS: *β* = {cgmv_lmm_d_res$estimate[2]}, ",
        "*p* = {cgmv_lmm_d_res$p.value[2]}<br>",
        "Interaction Effect: *β* = {cgmv_lmm_d_res$estimate[3]}, ",
        "*p* = {cgmv_lmm_d_res$p.value[3]}<br>",
        "Replicated in Holdout: Yes"
    )
)

master_srs_brain_df %>%
    filter(dataset == "Discovery") %>% 
    pull(avg_srs) %>% 
    quantile(probs = c(0.05, 0.5, 0.95))

cgmv_d_int_fig <- cgmv_lmm_d_fit %>%
    sjPlot::plot_model(
        type = "pred",
        terms = c("interview_age", "avg_srs [-1.44, -0.014, 1.52]")
    ) +
    scale_y_continuous(limits = c(-0.1, 1.25), breaks = seq(0, 1.2, 0.4)) +
    scale_color_manual(
        values = c("#E69F00", "gray10", "#E41A1C"),
        labels = c(
            "High Resilience (5%)",
            "Normative Response (50%)",
            "High Susceptibility (95%)"
        )
    ) +
    scale_fill_manual(
        values = c("#E69F00", "gray10", "#E41A1C"),
        labels = c(
            "High Resilience (5%)",
            "Normative Response (50%)",
            "High Susceptibility (95%)"
        )
    ) +
    labs(
        x = "Interview Age (Months)",
        y = "Total Cortical Gray Matter Volume (Z-score)",
        color = "",
        fill = "",
        title = ""
    ) +
    theme_pander() +
    theme(
        legend.position = "top",
        plot.margin = margin(5, 5, 5, 5, "mm")
    ) +
    geom_richtext(
        data = tibble(x = 100, y = 0.27),
        mapping = aes(x = x, y = y),
        label = annotation_text,
        hjust = 0,
        vjust = 1,
        size = 4,
        fill = NA,
        label.color = NA,
        color = "gray30",
        family = "Helvetica Neue",
        inherit.aes = FALSE
    )
cgmv_d_int_fig


sgmv_lmm_d_fit <- master_srs_brain_df %>%
    filter(dataset == "Discovery") %>% 
    lmerTest::lmer(
        formula = smri_vol_scs_subcorticalgv ~ scale(interview_age) * avg_srs +
            demo_sex_v2 + race_ethnicity_3 + (1 | site_id_l) +
            (1 | src_subject_id)
    )
sgmv_lmm_d_res <- broom.mixed::tidy(sgmv_lmm_d_fit) %>%
    mutate(estimate = round(estimate, 3)) %>%
    mutate(p.value = case_when(
        p.value == 0 ~ "< 2e-16",
        p.value >= 0.001 ~ as.character(round(p.value, 3)),
        p.value < 0.001 ~ paste(formatC(p.value, format = "e", digits = 2)))
    ) %>% 
    slice(2, 3, 7) %>% 
    select(-c(1:2)) %>% 
    mutate(dataset = "Discovery", .before = term) %>% 
    mutate(modality = "Subcortical GMV", .before = term)

sgmv_lmm_h_fit <- master_srs_brain_df %>%
    filter(dataset == "Holdout") %>% 
    lmerTest::lmer(
        formula = smri_vol_scs_subcorticalgv ~ scale(interview_age) * avg_srs +
            demo_sex_v2 + race_ethnicity_3 + (1 | site_id_l) +
            (1 | src_subject_id)
    )
sgmv_lmm_h_res <- broom.mixed::tidy(sgmv_lmm_h_fit) %>%
    mutate(estimate = round(estimate, 3)) %>%
    mutate(p.value = case_when(
        p.value == 0 ~ "< 2e-16",
        p.value >= 0.001 ~ as.character(round(p.value, 3)),
        p.value < 0.001 ~ paste(formatC(p.value, format = "e", digits = 2)))
    ) %>% 
    slice(2, 3, 7) %>% 
    select(-c(1:2)) %>% 
    mutate(dataset = "Holdout", .before = term) %>% 
    mutate(modality = "Subcortical GMV", .before = term)

annotation_text <- glue::glue(
    paste0(
        "Effect of Age: *β* = {sgmv_lmm_d_res$estimate[1]}, ",
        "*p* {sgmv_lmm_d_res$p.value[1]}<br>",
        "Effect of SRS: *β* = {sgmv_lmm_d_res$estimate[2]}, ",
        "*p* = {sgmv_lmm_d_res$p.value[2]}<br>",
        "Interaction Effect: *β* = {sgmv_lmm_d_res$estimate[3]}, ",
        "*p* = {sgmv_lmm_d_res$p.value[3]}<br>",
        "Replicated in Holdout: Yes"
    )
)

sgmv_d_int_fig <- sgmv_lmm_d_fit %>%
    sjPlot::plot_model(
        type = "pred",
        terms = c("interview_age", "avg_srs [-1.44, -0.014, 1.52]")
    ) +
    scale_y_continuous(limits = c(-0.1, 1.2), breaks = seq(0, 1.2, 0.4)) +
    scale_color_manual(
        values = c("#E69F00", "gray10", "#E41A1C"),
        labels = c(
            "High Resilience (5%)",
            "Normative Response (50%)",
            "High Susceptibility (95%)"
        )
    ) +
    scale_fill_manual(
        values = c("#E69F00", "gray10", "#E41A1C"),
        labels = c(
            "High Resilience (5%)",
            "Normative Response (50%)",
            "High Susceptibility (95%)"
        )
    ) +
    labs(
        x = "Interview Age (Months)",
        y = "Total Subcortical Gray Matter Volume (Z-score)",
        color = "",
        fill = "",
        title = ""
    ) +
    theme_pander() +
    theme(
        legend.position = "top",
        plot.margin = margin(5, 5, 5, 5, "mm")
    ) +
    geom_richtext(
        data = tibble(x = 100, y = 0.27),
        mapping = aes(x = x, y = y),
        label = annotation_text,
        hjust = 0,
        vjust = 1,
        size = 4,
        fill = NA,
        label.color = NA,
        color = "gray30",
        family = "Helvetica Neue",
        inherit.aes = FALSE
    )
sgmv_d_int_fig

interaction_plot <- cgmv_d_int_fig +
    sgmv_d_int_fig +
    plot_layout(guides = "collect") &
    theme(legend.position = "top")
interaction_plot
ggsave(
    plot = interaction_plot,
    filename = here("outputs", "figures", "gmv_discovery_interaction.pdf"),
    height = 5,
    width = 10,
    device = cairo_pdf
)

bind_rows(
    cgmv_lmm_d_res,
    sgmv_lmm_d_res,
    cgmv_lmm_h_res,
    sgmv_lmm_h_res
) %>% 
    mutate(term = recode(
        term, 
        "scale(interview_age)" = "Age",
        "avg_srs" = "aSRS",
        "scale(interview_age):avg_srs" = "aSRS by Age"
    )) %>% 
    select(-df) %>% 
    mutate(
        std.error = round(std.error, 3),
        statistic = round(statistic, 2),
    ) %>% 
    rename_all(
        ~c("Dataset", "Modality", "Term", "Std. β", "SE", "t-value", "p")
    ) %>% 
    flextable() %>% 
    bold(part = "header") %>% 
    italic(j = c(3, 7), part = "header") %>%
    autofit() %>% 
    save_as_docx(
        path = here("outputs", "tables", "table_globalGMV_development.docx"),
        pr_section = sect_properties
    )


# Regional ----------------------------------------------
gmv_region_labels <- gmv_variables_df$variable_name[c(1:68, 70:83)]

region_d_results <- list()
for (ith in seq_along(gmv_region_labels)) {
    region_label <- gmv_region_labels[ith]
    formula_text <- glue::glue(
        "{region_label} ~ scale(interview_age) * avg_srs +
            demo_sex_v2 + race_ethnicity_3 + smri_vol_scs_intracranialv +
            (1 | site_id_l) + (1 | src_subject_id)"
    )
    region_d_results[[ith]] <- master_srs_brain_df %>%
        filter(dataset == "Discovery") %>% 
        lmerTest::lmer(
            formula = as.formula(formula_text),
            control = lmerControl(optimizer = "nloptwrap", calc.derivs = FALSE)
        ) %>%
        broom.mixed::tidy() %>%
        slice(2, 3, 8) %>%
        mutate(dataset = "Discovery", .before = term) %>% 
        mutate(region = region_label, .before = term) 
}
gmv_region_d_res <- region_d_results %>%
    reduce(bind_rows) %>%
    rename(
        variable_name = region
    ) %>% 
    select(-c(1, 2, 9))

srs_cgmv_d_effects <- gmv_region_d_res %>%
    filter(term == "avg_srs") %>%
    slice(1:68) %>% 
    mutate(p.adj = p.adjust(p.value, method = "fdr")) %>%
    left_join(gmv_variables_df, by = "variable_name") %>%
    relocate(ggseg_label, .before = term) %>% 
    rename(label = ggseg_label)
.save_effect2docx(
    srs_cgmv_d_effects, 
    pathname = "table_GMV_regional_cortical_d_development"
)

srs_sgmv_d_effects <- gmv_region_d_res %>%
    filter(term == "avg_srs") %>%
    slice(69:82) %>% 
    mutate(p.adj = p.adjust(p.value, method = "fdr")) %>%
    left_join(gmv_variables_df, by = "variable_name") %>%
    relocate(ggseg_label, .before = term) %>% 
    rename(label = ggseg_label)
.save_effect2docx(
    srs_sgmv_d_effects, 
    pathname = "table_GMV_regional_subcortical_d_development"
)

age_gmv_d_effects <- gmv_region_d_res %>%
    filter(term == "scale(interview_age)") %>%
    mutate(p.adj = p.adjust(p.value, method = "fdr")) %>%
    left_join(gmv_variables_df, by = "variable_name") %>%
    rename(label = ggseg_label)
   
int_gmv_d_effects <- gmv_region_d_res %>%
    filter(term == "scale(interview_age):avg_srs") %>%
    mutate(p.adj = p.adjust(p.value, method = "fdr")) %>%
    left_join(gmv_variables_df, by = "variable_name") %>%
    rename(label = ggseg_label)


region_h_results <- list()
for (ith in seq_along(gmv_region_labels)) {
    region_label <- gmv_region_labels[ith]
    formula_text <- glue::glue(
        "{region_label} ~ scale(interview_age) * avg_srs +
            demo_sex_v2 + race_ethnicity_3 + smri_vol_scs_intracranialv +
            (1 | site_id_l) + (1 | src_subject_id)"
    )
    region_h_results[[ith]] <- master_srs_brain_df %>%
        filter(dataset == "Holdout") %>% 
        lmerTest::lmer(
            formula = as.formula(formula_text),
            control = lmerControl(optimizer = "nloptwrap", calc.derivs = FALSE)
        ) %>%
        broom.mixed::tidy() %>%
        slice(2, 3, 8) %>%
        mutate(dataset = "Discovery", .before = term) %>% 
        mutate(region = region_label, .before = term) 
}

gmv_region_h_res <- region_h_results %>%
    reduce(bind_rows) %>%
    rename(variable_name = region) %>% 
    select(-c(1, 2, 9))
srs_cgmv_h_effects <- gmv_region_h_res %>%
    filter(term == "avg_srs") %>%
    slice(1:68) %>% 
    mutate(p.adj = p.adjust(p.value, method = "fdr")) %>%
    left_join(gmv_variables_df, by = "variable_name") %>%
    relocate(ggseg_label, .before = term) %>% 
    rename(label = ggseg_label)
.save_effect2docx(
    srs_cgmv_h_effects,
    pathname = "table_GMV_regional_cortical_h_development"
)
srs_sgmv_h_effects <- gmv_region_h_res %>%
    filter(term == "avg_srs") %>%
    slice(69:82) %>% 
    mutate(p.adj = p.adjust(p.value, method = "fdr")) %>%
    left_join(gmv_variables_df, by = "variable_name") %>%
    relocate(ggseg_label, .before = term) %>% 
    rename(label = ggseg_label)
.save_effect2docx(
    srs_sgmv_h_effects,
    pathname = "table_GMV_regional_subcortical_h_development"
)

appear_both <- intersect(
    srs_cgmv_d_effects %>% 
        filter(p.adj < 0.05) %>% 
        pull(label),
    srs_cgmv_h_effects %>% 
        filter(p.adj < 0.05) %>% 
        pull(label)
)


fig_a <- plot_ggseg_brain(
    dat = age_gmv_d_effects %>% filter(p.adj < 0.05),
    atlas = "dk",
    fill = "statistic",
    min = -120,
    max = 120,
    break_int = 30,
    title = "Effect of Age"
)
fig_b <- plot_ggseg_brain(
    dat = srs_cgmv_d_effects %>% filter(p.adj < 0.05),
    atlas = "dk",
    fill = "statistic",
    min = -6, 
    max = 6, 
    break_int = 3,
    title = "Effect of aSRS",
    highlight_regions = appear_both
)
fig_c <- plot_ggseg_brain(
    dat = int_gmv_d_effects %>% filter(p.adj < 0.05),
    atlas = "dk",
    fill = "statistic",
    min = -6,
    max = 6,
    break_int = 3,
    title = "Interaction Effect (Age by aSRS)"
)
fig4 <- fig_a / fig_b / fig_c
fig4
ggsave(
    plot = fig4,
    filename = here("outputs", "figures", "cgmv_d_brain.svg"),
    height = 5,
    width = 10,
    dpi = 600
)

fig_a <- plot_ggseg_brain(
    dat = age_gmv_d_effects %>% filter(p.adj < 0.05),
    atlas = "aseg",
    fill = "statistic",
    min = -150,
    max = 150,
    break_int = 60,
    title = " "
)
fig_b <- plot_ggseg_brain(
    dat = srs_sgmv_d_effects %>% filter(p.adj < 0.05),
    atlas = "aseg",
    fill = "statistic",
    min = -6,
    max = 6,
    break_int = 3,
    title = " "
)
fig_c <- plot_ggseg_brain(
    dat = int_gmv_d_effects %>% filter(p.adj < 0.05),
    atlas = "aseg",
    fill = "statistic",
    min = -6,
    max = 6,
    break_int = 4,
    title = " "
)
fig5 <- fig_a / fig_b / fig_c
fig5
ggsave(
    plot = fig5,
    filename = here("outputs", "figures", "sgmv_d_brain.svg"),
    height = 5,
    width = 5,
    dpi = 600
)