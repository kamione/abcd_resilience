# Environment ------------------------------------------------------------------
library(tidyverse)
library(here)
library(glue)
library(tidymodels)
library(ggthemes)
library(patchwork)
library(ggtext)
library(gtsummary)
library(flextable)
library(officer)

sect_properties <- prop_section(
    page_size = page_size(
        orient = "landscape",
        width = 8.3, height = 11.7
    ),
    type = "continuous",
    page_margins = page_mar(gutter = 0)
)

sect_properties_portrait <- prop_section(
    page_size = page_size(
        orient = "portrait",
        width = 11.7, height = 8.3
    ),
    type = "continuous",
    page_margins = page_mar(gutter = 0)
)

# Data I/O ---------------------------------------------------------------------
master_srs_df <- read_rds(
    file = here("outputs", "caches", "master_srs_df.rds")
) %>%
    group_by(src_subject_id) %>%
    slice(1) %>%
    ungroup() %>% 
    select(-eventname)

table_dataset_des <- master_srs_df %>% 
    select(dataset, interview_age_baseline, demo_sex_v2, race_ethnicity_3, 
           avg_srs) %>% 
    tbl_summary(
        by = dataset,
        statistic = list(all_continuous() ~ "{mean} ({sd})"),
        digits = all_continuous() ~ 2,
        label = list(
            interview_age_baseline ~ "Age at Baseline (months)",
            demo_sex_v2 ~ "Sex at Birth",
            race_ethnicity_3 ~ "Race/Ethnicity",
            avg_srs ~ "Average SR Scores"
        )
    ) %>% 
    add_overall() %>% 
    add_p() %>% 
    add_q() %>% 
    as_flex_table() %>% 
    autofit()
table_dataset_des
table_dataset_des %>% 
    save_as_docx(
        path = here("outputs", "tables", "table_dataset_desc.docx"),
        pr_section = sect_properties
    )

discovery_df <- filter(master_srs_df, dataset == "Discovery")
holdout_df <- filter(master_srs_df, dataset == "Holdout")

fa_labels <- discovery_df %>%
    select(`social services`:`family SES`) %>%
    colnames()

discovery_results <- list()
for (ith in seq_along(fa_labels)) {
    formula_text <- glue::glue(
        "avg_srs ~ `{fa_labels[ith]}` + interview_age_baseline + demo_sex_v2 +
            race_ethnicity_3 + (1 | site_id_l)"
    )
    discovery_results[[ith]] <- lmerTest::lmer(
        formula = as.formula(formula_text),
        data = discovery_df
    ) %>%
        broom.mixed::tidy(conf.int = TRUE) %>%
        slice(2)
}
discovery_results_df <- reduce(discovery_results, bind_rows) %>%
    mutate(p.adj = p.adjust(p.value, method = "fdr")) %>% 
    mutate(term = str_replace_all(term, "`", "")) %>% 
    select(3:5, 9, 10, 6, 8, 11)

set_e_notation <- function(x) {
    formatC(x, format = "e", digits = 2)
}

univariate_reg_results_d <- discovery_results_df %>% 
    arrange(desc(abs(estimate))) %>% 
    rename_all(~c("Environmental Factor", "β", "SE", "5% CI", "95% CI", "t-value", "p", "p(FDR)")) %>% 
    flextable() %>% 
    colformat_double(
        j = 2:5,
        digits = 3
    ) %>% 
    colformat_double(
        j = 6,
        digits = 2
    ) %>% 
    set_formatter(p = set_e_notation) %>% 
    set_formatter(`p(FDR)` = set_e_notation) %>%
    align(align = "center", part = "header") %>% 
    bold(part = "header") %>% 
    italic(j = c(2, 7, 8), italic = TRUE, part = "header") %>% 
    fontsize(size = 8, part = "all") %>% 
    autofit()
univariate_reg_results_d 

univariate_reg_results_d %>% save_as_docx(
    path = here("outputs", "tables", "univariate_reg_results_d.docx"),
    pr_section = sect_properties_portrait
)


discovery_sig_labels <- discovery_results_df %>%
    filter(p.adj < 0.05) %>%
    arrange(desc(abs(estimate))) %>% 
    pull(term)

write_rds(
    discovery_sig_labels,
    here("outputs", "caches", "discovery_sig_labels.rds")
)

results_holdout <- list()
for (ith in seq_along(fa_labels)) {
    formula_text <- glue::glue(
        "avg_srs ~ `{fa_labels[ith]}` + interview_age_baseline + demo_sex_v2 +
            race_ethnicity_3 + (1 | site_id_l)"
    )
    results_holdout[[ith]] <- lmerTest::lmer(
        formula = as.formula(formula_text),
        data = holdout_df
    ) %>%
        broom.mixed::tidy(conf.int = TRUE) %>%
        slice(2)
}
holdout_results_df <- reduce(results_holdout, bind_rows) %>%
    mutate(p.adj = p.adjust(p.value, method = "fdr")) %>%
    mutate(term = str_replace_all(term, "`", "")) %>% 
    select(3:5, 9, 10, 6, 8, 11)

univariate_reg_results_h <- holdout_results_df %>% 
    arrange(desc(abs(estimate))) %>% 
    rename_all(~c("Environmental Factor", "β", "SE", "5% CI", "95% CI", "t-value", "p", "p(FDR)")) %>% 
    flextable() %>% 
    colformat_double(
        j = 2:5,
        digits = 3
    ) %>% 
    colformat_double(
        j = 6,
        digits = 2
    ) %>% 
    set_formatter(p = set_e_notation) %>% 
    set_formatter(`p(FDR)` = set_e_notation) %>%
    align(align = "center", part = "header") %>% 
    bold(part = "header") %>% 
    italic(j = c(2, 7, 8), italic = TRUE, part = "header") %>% 
    fontsize(size = 8, part = "all") %>% 
    autofit()
univariate_reg_results_h

univariate_reg_results_h %>% save_as_docx(
    path = here("outputs", "tables", "univariate_reg_results_h.docx"),
    pr_section = sect_properties_portrait
)

holdout_falabels <- holdout_results_df %>% 
    filter(p.adj < 0.05) %>%
    pull(term)

intersect(holdout_falabels, discovery_sig_labels)

univariate_results_fig <- discovery_results_df %>%
    left_join()
    mutate(
        sig_holdout = if_else(term %in% holdout_falabels, TRUE, FALSE)
    ) %>%
    ggplot(aes(x = term, y = estimate)) +
    geom_point(
        aes(size = -log(p.adj)),
        color = "gray40"
    ) +
    geom_point(
        data = . %>% filter(sig_holdout == TRUE),
        aes(size = -log(p.adj) + 10),
        color = "tomato3",
        shape = 1
    ) +
    scale_x_discrete(labels = function(x) gsub("`", "", x)) +
    scale_y_continuous(
        limits = c(-0.41, 0.45),
        breaks = seq(-0.5, 0.5, 0.1),
        sec.axis = sec_axis(
            ~.,  # Right y-axis mirroring
            breaks = c(-0.25, 0.25), 
            labels = c("<-- more resilient", "more susceptible -->")
        )
    ) +
    scale_size_continuous(
        range = c(1, 7),
        breaks = c(2, 20, 70),
        labels = c("p < 0.05", "p < 0.001", "p < 1e-10")
    ) +
    geom_hline(yintercept = 0, color = "gray70") +
    labs(
        x = "",
        y = "standardized *<span>&#x03B2;</span>* coefficient",
        size = ""
    ) +
    gghighlight::gghighlight(
        p.adj < 0.05,
        use_direct_label = FALSE
    ) +
    ggthemes::theme_pander() +
    theme(
        axis.text.x = element_text(angle = 60, hjust = 1),
        axis.title.y = element_markdown(),
        axis.text.y.right = element_markdown(angle = 90, hjust = 0.5),
        plot.margin = margin(10, 10, 10, 10, "mm"),
        legend.text = element_text(size = 12),
        legend.position = c(0.07, 0.9),
        legend.background = element_rect(fill = NA, color = NA),
        text = element_text(family = "Helvetica Neue")
    )
univariate_results_fig
ggsave(
    plot = univariate_results_fig,
    filename = here("outputs", "figures", "univariate_results_fig.pdf"),
    height = 6,
    width = 10,
    device = cairo_pdf
)
