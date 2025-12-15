# Environment ------------------------------------------------------------------
library(tidyverse)
library(here)
library(glue)
library(tidymodels)
library(ggthemes)
library(patchwork)
library(lme4)


# Data I/O ---------------------------------------------------------------------
# ABCD data folder path
# using data release 5.1
abcd_folder_path <- "path_to_abcd_folder" %>%
    here("abcd-data-release-5.1/core")

demo_y_age <- here(abcd_folder_path, "abcd-general", "abcd_y_lt.csv") %>%
    read_csv(show_col_types = FALSE) %>%
    select(src_subject_id, eventname, interview_age)
mh_p_le <- here(abcd_folder_path, "mental-health", "mh_p_le.csv") %>%
    read_csv(show_col_types = FALSE)
mh_p_ksads_ptsd <- here(
    abcd_folder_path, "mental-health", "mh_p_ksads_ptsd.csv"
) %>%
    read_csv(show_col_types = FALSE)
mh_y_le <- here(
    abcd_folder_path, "mental-health", "mh_y_le.csv"
) %>%
    read_csv(show_col_types = FALSE)
mh_p_cbcl <- here(abcd_folder_path, "mental-health", "mh_p_cbcl.csv") %>%
    read_csv(show_col_types = FALSE) %>%
    filter(eventname != "baseline_year_1_arm_1")
master_fascores_df <- read_rds(
    file = here("outputs", "caches", "master_fascores_df.rds")
) %>%
    rename(interview_age_baseline = interview_age)

master_df <- mh_y_le %>%
    left_join(mh_p_le, by = c("src_subject_id", "eventname")) %>%
    left_join(mh_p_cbcl, by = c("src_subject_id", "eventname")) %>%
    left_join(master_fascores_df, by = "src_subject_id") %>%
    drop_na(dataset) %>%
    mutate(
        n_y_le = rowSums(pick(contains("past_yr_y")), na.rm = TRUE),
        n_p_le = rowSums(pick(contains("past_yr_p")), na.rm = TRUE)
    ) %>%
    select(
        src_subject_id, eventname, site_id_l, dataset,
        cbcl_scr_syn_totprob_t, n_y_le, n_p_le
    ) %>%
    drop_na()

discovery_df <- filter(master_df, dataset == "Discovery")
holdout_df <- filter(master_df, dataset == "Holdout")

discovery_std_df <- discovery_df %>%
    mutate(across(where(is.numeric), ~ scale(.))) %>%
    mutate(
        n_sle = (n_y_le + n_p_le) / 2
    )

ep_fit <- discovery_std_df %>%
    lmerTest::lmer(
        formula =
            cbcl_scr_syn_totprob_t ~ poly(n_sle, 1, raw = TRUE) +
            (1 + n_sle | src_subject_id) + (1 | site_id_l),
        control = lmerControl(optimizer = "nloptwrap", calc.derivs = FALSE)
    )

ep_fit_poly <- discovery_std_df %>%
    lmerTest::lmer(
        formula =
            cbcl_scr_syn_totprob_t ~ poly(n_sle, 2, raw = TRUE) +
            (1 + n_sle | src_subject_id) + (1 | site_id_l),
        control = lmerControl(optimizer = "nloptwrap", calc.derivs = FALSE)
    )
# does not have any improvement using quadratic
anova(ep_fit, ep_fit_poly)

broom.mixed::tidy(ep_fit)
r2_values <- performance::r2(ep_fit)
r2_values
correlation::correlation(
    data = discovery_std_df
)

discovery_std_df %>%
    group_by(src_subject_id) |> 
    mutate(
        totprob_sum = sum(cbcl_scr_syn_totprob_t),
        n_sle_sum = sum(n_sle)
    ) |> 
    slice(1) |> 
    ungroup() |> 
    lmerTest::lmer(
        formula =
            totprob_sum ~ poly(n_sle_sum, 1, raw = TRUE) + (1 | site_id_l),
        control = lmerControl(optimizer = "nloptwrap", calc.derivs = FALSE)
    ) |> 
    broom.mixed::tidy()


ep_fit_norandom <- discovery_std_df %>%
    lm(
        formula =
            cbcl_scr_syn_totprob_t ~ n_sle
    )

# use only fixed effects
discovery_srs_df <- discovery_std_df %>%
    mutate(
        pred_cbcl = predict(ep_fit, re.form = NA),
        srs = cbcl_scr_syn_totprob_t - pred_cbcl
    ) %>%
    group_by(src_subject_id) %>%
    mutate(
        avg_srs = mean(srs, na.rm = TRUE),
    ) %>%
    ungroup()

holdout_srs_df <- holdout_df %>%
    mutate(
        cbcl_scr_syn_totprob_t = scale(
            cbcl_scr_syn_totprob_t,
            center = attr(
                discovery_srs_df$cbcl_scr_syn_totprob_t, "scaled:center"
            ),
            scale = attr(
                discovery_srs_df$cbcl_scr_syn_totprob_t, "scaled:scale"
            )
        ),
        n_y_le = scale(
            n_y_le,
            center = attr(discovery_srs_df$n_y_le, "scaled:center"),
            scale = attr(discovery_srs_df$n_y_le, "scaled:scale")
        ),
        n_p_le = scale(
            n_p_le,
            center = attr(discovery_srs_df$n_p_le, "scaled:center"),
            scale = attr(discovery_srs_df$n_p_le, "scaled:scale")
        ),
        n_sle = (n_y_le + n_p_le) / 2
    ) %>%
    mutate(
        pred_cbcl = predict(ep_fit, newdata = ., re.form = NA),
        srs = cbcl_scr_syn_totprob_t - pred_cbcl
    ) %>%
    group_by(src_subject_id) %>%
    mutate(
        avg_srs = mean(srs, na.rm = TRUE)
    ) %>%
    ungroup()


# Plotting the relationship ----------------------------------------------------
ep_fit_labels <- c(
    paste("italic(beta) ==", round(broom.mixed::tidy(ep_fit)$estimate[2], 3)),
    paste("SE ==", round(broom.mixed::tidy(ep_fit)$std.error[2], 3)),
    paste("italic(p) ==", scales::label_scientific(digits = 3)(
        broom.mixed::tidy(ep_fit)$p.value[2])
    )
)

fig_a <- discovery_srs_df |> 
    ggplot(aes(x = n_sle, y = cbcl_scr_syn_totprob_t)) +
    geom_hex(binwidth = c(0.75, 0.75)) +
    geom_smooth(method = "lm", color = "tomato1", fill = "tomato3") +
    scale_fill_gradient(low = "gray20", high = "gray70") +
    labs(
        x = "Stressful Life Events (Z-score)",
        y = "CBCL Total Problems (Z-score)",
        fill = "",
        title = "A. SLE is correlated with CBCL Total Problems"
    ) +
    annotate(
        geom = "text",
        x = 4.5, 
        y = c(-1.25, -1.6, -1.9),
        label = ep_fit_labels,
        parse = TRUE,
        hjust = 0,
        vjust = 1,
        size = 5,
        family = "Helvetica Neue"
    ) +
    theme_pander() +
    theme(
        plot.margin = margin(5, 5, 5, 5, "mm"),
        legend.background = element_rect(fill = NA, color = NA),
        plot.title.position = "plot",
        text = element_text(family = "Helvetica Neue")
    )

fig_b <- discovery_srs_df |>
    group_by(src_subject_id) |> 
    slice(1) |> 
    ungroup() |> 
    ggplot(aes(x = avg_srs)) +
    geom_histogram(aes(y = after_stat(density)), binwidth = 0.25) +
    scale_x_continuous(
        labels = c("<- More Resilient", "0", "More Susceptible ->"),
        breaks = c(-1.8, 0, 2)
    ) +
    geom_density(adjust = 1.5, color = NA, fill = "gray30", alpha = 0.5) +
    theme_pander() +
    labs(
        x = "Stressor Reactivity Score (SRS)",
        y = "",
        title = "B. Distribution of SRS in Discovery Dataset"
    ) +
    theme(
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        plot.title.position = "plot",
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        text = element_text(family = "Helvetica Neue")
    )

fig_c <- holdout_srs_df |>
    group_by(src_subject_id) |> 
    slice(1) |> 
    ungroup() |> 
    ggplot(aes(x = avg_srs)) +
    geom_histogram(aes(y = after_stat(density)), binwidth = 0.2) +
    scale_x_continuous(
        labels = c("<- More Resilient", "0", "More Susceptible ->"),
        breaks = c(-1.8, 0, 2)
    ) +
    geom_density(adjust = 1.5, color = NA, fill = "gray50", alpha = 0.5) +
    theme_pander() +
    labs(
        x = "Stressor Reactivity Score (SRS)",
        y = "",
        title = "C. Distribution of SRS in Holdout Dataset"
    ) +
    theme(
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        plot.title.position = "plot",
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        text = element_text(family = "Helvetica Neue")
    )

fig <- free(fig_a) / (fig_b + fig_c) +
    plot_layout(heights = c(1.5, 1))
fig
ggsave(
    here("outputs", "figures", "srs_distribution.png"),
    fig,
    width = 10,
    height = 10,
    dpi = 600
)

baseline_age_df <- demo_y_age %>% 
    filter(eventname == "baseline_year_1_arm_1") %>% 
    rename(interview_age_baseline = interview_age) %>% 
    select(-eventname)

master_srs_df <- bind_rows(discovery_srs_df, holdout_srs_df) %>%
    left_join(demo_y_age, by = c("src_subject_id", "eventname")) %>% 
    left_join(
        master_fascores_df,
        by = c("src_subject_id", "site_id_l", "dataset")
    ) %>%
    mutate_if(is.numeric, as.numeric)
write_rds(
    x = master_srs_df,
    file = here("outputs", "caches", "master_srs_df.rds")
)
