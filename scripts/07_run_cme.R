# Environment ------------------------------------------------------------------
library(tidyverse)
library(here)
library(glue)
library(tidymodels)
library(ggthemes)
library(patchwork)
library(lme4)
library(lavaan)


# Data I/O ---------------------------------------------------------------------
env_labels <- read_rds(
    file = here("outputs", "caches", "discovery_sig_labels.rds")
)

master_srs_brain_df <- read_rds(
    here("outputs", "caches", "master_srs_brain_df.rds")
)


# Preprocessing ----------------------------------------------------------------
srs_env_brain_baseline_df <- master_srs_brain_df %>%
    filter(eventname == "baseline_year_1_arm_1") %>%
    mutate(
        pfc_volume =
            smri_vol_cdk_sufrlh + smri_vol_cdk_cdmdfrlh + 
            smri_vol_cdk_rracatelh + smri_vol_cdk_parstgrislh + 
            smri_vol_cdk_parsobislh + smri_vol_cdk_parsopclh +
            smri_vol_cdk_lobfrlh + smri_vol_cdk_mobfrlh + 
            smri_vol_cdk_frpolelh + smri_vol_cdk_rracatelh + 
            smri_vol_cdk_cdacatelh +
            smri_vol_cdk_sufrrh + smri_vol_cdk_cdmdfrrh + 
            smri_vol_cdk_rracaterh + smri_vol_cdk_parstgrisrh + 
            smri_vol_cdk_parsobisrh + smri_vol_cdk_parsopcrh +
            smri_vol_cdk_lobfrrh + smri_vol_cdk_mobfrrh + 
            smri_vol_cdk_frpolerh + smri_vol_cdk_rracaterh + 
            smri_vol_cdk_cdacaterh,
        amy_volume = smri_vol_scs_amygdalalh + smri_vol_scs_amygdalarh,
        vs_volume = smri_vol_scs_caudatelh + smri_vol_scs_caudaterh +
            smri_vol_scs_putamenlh + smri_vol_scs_putamenrh +
            smri_vol_scs_aal + smri_vol_scs_aar,
        hc_volume = smri_vol_scs_hpuslh + smri_vol_scs_hpusrh
    ) |> 
    select(
        src_subject_id, dataset, demo_sex_v2, race_ethnicity_3, site_id_l, 
        eventname, avg_srs, interview_age, `social services`:`family SES`, 
        pfc_volume, amy_volume, vs_volume, hc_volume, smri_vol_cdk_total,
        smri_vol_scs_subcorticalgv, smri_vol_scs_intracranialv
    )


# Function ---------------------------------------------------------------------
run_causal_mediation <- function(
    df, 
    exposure_vars, 
    mediator_var, 
    outcome_var = "avg_srs",
    covars,
    random_effect = "site_id_l",
    subset_type,
    standardize = TRUE,
    sims = 1000
) {
    
    # 1. Filter Data based on subset_type
    if (subset_type == "Discovery") {
        data_subset <- df %>% filter(dataset == "Discovery")
        message("--- Running on Discovery dataset ---")
    } else if (subset_type == "Holdout") {
        data_subset <- df %>% filter(dataset == "Holdout")
        message("--- Running on HOLDOUT dataset ---")
    } else if (subset_type == "Whole") {
        data_subset <- df
        message("--- Running on FULL dataset ---")
    } else {
        stop("Please indicate your dataset")
    }
    
    # 2. Combine covariates and random effects into strings
    covars_formula <- paste(covars, collapse = " + ")
    re_formula <- paste0("(1 | ", random_effect, ")")
    
    # 3. Loop through each exposure variable
    results <- purrr::map_dfr(exposure_vars, function(env) {
        
        message(glue::glue("Analyzing: Exposure='{env}' -> Mediator='{mediator_var}'"))
        
        # --- Step A: Standardization (Optional) ---
        # We do this inside the loop to ensure we only scale the specific vars used in this iteration
        analysis_data <- data_subset
        
        if (standardize) {
            # Identify all variables involved in this specific model
            vars_in_model <- c(outcome_var, mediator_var, env, covars)
            
            # Find which of these are numeric (to avoid scaling Factors/Categorical vars)
            numeric_cols <- analysis_data %>%
                select(all_of(vars_in_model)) %>%
                select(where(is.numeric)) %>%
                names()
            
            # Apply Z-score scaling
            # We use as.vector() to strip attributes, preventing issues with lme4
            analysis_data <- analysis_data %>%
                mutate(across(all_of(numeric_cols), ~ as.vector(scale(.))))
        }
        
        # --- Step B: Build Formulas ---
        # Model M: Mediator ~ Exposure + Covariates + (1|Site)
        f_mediator <- as.formula(glue::glue("{mediator_var} ~ `{env}` + {covars_formula} + {re_formula}"))
        
        # Model Y: Outcome ~ Exposure + Mediator + Covariates + (1|Site)
        f_outcome <- as.formula(glue::glue("{outcome_var} ~ `{env}` + {mediator_var} + {covars_formula} + {re_formula}"))
        
        # --- Step C: Fit Mixed Models (lmer) ---
        tryCatch({
            model_m <- lme4::lmer(
                f_mediator, 
                data = analysis_data, 
                control = lmerControl(
                    optimizer = "nloptwrap", calc.derivs = FALSE
                )
            )
            model_y <- lme4::lmer(
                f_outcome, 
                data = analysis_data, 
                control = lmerControl(
                    optimizer = "nloptwrap", calc.derivs = FALSE
                )
            )
            
            # --- Step D: Run Causal Mediation (Simulations) ---
            med_out <- mediation::mediate(
                model.m = model_m,
                model.y = model_y,
                treat = env,
                mediator = mediator_var,
                sims = sims
            )
            
            # --- Step E: Extract Results ---
            med_summary <- summary(med_out)
            
            tibble(
                exposure = env,
                mediator = mediator_var,
                outcome = outcome_var,
                type = if_else(standardize, "Standardized", "Unstandardized"),
                effect_type = c("ACME (Average Causal Mediation)", "ADE (Average Direct Effect)", "Prop. Mediated"),
                estimate = c(med_summary$d0, med_summary$z0, med_summary$n0),
                p_value = c(med_summary$d0.p, med_summary$z0.p, med_summary$n0.p),
                ci_lower = c(med_summary$d0.ci[1], med_summary$z0.ci[1], med_summary$n0.ci[1]),
                ci_upper = c(med_summary$d0.ci[2], med_summary$z0.ci[2], med_summary$n0.ci[2])
            )
            
        }, error = function(e) {
            message(glue::glue("Error for {env}: {e$message}"))
            return(tibble()) # Return empty if model fails
        })
    })
    
    return(results)
}

plot_validation_forest <- function(
    discovery_res,
    holdout_res,
    effect_name = "ACME",
    plot_title = NULL,
    x_limit = NULL
) {
    
    # --- Step 1: Process Holdout ---
    # identify which variables are significant (FDR < 0.05) in the Holdout set
    holdout_sig_status <- holdout_res |> 
        filter(str_detect(effect_type, effect_name)) |> 
        mutate(padj = p.adjust(p_value, "fdr")) |> 
        mutate(
            is_replicated = if_else(padj < 0.05, "Yes", "No"),
            is_replicated = factor(is_replicated, levels = c("Yes", "No"))
        ) |> 
        select(exposure, is_replicated)
    
    # --- Step 2: Process Discovery ---
    plot_data <- discovery_res |> 
        filter(str_detect(effect_type, effect_name)) |> 
        mutate(padj = p.adjust(p_value, "fdr")) |> 
        mutate(is_significant_discovery = padj < 0.05) |> 
        left_join(holdout_sig_status, by = "exposure")
    
    # --- Step 3: Generate Plot ---
    p <- plot_data |> 
        ggplot(aes(
            x = estimate, 
            y = exposure,
            color = is_replicated,
            alpha = is_significant_discovery
        )) +
        geom_vline(xintercept = 0, color = "grey70", linetype = "dashed") +
        geom_errorbarh(aes(xmin = ci_lower, xmax = ci_upper), height = 0.2) +
        geom_point(size = 3) +
        scale_y_discrete(limits = rev) +
        scale_color_manual(
            values = c("Yes" = "tomato2", "No" = "gray30"),
            limits = c("Yes", "No"),
            drop = FALSE
        ) +
        scale_alpha_manual(
            values = c(`TRUE` = 1, `FALSE` = 0.3), guide = "none"
        ) +
        labs(
            title = plot_title,
            x = "Standardized Estimate (Beta)",
            y = NULL,
            color = "Replicated in Holdout?"
        ) +
        theme_minimal() +
        theme(
            legend.position = "bottom",
            panel.grid.minor = element_blank(),
            plot.title = element_text(hjust = 0.5),
            axis.text.y = element_text(size = 10, color = "black")
        )
    
    # Apply X limits if provided
    if (!is.null(x_limit)) {
        p <- p + scale_x_continuous(limits = x_limit)
    }
    
    return(p)
}

extract_to_doc <- function(df, effect_name, filename = NULL) {
    
    sect_properties <- prop_section(
        page_size = page_size(
            orient = "portrait",
            width = 8.3, height = 11.7
        ),
        type = "continuous",
        page_margins = page_mar(gutter = 0)
    )
    
    extracted_df <- df |> 
        filter(str_detect(effect_type, effect_name)) |> 
        mutate(padj = p.adjust(p_value, method = "fdr")) |> 
        select(exposure, estimate, ci_lower, ci_upper, p_value, padj) |> 
        arrange(desc(abs(estimate))) |> 
        mutate_if(is.numeric, ~ round(., 3)) |> 
        rename_with(
            ~ c("Environmental Factors", "Estimate", "Lower 95% CI",
                "Upper 95% CI", "p", "p_fdr")
        )
    
    if (!is.null(filename)) {
        
        extracted_df |> 
            flextable() |> 
            autofit() |> 
            bold(part = "header") |> 
            compose(
                j = "p_fdr", 
                part = "header", 
                value = as_paragraph("p", as_sub("fdr"))
            ) |>
            align(j = 2:6, align = "center", part = "all") |> 
            save_as_docx(
                path = here("outputs", "tables", paste0(filename, ".docx")),
                pr_section = sect_properties
            )
    }
    
    return(extracted_df) 
}


# Causal Mediation Analysis ----------------------------------------------------
covariates_global <- c(
    "interview_age",
    "demo_sex_v2",
    "race_ethnicity_3"
)

cme_d_cgmv <- run_causal_mediation(
    df = srs_env_brain_baseline_df,
    exposure_vars = env_labels,
    covars = covariates_global,
    subset_type = "Discovery",
    mediator_var = "smri_vol_cdk_total"
)
cme_h_cgmv <- run_causal_mediation(
    df = srs_env_brain_baseline_df,
    exposure_vars = env_labels,
    covars = covariates_global,
    subset_type = "Holdout",
    mediator_var = "smri_vol_cdk_total"
)
global_a <- plot_validation_forest(
    cme_d_cgmv, cme_h_cgmv, x_limit = c(-0.015, 0.015),
    plot_title = "Cortical GMV"
)

extract_to_doc(cme_d_cgmv, effect_name = "ACME", filename = "cme_d_cgmv")
extract_to_doc(cme_h_cgmv, effect_name = "ACME", filename = "cme_h_cgmv")

cme_d_sgmv <- run_causal_mediation(
    df = srs_env_brain_baseline_df,
    exposure_vars = env_labels,
    covars = covariates_global,
    subset_type = "Discovery",
    mediator_var = "smri_vol_scs_subcorticalgv"
)
cme_h_sgmv <- run_causal_mediation(
    df = srs_env_brain_baseline_df,
    exposure_vars = env_labels,
    covars = covariates_global,
    subset_type = "Holdout",
    mediator_var = "smri_vol_scs_subcorticalgv"
)

extract_to_doc(cme_d_sgmv, effect_name = "ACME", filename = "cme_d_sgmv")
extract_to_doc(cme_h_sgmv, effect_name = "ACME", filename = "cme_h_sgmv")

global_b <- plot_validation_forest(
    cme_d_sgmv, cme_h_sgmv, x_limit = c(-0.015, 0.015),
    plot_title = "Subcortical GMV"
) +
    theme(
        axis.text.y = element_blank(),
        legend.position = "None"
    )

total_gmv_plot <- (global_a + global_b)
total_gmv_plot
ggsave(
    plot = total_gmv_plot,
    filename = here("outputs", "figures", "total_gmv_plot.pdf"),
    device = cairo_pdf,
    width = 8,
    height = 5.5
)


covariates_region <- c(
    "interview_age",
    "demo_sex_v2",
    "race_ethnicity_3",
    "smri_vol_scs_intracranialv"
)

cme_d_pfc <- run_causal_mediation(
    df = srs_env_brain_baseline_df,
    covars = covariates_region,
    exposure_vars = env_labels,
    subset_type = "Discovery",
    mediator_var = "pfc_volume"
)

extract_to_doc(cme_d_pfc, effect_name = "ACME", filename = "cme_d_pfc")

cme_h_pfc <- run_causal_mediation(
    df = srs_env_brain_baseline_df,
    covars = covariates_region,
    exposure_vars = env_labels,
    subset_type = "Holdout",
    mediator_var = "pfc_volume"
)

extract_to_doc(cme_h_pfc, effect_name = "ACME", filename = "cme_h_pfc")

plot_validation_forest(cme_d_pfc, cme_h_pfc)
    
cme_d_hc <- run_causal_mediation(
    df = srs_env_brain_baseline_df,
    covars = covariates_region,
    exposure_vars = env_labels,
    subset_type = "Discovery",
    mediator_var = "hc_volume"
)

cme_h_hc <- run_causal_mediation(
    df = srs_env_brain_baseline_df,
    covars = covariates_region,
    exposure_vars = env_labels,
    subset_type = "Holdout",
    mediator_var = "hc_volume"
)

extract_to_doc(cme_d_hc, effect_name = "ACME", filename = "cme_d_hc")
extract_to_doc(cme_h_hc, effect_name = "ACME", filename = "cme_h_hc")

plot_validation_forest(cme_d_hc, cme_h_hc)

cme_d_amy <- run_causal_mediation(
    df = srs_env_brain_baseline_df,
    covars = covariates_region,
    exposure_vars = env_labels,
    subset_type = "Discovery",
    mediator_var = "amy_volume"
)

cme_h_amy <- run_causal_mediation(
    df = srs_env_brain_baseline_df,
    covars = covariates_region,
    exposure_vars = env_labels,
    subset_type = "Holdout",
    mediator_var = "amy_volume"
)

extract_to_doc(cme_d_amy, effect_name = "ACME", filename = "cme_d_amy")
extract_to_doc(cme_h_amy, effect_name = "ACME", filename = "cme_h_amy")

plot_validation_forest(cme_d_amy, cme_h_amy)

cme_d_vs <- run_causal_mediation(
    df = srs_env_brain_baseline_df,
    covars = covariates_region,
    exposure_vars = env_labels,
    subset_type = "Discovery",
    mediator_var = "vs_volume"
)

cme_h_vs <- run_causal_mediation(
    df = srs_env_brain_baseline_df,
    covars = covariates_region,
    exposure_vars = env_labels,
    subset_type = "Holdout",
    mediator_var = "vs_volume"
)

extract_to_doc(cme_d_vs, effect_name = "ACME", filename = "cme_d_vs")
extract_to_doc(cme_h_vs, effect_name = "ACME", filename = "cme_h_vs")

plot_validation_forest(cme_d_vs, cme_h_vs)