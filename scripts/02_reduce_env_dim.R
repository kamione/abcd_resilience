# Environment ------------------------------------------------------------------
library(here)
library(tidyverse)
library(ggpubr)

# Data I/O ---------------------------------------------------------------------
master_cleaned_imputated_df <- read_rds(
    here("data", "processed", "master_cleaned_imputated_df.rds")
)
# exposome variables
exposome_variables_df <- here("data", "raw", "included_variables.xlsx") %>%
    readxl::read_excel(sheet = "Exposome Raw")


# Helper Function --------------------------------------------------------------
.gg_parallel <- function(parallel_results) {
    
    # 1. Extract the relevant vectors from the psych object
    df_list <- list()
    if (!is.null(parallel_results$fa.values)) {
        df_list$fa_obs <- data.frame(
            x = seq_along(parallel_results$fa.values), 
            value = parallel_results$fa.values, 
            type = "Observed" 
        )
    }
    if (!is.null(parallel_results$fa.sim)) {
        df_list$fa_sim <- data.frame(
            x = seq_along(parallel_results$fa.sim), 
            value = parallel_results$fa.sim, 
            type = "Simulated (95th Percentile)"
        )
    }
    
    plot_data <- dplyr::bind_rows(df_list)
    
    # 2. Create the Plot
    p <- plot_data |> 
        ggplot(aes(x = x, y = value, group = type, color = type, linetype = type)) +
        # Add Kaiser criterion line (eigenvalue = 1)
        geom_hline(yintercept = 1, linetype = "dotted", color = "gray50") +
        geom_line(size = 0.8) +
        geom_point(size = 2) +
        scale_color_manual(values = c(
            "Observed" = "grey30",   
            "Simulated (95th Percentile)" = "grey70"
        )) +
        scale_linetype_manual(values = c(
            "Observed" = "solid", 
            "Simulated (95th Percentile)" = "dashed"
        )) +
        scale_x_continuous(breaks = function(x) unique(floor(pretty(x)))) +
        labs(
            title = "A. Parallel Analysis",
            x = "Number of Factor",
            y = "Eigenvalue",
            color = NULL,
            linetype = NULL,
        ) +
        theme_minimal() +
        theme(
            legend.position = "bottom",
            plot.title.position = "plot",
            plot.title = element_text(face = "bold"),
        )
    
    return(p)
}

.gg_heatmap_loadings <- function(fa_res, fa_labels) {
    loadings_obj <- fa_res$loadings
    loadings_df <- as.data.frame(unclass(loadings_obj)) |>
        rename_all(~ c(fa_labels)) |> 
        tibble::rownames_to_column("Variable") |> 
        pivot_longer(cols = -Variable, names_to = "Factor", values_to = "Loading")
    
    p <- loadings_df |> 
        ggplot(aes(x = Factor, y = Variable, fill = Loading)) +
        geom_tile(color = "white") +
        geom_text(aes(label = round(Loading, 2)), color = "gray30", size = 4) +
        scale_fill_gradient2(
            low = "royalblue", mid = "white", high = "tomato3", 
            midpoint = 0, limit = c(-1.15, 1.15)
        ) +
        labs(title = "B. Factor Loadings", x = NULL, y = NULL) +
        theme_minimal() +
        theme(
            plot.title.position = "plot",
            plot.title = element_text(face = "bold"),
            axis.text.x = element_text(size = 12, angle = 45, hjust = 1),
            axis.text.y = element_text(size = 12),
            panel.grid = element_blank()
        )
    
    return(p)
}

.gg_heatmap_lower_tri <- function(fa_res, fa_labels) {
    
    phi_masked <- fa_res$Phi
    phi_masked[upper.tri(phi_masked, diag = FALSE)] <- NA
    
    colnames(phi_masked) <- fa_labels
    rownames(phi_masked) <- fa_labels
    
    cor_df <- as.data.frame(unclass(phi_masked)) |> 
        rownames_to_column("Factor1") |> 
        pivot_longer(cols = -Factor1, names_to = "Factor2", values_to = "r") |> 
        filter(!is.na(r))
    
    factor_names <- colnames(phi_masked)
    cor_df$Factor1 <- factor(cor_df$Factor1, levels = rev(factor_names))
    cor_df$Factor2 <- factor(cor_df$Factor2, levels = factor_names)

    p <- cor_df |> 
        ggplot(aes(x = Factor2, y = Factor1, fill = r)) +
        geom_tile(color = "white", width = 0.95, height = 0.95) +
        geom_text(aes(label = round(r, 2)), color = "black", size = 4) +
        scale_fill_gradient2(
            low = "royalblue", mid = "white", high = "tomato3", 
            midpoint = 0, limit = c(-1.15, 1.15)
        ) +
        # Formatting
        labs(title = "C. Factor Correlations", x = NULL, y = NULL, fill = "Phi") +
        theme_minimal() +
        theme(
            plot.title.position = "plot",
            plot.title = element_text(face = "bold"),
            axis.text.x = element_text(size = 12, angle = 45, hjust = 1),
            axis.text.y = element_text(size = 12),
            panel.grid = element_blank()
        )
    
    return(p)
}

.plot_fa <- function(parallel_results, fa_res, fa_labels) {
    
    p_parallel <- patchwork::wrap_elements(.gg_parallel(parallel_results))
    
    if (length(fa_labels) > 1) {
        p_bottom <- .gg_heatmap_loadings(fa_res, fa_labels) + 
            .gg_heatmap_lower_tri(fa_res, fa_labels)
        
        p <- free(p_parallel) / p_bottom
        
    } else {
        p <- free(p_parallel) / .gg_heatmap_loadings(fa_res, fa_labels)
    }
    
    return(p)
}


# Dimensionality Reduction -----------------------------------------------------
# using factor analysis to reduce to the minimal numbers of scores for each
# category
discovery_df <- master_cleaned_imputated_df %>%
    filter(dataset == "Discovery")
holdout_df <- master_cleaned_imputated_df %>%
    filter(dataset == "Holdout")

discovery_fascores_df <- discovery_df
holdout_fascores_df <- holdout_df

n_obs <- dim(discovery_df)[1]
n_vars <- dim(discovery_df)[2]
cat(
    "There are", n_obs, "participants and", n_vars,
    "variables in the discovery dataset"
)

# 1. Amenities and Services
amenities_vars <- exposome_variables_df %>%
    filter(category == "Amenities and Services") %>%
    pull(variable_name)

amenities_raw <- discovery_df %>%
    select(any_of(amenities_vars)) %>%
    mutate_all(as.numeric)
amenities_xcor <- amenities_raw %>%
    psych::mixedCor(ncat = 5) %>%
    pluck("rho")

amenities_par <- psych::fa.parallel(amenities_xcor, fa = "fa", n.obs = n_obs, n.iter = 1000)
amenities_par
for (ith in 2:4) {
    psych::fa.sort(
        psych::fa(amenities_xcor, nfactors = ith, rotate = "promax")
    ) %>%
        pluck("loadings") %>%
        print(cutoff = 0.2)
}

amenities_res <- psych::fa.sort(
    psych::fa(amenities_raw, nfactors = 4, rotate = "promax")
)
amenities_fa_labels <- c(
    "social services",
    "public organizations",
    "park availability",
    "no internet access"
)
discovery_fascores_df <- discovery_fascores_df %>%
    bind_cols(
        amenities_res$scores %>%
            as_tibble() %>%
            set_names(amenities_fa_labels)
    )
holdout_fascores_df <- holdout_fascores_df %>%
    bind_cols(
        holdout_df %>%
            select(any_of(amenities_vars)) %>%
            mutate_all(as.numeric) %>%
            predict(amenities_res, data = .) %>%
            as_tibble() %>%
            set_names(amenities_fa_labels)
    )

amenities_fa_fig <- .plot_fa(amenities_par, amenities_res, amenities_fa_labels)
ggsave(
    plot = amenities_fa_fig,
    filename = here("outputs", "figures", "amenities_fa_fig.png"),
    height = 8,
    width = 10,
    dpi = 600
)

# Community Health Care
comhealth_vars <- exposome_variables_df %>%
    filter(category == "Community Health Care") %>%
    pull(variable_name)

comhealth_raw <- discovery_df %>%
    select(any_of(comhealth_vars)) %>%
    mutate_all(as.numeric) # 29
comhealth_xcor <- comhealth_raw %>%
    psych::mixedCor(ncat = 5) %>%
    pluck("rho")

comhealth_par <- psych::fa.parallel(comhealth_xcor, fa = "fa", n.obs = n_obs, n.iter = 1000)
for (ith in 7:9) {
    psych::fa.sort(
        psych::fa(comhealth_xcor, nfactors = ith, rotate = "promax")
    ) %>%
        pluck("loadings") %>%
        print(cutoff = 0.5)
}

comhealth_res <- psych::fa.sort(
    psych::fa(comhealth_raw, nfactors = 7, rotate = "promax")
)
comhealth_fa_labels <- c(
    "healthcare access", "community chronic conditions", "poorer community wellbeing",
    "preventive screening", "healthcare service availability",
    "elderly care", "binge drinking population"
)
discovery_fascores_df <- discovery_fascores_df %>%
    bind_cols(
        comhealth_res$scores %>%
            as_tibble() %>%
            set_names(comhealth_fa_labels)
    )
holdout_fascores_df <- holdout_fascores_df %>%
    bind_cols(
        holdout_df %>%
            select(any_of(comhealth_vars)) %>%
            mutate_all(as.numeric) %>%
            predict(comhealth_res, data = .) %>%
            as_tibble() %>%
            set_names(comhealth_fa_labels)
    )

comhealth_fa_fig <- .plot_fa(comhealth_par, comhealth_res, comhealth_fa_labels)
ggsave(
    plot = comhealth_fa_fig,
    filename = here("outputs", "figures", "comhealth_fa_fig.png"),
    height = 14,
    width = 16,
    dpi = 600
)

# 3. Developmental History
devhist_raw <- discovery_df %>%
    select(matches("devhx"), birth_weight_kg, ksads_ptsd_sum) %>%
    mutate_all(as.numeric) # 16
devhist_xcor <- devhist_raw %>%
    psych::mixedCor(ncat = 3) %>%
    pluck("rho")

devhist_par <- psych::fa.parallel(devhist_xcor, fa = "fa", n.obs = n_obs, n.iter = 1000)
devhist_par
for (ith in 6:8) {
    psych::fa.sort(psych::fa(devhist_raw, ith, rotate = "promax")) %>%
        pluck("loadings") %>%
        print(cutoff = 0.1)
}

devhist_res <- psych::fa.sort(psych::fa(devhist_raw, 8, rotate = "promax"))

devhist_fa_labels <- c(
    "parental ages at childbirth", "birth weight",
    "developmental delays", "perinatal conditions", "pregnancy condition",
    "prenatal & postnatal care", "prenatal substance use",
    "prenatal & postnatal adversity"
)
discovery_fascores_df <- discovery_fascores_df %>%
    bind_cols(
        devhist_res$scores %>%
            as_tibble() %>%
            set_names(devhist_fa_labels)
    )
holdout_fascores_df <- holdout_fascores_df %>%
    bind_cols(
        holdout_df %>%
            select(matches("devhx"), birth_weight_kg, ksads_ptsd_sum) %>%
            mutate_all(as.numeric) %>%
            predict(devhist_res, data = .) %>%
            as_tibble() %>%
            set_names(devhist_fa_labels)
    )

devhist_fa_fig <- .plot_fa(devhist_par, devhist_res, devhist_fa_labels)
ggsave(
    plot = devhist_fa_fig,
    filename = here("outputs", "figures", "devhist_fa_fig.png"),
    height = 14,
    width = 16,
    dpi = 600
)

# Family Values
famvalue_vars <- exposome_variables_df %>%
    filter(category == "Family Values") %>%
    pull(variable_name)

famvalue_raw <- discovery_df %>%
    select(any_of(famvalue_vars)) %>%
    mutate_all(as.numeric) # 13
famvalue_xcor <- famvalue_raw %>%
    psych::mixedCor(ncat = 3) %>%
    pluck("rho")

famvalue_par <- psych::fa.parallel(famvalue_xcor, fa = "fa", n.obs = n_obs, n.iter = 1000)
for (ith in 4:6) {
    psych::fa.sort(
        psych::fa(
            famvalue_xcor, nfactors = ith, rotate = "promax"
        )
    ) %>%
        pluck("loadings") %>%
        print(cutoff = 0.2)
}

famvalue_res <- psych::fa.sort(
    psych::fa(famvalue_raw, nfactors = 4, rotate = "promax")
)
famvalue_fa_labels <- c(
    "substance accessibility", "substance rules",
    "parental monitoring & support", "alcohol accessibility"
)
discovery_fascores_df <- discovery_fascores_df %>%
    bind_cols(
        famvalue_res$scores %>%
            as_tibble() %>%
            set_names(famvalue_fa_labels)
    )
holdout_fascores_df <- holdout_fascores_df %>%
    bind_cols(
        holdout_df %>%
            select(any_of(famvalue_vars)) %>%
            mutate_all(as.numeric) %>%
            predict(famvalue_res, data = .) %>%
            as_tibble() %>%
            set_names(famvalue_fa_labels)
    )

famvalue_fa_fig <- .plot_fa(famvalue_par, famvalue_res, famvalue_fa_labels)
ggsave(
    plot = famvalue_fa_fig,
    filename = here("outputs", "figures", "famvalue_fa_fig.png"),
    height = 12,
    width = 14,
    dpi = 600
)

 # Laws and Policies
law_vars <- exposome_variables_df %>%
    filter(category == "Laws and Policies") %>%
    pull(variable_name)

law_raw <- discovery_df %>%
    select(any_of(law_vars)) %>%
    mutate_all(as.numeric) # 13
law_xcor <- law_raw %>%
    psych::mixedCor(ncat = 3) %>%
    pluck("rho")

law_par <- psych::fa.parallel(law_xcor, fa = "fa", n.obs = n_obs, n.iter = 1000)
for (ith in 4:7) {
    psych::fa.sort(psych::fa(law_xcor, nfactors = ith, rotate = "promax")) %>%
        pluck("loadings") %>%
        print(cutoff = 0.2)
}

law_res <- psych::fa.sort(psych::fa(law_raw, nfactors = 5, rotate = "promax"))
law_fa_labels <- c(
    "drug policy & reforms", "naloxone access & regulations",
    "Prescription opioid safety policies", "Good Samaritan law",
    "drug policy environment"
)
discovery_fascores_df <- discovery_fascores_df %>%
    bind_cols(
        law_res$scores %>%
            as_tibble() %>%
            set_names(law_fa_labels)
    )
holdout_fascores_df <- holdout_fascores_df %>%
    bind_cols(
        holdout_df %>%
            select(any_of(law_vars)) %>%
            mutate_all(as.numeric) %>%
            predict(law_res, data = .) %>%
            as_tibble() %>%
            set_names(law_fa_labels)
    )

law_fa_fig <- .plot_fa(law_par, law_res, law_fa_labels)
ggsave(
    plot = law_fa_fig,
    filename = here("outputs", "figures", "law_fa_fig.png"),
    height = 12,
    width = 14,
    dpi = 600
)



# School Environment
school_vars <- exposome_variables_df %>%
    filter(category == "School Environment") %>%
    pull(variable_name)

school_raw <- discovery_df %>%
    select(any_of(school_vars)) %>%
    mutate_all(as.numeric) # 12
school_xcor <- school_raw %>%
    psych::mixedCor(ncat = 5) %>%
    pluck("rho")

school_par <- psych::fa.parallel(school_xcor, fa = "fa", n.obs = n_obs, n.iter = 1000)
for (ith in 4:6) {
    psych::fa.sort(psych::fa(school_raw, nfactors = ith, rotate = "promax")) %>%
        pluck("loadings") %>%
        print(cutoff = 0.2)
}

school_res <- psych::fa.sort(
    psych::fa(school_raw, nfactors = 5, rotate = "promax")
)
school_fa_labels <- c(
    "school engagement & environment", "early education access",
    "academic performance", "school resource availability", "education enrollment"
)
discovery_fascores_df <- discovery_fascores_df %>%
    bind_cols(
        school_res$scores %>%
            as_tibble() %>%
            set_names(school_fa_labels)
    )
holdout_fascores_df <- holdout_fascores_df %>%
    bind_cols(
        holdout_df %>%
            select(any_of(school_vars)) %>%
            mutate_all(as.numeric) %>%
            predict(school_res, data = .) %>%
            as_tibble() %>%
            set_names(school_fa_labels)
    )

school_fa_fig <- .plot_fa(school_par, school_res, school_fa_labels)
ggsave(
    plot = school_fa_fig,
    filename = here("outputs", "figures", "school_fa_fig.png"),
    height = 12,
    width = 14,
    dpi = 600
)


# Neighborhood Environment
negbenv_vars <- exposome_variables_df %>%
    filter(category == "Neighborhood Environment") %>%
    pull(variable_name)

negbenv_raw <- discovery_df %>%
    select(any_of(negbenv_vars)) %>%
    mutate_all(as.numeric) # 22
negbenv_xcor <- negbenv_raw %>%
    psych::mixedCor(ncat = 5) %>%
    pluck("rho")

negbenv_par <- psych::fa.parallel(negbenv_xcor, fa = "fa", n.obs = n_obs, n.iter = 1000)
for (ith in 7:9) {
    psych::fa.sort(
        psych::fa(negbenv_xcor, nfactors = ith, rotate = "promax")
    ) %>%
        pluck("loadings") %>%
        print(cutoff = 0.2)
}

negbenv_res <- psych::fa.sort(
    psych::fa(negbenv_raw, nfactors = 8, rotate = "promax")
)
negbenv_fa_labels <- c(
    "ethnic & racial diversity", "concentrated disadvantage",
    "Asian assimilation", "homogeneous Hispanic clustering", 
    "community racial/ethnic clustering", "homogeneous Asian clustering",
    "American Indian population", "Russian population"
)
discovery_fascores_df <- discovery_fascores_df %>%
    bind_cols(
        negbenv_res$scores %>%
            as_tibble() %>%
            set_names(negbenv_fa_labels)
    )
holdout_fascores_df <- holdout_fascores_df %>%
    bind_cols(
        holdout_df %>%
            select(any_of(negbenv_vars)) %>%
            mutate_all(as.numeric) %>%
            predict(negbenv_res, data = .) %>%
            as_tibble() %>%
            set_names(negbenv_fa_labels)
    )

negbenv_fa_fig <- .plot_fa(negbenv_par, negbenv_res, negbenv_fa_labels)
ggsave(
    plot = negbenv_fa_fig,
    filename = here("outputs", "figures", "negbenv_fa_fig.png"),
    height = 14,
    width = 16,
    dpi = 600
)


# Neighborhood Socioeconomic Status
negbses_vars <- exposome_variables_df %>%
    filter(category == "Neighborhood Socioeconomic Status") %>%
    pull(variable_name)

negbses_raw <- discovery_df %>%
    select(any_of(negbses_vars)) %>%
    mutate_all(as.numeric) # 34
negbses_xcor <- negbses_raw %>%
    psych::mixedCor(ncat = 5) %>%
    pluck("rho")

negbses_par <- psych::fa.parallel(negbses_xcor, fa = "fa", n.obs = n_obs, n.iter = 1000)
for (ith in 5:9) {
    psych::fa.sort(
        psych::fa(negbses_xcor, nfactors = ith, rotate = "promax")
    ) %>%
        pluck("loadings") %>%
        print(cutoff = 0.2)
}

negbses_res <- psych::fa.sort(
    psych::fa(negbses_raw, nfactors = 5, rotate = "promax")
)
negbses_fa_labels <- c(
    "neighborhood urbanicity", "neighbourhood wealth",
    "less education & crowding", "job density", "rent burden"
)
discovery_fascores_df <- discovery_fascores_df %>%
    bind_cols(
        negbses_res$scores %>%
            as_tibble() %>%
            set_names(negbses_fa_labels)
    )
holdout_fascores_df <- holdout_fascores_df %>%
    bind_cols(
        holdout_df %>%
            select(any_of(negbses_vars)) %>%
            mutate_all(as.numeric) %>%
            predict(negbses_res, data = .) %>%
            as_tibble() %>%
            set_names(negbses_fa_labels)
    )
negbses_fa_fig <- .plot_fa(negbses_par, negbses_res, negbses_fa_labels)
ggsave(
    plot = negbses_fa_fig,
    filename = here("outputs", "figures", "negbses_fa_fig.png"),
    height = 16,
    width = 16,
    dpi = 600
)


# Pollution
pollution_vars <- exposome_variables_df %>%
    filter(category == "Pollution") %>%
    pull(variable_name)

pollution_raw <- discovery_df %>%
    select(any_of(pollution_vars)) %>%
    mutate_all(as.numeric) # 29
pollution_xcor <- pollution_raw %>%
    psych::mixedCor(ncat = 5) %>%
    pluck("rho")

pollution_par <- psych::fa.parallel(pollution_xcor, fa = "fa", n.obs = n_obs, n.iter = 1000)
for (ith in 6:10) {
    psych::fa.sort(
        psych::fa(pollution_xcor, nfactors = ith, rotate = "promax")
    ) %>%
        pluck("loadings") %>%
        print(cutoff = 0.5)
}

pollution_res <- psych::fa.sort(
    psych::fa(pollution_raw, nfactors = 8, rotate = "promax")
)
pollution_fa_labels <- c(
    "environmental pollutants", "air quality", "soil fertility",
    "residential air pollutants", "heavy metal", "lead exposure risk",
    "vanadium exposure", "airborne toxins"
)
discovery_fascores_df <- discovery_fascores_df %>%
    bind_cols(
        pollution_res$scores %>%
            as_tibble() %>%
            set_names(pollution_fa_labels)
    )
holdout_fascores_df <- holdout_fascores_df %>%
    bind_cols(
        holdout_df %>%
            select(any_of(pollution_vars)) %>%
            mutate_all(as.numeric) %>%
            predict(pollution_res, data = .) %>%
            as_tibble() %>%
            set_names(pollution_fa_labels)
    )


pollution_fa_fig <- .plot_fa(pollution_par, pollution_res, pollution_fa_labels)
ggsave(
    plot = pollution_fa_fig,
    filename = here("outputs", "figures", "pollution_fa_fig.png"),
    height = 16,
    width = 16,
    dpi = 600
)

# Urbanization
urban_vars <- exposome_variables_df %>%
    filter(category == "Urbanization") %>%
    pull(variable_name)

urban_raw <- discovery_df %>%
    select(any_of(urban_vars)) %>%
    mutate_all(as.numeric) # 27
urban_xcor <- urban_raw %>%
    psych::mixedCor(ncat = 5) %>%
    pluck("rho")

urban_par <- psych::fa.parallel(urban_xcor, fa = "fa", n.obs = n_obs, n.iter = 1000)
for (ith in 9:13) {
    psych::fa.sort(psych::fa(urban_xcor, nfactors = ith, rotate = "promax")) %>%
        pluck("loadings") %>%
        print(cutoff = 0.2)
}

urban_res <- psych::fa.sort(
    psych::fa(urban_raw, nfactors = 9, rotate = "promax")
)
urban_fa_labels <- c(
    "green space", "urban density", "wood cover", "open space",
    "grass percentage", "agricultural land use", "urban development intensity",
    "vegetation cover", "presence of wetland"
)
discovery_fascores_df <- discovery_fascores_df %>%
    bind_cols(
        urban_res$scores %>%
            as_tibble() %>%
            set_names(urban_fa_labels)
    )
holdout_fascores_df <- holdout_fascores_df %>%
    bind_cols(
        holdout_df %>%
            select(any_of(urban_vars)) %>%
            mutate_all(as.numeric) %>%
            predict(urban_res, data = .) %>%
            as_tibble() %>%
            set_names(urban_fa_labels)
    )

urban_fa_fig <- .plot_fa(urban_par, urban_res, urban_fa_labels)
ggsave(
    plot = urban_fa_fig,
    filename = here("outputs", "figures", "urban_fa_fig.png"),
    height = 16,
    width = 16,
    dpi = 600
)

# Parental Psychopathology
parpsy_vars <- exposome_variables_df %>%
    filter(category == "Parental Psychopathology") %>%
    pull(variable_name)

parpsy_raw <- discovery_df %>%
    select(
        any_of(parpsy_vars), famhx_ss_fath_mh_history, famhx_ss_moth_mh_history
    ) %>%
    mutate_all(as.numeric) # 11
parpsy_xcor <- parpsy_raw %>%
    psych::mixedCor(ncat = 5) %>%
    pluck("rho")

parpsy_par <- psych::fa.parallel(parpsy_xcor, fa = "fa", n.obs = n_obs, n.iter = 1000)
for (ith in 2:6) {
    psych::fa.sort(
        psych::fa(parpsy_xcor, nfactors = ith, rotate = "promax")
    ) %>%
        pluck("loadings") %>%
        print(cutoff = 0.2)
}

parpsy_res <- psych::fa.sort(
    psych::fa(parpsy_raw, nfactors = 4, rotate = "promax")
)
parpsy_fa_labels <- c(
    "parental internalizing", "parental externalizing",
    "parental mental health history", "parental strength"
)
discovery_fascores_df <- discovery_fascores_df %>%
    bind_cols(
        parpsy_res$scores %>%
            as_tibble() %>%
            set_names(parpsy_fa_labels)
    )
holdout_fascores_df <- holdout_fascores_df %>%
    bind_cols(
        holdout_df %>%
            select(
                any_of(parpsy_vars), famhx_ss_fath_mh_history,
                famhx_ss_moth_mh_history
            ) %>%
            mutate_all(as.numeric) %>%
            predict(parpsy_res, data = .) %>%
            as_tibble() %>%
            set_names(parpsy_fa_labels)
    )

parpsy_fa_fig <- .plot_fa(parpsy_par, parpsy_res, parpsy_fa_labels)
ggsave(
    plot = parpsy_fa_fig,
    filename = here("outputs", "figures", "parpsy_fa_fig.png"),
    height = 14,
    width = 14,
    dpi = 600
)

# Neighborhood Safety
neigsafe_vars <- exposome_variables_df %>%
    filter(category == "Neighborhood Safety") %>%
    pull(variable_name)

neigsafe_raw <- discovery_df %>%
    select(any_of(neigsafe_vars)) %>%
    mutate_all(as.numeric) # 3
neigsafe_xcor <- neigsafe_raw %>%
    psych::mixedCor(ncat = 5) %>%
    pluck("rho")

neigsafe_par <- psych::fa.parallel(neigsafe_xcor, fa = "fa", n.obs = n_obs, n.iter = 1000)

neigsafe_res <- psych::fa.sort(
    psych::fa(neigsafe_raw, nfactors = 1, rotate = "promax")
)
discovery_fascores_df <- discovery_fascores_df %>%
    bind_cols(
        neigsafe_res$scores %>%
            as_tibble() %>%
            set_names("neighborhood safety")
    )
holdout_fascores_df <- holdout_fascores_df %>%
    bind_cols(
        holdout_df %>%
            select(any_of(neigsafe_vars)) %>%
            mutate_all(as.numeric) %>%
            predict(neigsafe_res, data = .) %>%
            as_tibble() %>%
            set_names("neighborhood safety")
    )

neigsafe_fa_fig <- .plot_fa(neigsafe_par, neigsafe_res, "neighborhood safety")
ggsave(
    plot = neigsafe_fa_fig,
    filename = here("outputs", "figures", "neigsafe_fa_fig.png"),
    height = 8,
    width = 6,
    dpi = 600
)

# Family Socioeconomic Status
ses_vars <- exposome_variables_df %>%
    filter(category == "Family Socioeconomic Status") %>%
    pull(variable_name)

ses_raw <- discovery_df %>%
    select(any_of(ses_vars)) %>%
    mutate_all(as.numeric) # 3
ses_xcor <- ses_raw %>%
    psych::mixedCor(ncat = 5) %>%
    pluck("rho")

ses_par <- psych::fa.parallel(ses_xcor, fa = "fa", n.obs = n_obs, n.iter = 1000)
for (ith in 1:2) {
    psych::fa.sort(psych::fa(ses_raw, nfactors = ith, rotate = "promax")) %>%
        pluck("loadings") %>%
        print(cutoff = 0.4)
}


ses_res <- psych::fa.sort(psych::fa(ses_raw, nfactors = 1, rotate = "promax"))
discovery_fascores_df <- discovery_fascores_df %>%
    bind_cols(
        ses_res$scores %>%
            as_tibble() %>%
            set_names("family SES")
    )
holdout_fascores_df <- holdout_fascores_df %>%
    bind_cols(
        holdout_df %>%
            select(any_of(ses_vars)) %>%
            mutate_all(as.numeric) %>%
            predict(ses_res, data = .) %>%
            as_tibble() %>%
            set_names("family SES")
    )

ses_fa_fig <- .plot_fa(ses_par, ses_res, "family SES")
ggsave(
    plot = ses_fa_fig,
    filename = here("outputs", "figures", "ses_fa_fig.png"),
    height = 8,
    width = 6,
    dpi = 600
)

# combine
master_fascores_df <- bind_rows(discovery_fascores_df, holdout_fascores_df)

write_rds(
    master_fascores_df,
    file = here("outputs", "caches", "master_fascores_df.rds")
)

