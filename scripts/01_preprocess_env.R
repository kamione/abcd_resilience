# Environment ------------------------------------------------------------------
library(tidyverse)
library(here)
library(glue)
library(tidymodels)
library(missForestPredict)
library(doParallel)

source(here("src", "R", "sample_with_seed.R"))
source(here("src", "R", "is_binary.R"))


# Data I/O ---------------------------------------------------------------------
# using data release 5.1
# please change the path to ABCD study folder
abcd_folder_path <- "path_to_abcd_folder" %>%
    here("abcd-data-release-5.1/core")

# general
abcd_y_lt <- here(abcd_folder_path, "abcd-general", "abcd_y_lt.csv") %>%
    read_csv(show_col_types = FALSE)
abcd_p_demo <- here(abcd_folder_path, "abcd-general", "abcd_p_demo.csv") %>%
    read_csv(show_col_types = FALSE) %>%
    select(-c(demo_comb_income_v2, demo_prnt_prtnr_v2, demo_prnt_ed_v2))

# exposome variables
exposome_variables_df <- here("data", "raw", "included_variables.xlsx") %>%
    readxl::read_excel(sheet = "Exposome Raw")
selected_filenames <- c(unique(pull(exposome_variables_df, table_name)))
selected_variables <- c(pull(exposome_variables_df, variable_name)) # 344

selected_files <- list()
for (ith in seq_along(selected_filenames)) {
    selected_files[[ith]] <- list.files(
        path = abcd_folder_path,
        full.names = TRUE,
        pattern = glue("{selected_filenames[ith]}.csv"),
        recursive = TRUE
    ) %>%
        read_csv(show_col_types = FALSE) %>%
        select(src_subject_id, eventname, any_of(selected_variables))
}

# combine
master_df <- selected_files %>%
    reduce(left_join, by = c("src_subject_id", "eventname")) %>%
    right_join(abcd_p_demo, by = c("src_subject_id", "eventname")) %>%
    right_join(abcd_y_lt, by = c("src_subject_id", "eventname"))


# Preprocessing ----------------------------------------------------------------
master_preprocessed_df <- master_df %>%
    filter(eventname == "baseline_year_1_arm_1") %>% # 11,868 subjects
    group_by(rel_family_id) %>%
    mutate(random_order_id = sample_with_seed(rel_family_id, n())) %>%
    filter(random_order_id == 1) %>% # 9,807 subjects
    ungroup() %>%
    mutate(
        across(where(is.numeric), ~na_if(., 777)),
        across(where(is.numeric), ~na_if(., 999)),
        # these variables should be numeric
        devhx_15 = as.numeric(devhx_15),
        devhx_16_p = as.numeric(devhx_16_p),
        devhx_16_p = if_else(devhx_16_p == 9990, NA, devhx_16_p),
        devhx_17_p = as.numeric(devhx_17_p)
    ) %>%
    filter(
        demo_sex_v2 %in% c(1, 2), # 9,804 subjects
    ) %>%
    mutate(
        site_id_l = str_remove(site_id_l, "site"),
        site_id_l = factor(as.numeric(site_id_l)),
        demo_sex_v2 = factor(demo_sex_v2, labels = c("Male", "Female")),
        race_white = if_else(demo_race_a_p___10 == 1, 1, 0),
        race_black = if_else(demo_race_a_p___11 == 1, 1, 0),
        race_asian = if_else(
            (demo_race_a_p___18 + demo_race_a_p___19 + demo_race_a_p___20 +
                 demo_race_a_p___21 + demo_race_a_p___22 + demo_race_a_p___23 +
                 demo_race_a_p___24) > 0, 1, 0
        ),
        race_aian = if_else(
            (demo_race_a_p___12 + demo_race_a_p___13) > 0, 1, 0
        ),
        race_nhpi = if_else(
            (demo_race_a_p___14 + demo_race_a_p___15 + demo_race_a_p___17 +
                 demo_race_a_p___17) > 0, 1, 0
        ),
        race_other = if_else(demo_race_a_p___25 == 1, 1, 0),
        race_mixed = if_else(
            (race_white + race_black + race_asian + race_aian + race_nhpi +
                 race_other) > 1, 1, 0
        ),
        race_ethnicity_3 = case_when(
            race_white == 1 & race_mixed == 0 & demo_ethn_v2 == 2 ~ 1,
            race_black == 1 & race_mixed == 0 ~ 2,
            demo_race_a_p___77 | demo_race_a_p___99 ~ NA,
            TRUE ~ 3
        ),
        race_ethnicity_3 = factor(
            race_ethnicity_3,
            levels = 1:3,
            label = c("White", "Black", "Other")
        ),
        demo_comb_income_v2 = case_when(
            demo_comb_income_v2 <= 4 ~ 1,
            demo_comb_income_v2 > 4 & demo_comb_income_v2 <= 6 ~ 2,
            demo_comb_income_v2 == 7 ~ 3,
            demo_comb_income_v2 == 8 ~ 4,
            demo_comb_income_v2 == 9 ~ 5,
            demo_comb_income_v2 == 10 ~ 6,
            TRUE ~ NA
        ),
        demo_comb_income_v2 = ordered(demo_comb_income_v2)
    ) %>%
    drop_na(
        src_subject_id, site_id_l, interview_age, demo_sex_v2, race_ethnicity_3
    ) %>% # 9,664 subjects
    select(
        src_subject_id, site_id_l, interview_age, demo_sex_v2, race_ethnicity_3,
        all_of(selected_variables)
    ) %>%
    mutate_if(is_binary, factor) # 9,664 x 348

# pull participants with more than 15% missing variables
remove_id <- master_preprocessed_df %>%
    is.na() %>%
    rowSums() %>%
    as_tibble() %>%
    mutate(
        src_subject_id = master_preprocessed_df$src_subject_id,
        prop_missing = value / dim(master_preprocessed_df)[2]
    ) %>%
    filter(prop_missing > 0.15) %>%
    pull(src_subject_id) # 580 subjects

# pull variables with more than 15% missing values
high_incomplete_variables <- master_preprocessed_df %>%
    select(6:348) %>%
    mutate_all(as.numeric) %>%
    summarise_all(~ mean(is.na(.))) %>%
    pivot_longer(
        cols = everything(),
        names_to = "variable",
        values_to = "prop_missing"
    ) %>%
    filter(prop_missing > 0.15) %>%
    pull(variable) # 12 variables

master_preprocessed_cleaned_df <- master_preprocessed_df %>%
    filter(!(src_subject_id %in% remove_id)) %>%
    select(
        -all_of(high_incomplete_variables)
    ) # 9,084 x 336 (env: 331 variables)

set.seed(1234)
df_split <- initial_split(
    master_preprocessed_cleaned_df,
    prop = 2 / 3
)
discovery_ids <- training(df_split) %>% pull(src_subject_id) # 6,056
holdout_ids <- testing(df_split) %>% pull(src_subject_id) # 3,028

master_preprocessed_cleaned_df <- master_preprocessed_cleaned_df %>%
    mutate(
        dataset = case_when(
            src_subject_id %in% discovery_ids ~ "Discovery",
            src_subject_id %in% holdout_ids ~ "Holdout"
        )
    )


# Missing Data Preprocessing ---------------------------------------------------
# this will take ~12 minutes per iteration with 6 cores
set.seed(1234)
n_cores <- detectCores() - 2 # avoid exhaustion of CPU cores
cl <- makeCluster(n_cores)
registerDoParallel(cl)
missforest_imputation <- master_preprocessed_cleaned_df %>%
    filter(dataset == "Discovery") %>%
    # remove character for imputation
    select_if(negate(is.character)) %>%
    as.data.frame() %>%
    missForestPredict::missForest(
        verbose = TRUE,
        save_models = TRUE,
        maxiter = 5,
        num.threads = n_cores
    )
stopCluster(cl)

write_rds(
    missforest_imputation,
    compress = c("gz"),
    file = here("outputs", "caches", "missforest_imputation.rds")
)

holdout_imputated_df <- missForestPredict::missForestPredict(
    missForestObj = missforest_imputation,
    newdata = filter(master_preprocessed_cleaned_df, dataset == "Holdout")
)

discovery_imputated_df <- master_preprocessed_cleaned_df %>%
    filter(dataset == "Discovery") %>%
    select_if(is.character) %>%
    bind_cols(missforest_imputation$ximp)

master_imputated_df <- bind_rows(
    discovery_imputated_df, holdout_imputated_df
)

write_rds(
    master_imputated_df,
    file = here("outputs", "caches", "master_imputated_df.rds")
)



# Variable Mutation ------------------------------------------------------------
master_processed_imputated_df <- master_imputated_df %>%
    mutate_at(vars(contains("famhx_ss_fath_prob_")), as.numeric) %>%
    mutate_at(vars(contains("famhx_ss_moth_prob_")), as.numeric) %>%
    mutate_at(vars(contains("devhx_10")), as.numeric) %>%
    mutate_at(vars(contains("devhx_14")), as.numeric) %>%
    mutate_at(vars(contains("ksads_ptsd_raw")), as.numeric) %>%
    # the following calculations are performed on a row-wise basis
    mutate(
        birth_weight_kg = birth_weight_lbs * 0.45359237 +
            birth_weight_oz * 0.0283495231,
        famhx_ss_fath_mh_history = famhx_ss_fath_prob_alc_p +
            famhx_ss_fath_prob_dg_p + famhx_ss_fath_prob_dprs_p +
            famhx_ss_fath_prob_ma_p + famhx_ss_fath_prob_vs_p +
            famhx_ss_fath_prob_trb_p + famhx_ss_fath_prob_nrv_p +
            famhx_ss_fath_prob_prf_p + famhx_ss_fath_prob_hspd_p +
            famhx_ss_fath_prob_scd_p,
        famhx_ss_moth_mh_history = famhx_ss_moth_prob_alc_p +
            famhx_ss_moth_prob_dg_p + famhx_ss_moth_prob_dprs_p +
            famhx_ss_moth_prob_ma_p + famhx_ss_moth_prob_vs_p +
            famhx_ss_moth_prob_trb_p + famhx_ss_moth_prob_nrv_p +
            famhx_ss_moth_prob_prf_p + famhx_ss_moth_prob_hspd_p +
            famhx_ss_moth_prob_scd_p,
        devhx_10sum3 = devhx_10a3_p + devhx_10b3_p + devhx_10c3_p +
            devhx_10d3_p + devhx_10e3_p + devhx_10f3_p + devhx_10g3_p +
            devhx_10h3_p + devhx_10i3_p + devhx_10j3_p + devhx_10k3_p +
            devhx_10l3_p + devhx_10m3_p,
        devhx_14sum3 = devhx_14a3_p + devhx_14b3_p + devhx_14c3_p +
            devhx_14d3_p + devhx_14e3_p + devhx_14f3_p + devhx_14g3_p +
            devhx_14h3_p,
        devhx_ss_cigs_per_day_p = (
            devhx_ss_8_cigs_per_day_p + devhx_ss_9_cigs_per_day_p
        ) / 2,
        devhx_ss_alcohol_avg_p = (
            devhx_ss_8_alcohol_avg_p + devhx_ss_9_alcohol_avg_p
        ) / 2,
        devhx_12 = ifelse(devhx_12a_p == 0, 0, devhx_ss_12_p),
        ksads_ptsd_sum = ksads_ptsd_raw_754_p + ksads_ptsd_raw_755_p +
            ksads_ptsd_raw_756_p + ksads_ptsd_raw_757_p + ksads_ptsd_raw_758_p +
            ksads_ptsd_raw_759_p + ksads_ptsd_raw_760_p + ksads_ptsd_raw_761_p +
            ksads_ptsd_raw_762_p + ksads_ptsd_raw_763_p + ksads_ptsd_raw_764_p +
            ksads_ptsd_raw_765_p + ksads_ptsd_raw_766_p + ksads_ptsd_raw_767_p +
            ksads_ptsd_raw_768_p + ksads_ptsd_raw_769_p + ksads_ptsd_raw_770_p -
            17
    ) %>%
    select(-c(
        birth_weight_lbs,
        birth_weight_oz,
        devhx_ss_8_cigs_per_day_p,
        devhx_ss_9_cigs_per_day_p,
        devhx_ss_8_alcohol_avg_p,
        devhx_ss_9_alcohol_avg_p,
        devhx_12a_p,
        devhx_ss_12_p,
        matches("devhx_19[abcd]_p"),
        matches("famhx_ss_fath_prob_"),
        matches("famhx_ss_moth_prob_"),
        matches("devhx_10[abcdefghijklm]3"),
        matches("devhx_14[abcdefgh]3"),
        matches("ksads_ptsd_raw_7")
    )) # 9,084 × 279

# pull variables with high correlation (r >= 0.9) in the discovery dataset
high_corr_vars <- master_processed_imputated_df %>%
    filter(dataset == "Discovery") %>% 
    select_if(negate(is.character)) %>%
    correlation::correlation() %>%
    as_tibble() %>%
    filter(r >= 0.9) %>%
    pluck("Parameter2") %>%
    unique() # 37 variables

# pull variables with near zero variance 
nzv_vars <- master_processed_imputated_df %>%
    caret::nearZeroVar(names = TRUE) # 18 variables

master_cleaned_imputated_df <- master_processed_imputated_df %>% 
    select(
        -any_of(high_corr_vars),
        -any_of(nzv_vars)
    ) # 9,084 × 224

write_rds(
    master_cleaned_imputated_df,
    file = here("data", "processed", "master_cleaned_imputated_df.rds")
)
