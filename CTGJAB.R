# Load necessary libraries
library(dplyr)
library(tidyr)
library(readr)
library(stringr)
library(purrr)
library(tibble)

# Load the data
study3 <- read_csv("../jtc/study3.csv")
source("../rfuncs/JAB.R")
# Clean the study3 data
clean_study3 <- function(data) {
  ids_to_remove <- c(
    "9f08f829a4d81c7b436b818dd972e1b4",
    "01bf66ca01ef01f21b9e3bfd18d4c81f",
    "25638018d5a7c1c3cd22fcddd0cf44d8",
    "b39984376c029bac02dc6de78af004cd"
  )
  data %>%
    mutate(analysisId = as.character(analysisId)) %>%                 # ensure character for matching
    filter(!analysisId %in% ids_to_remove) %>%                        # drop specified analyses
    mutate(has_less_than = grepl("<", pValue)) %>%
    mutate(pValue = str_extract(pValue, "[0-9]*\\.?[0-9]+")) %>%
    filter(!is.na(pValue) & pValue != "") %>%
    mutate(pValue = as.numeric(pValue)) %>%
    filter(pValue > 0 & pValue < 1 & value != 0)
}

study3_clean <- clean_study3(study3)

# Further clean and categorize the data
categorize_study3 <- function(data) {
  data %>%
    mutate(
      originalMethod = statisticalMethod,
    ) %>%
    mutate(
      statisticalMethod = case_when(
        # --- T-TESTS ---
        str_detect(originalMethod, regex("\\bt[- ]?test", ignore_case = TRUE)) &
          str_detect(originalMethod, regex("one[- ]?sample", ignore_case = TRUE)) ~ "One-sample t-test",

        str_detect(originalMethod, regex("\\bt[- ]?test", ignore_case = TRUE)) &
          str_detect(originalMethod, regex("paired", ignore_case = TRUE)) ~ "One-sample t-test",

        str_detect(originalMethod, regex("\\bt[- ]?test", ignore_case = TRUE)) ~ "Two-sample t-test",

        # --- CONDITIONAL LOGISTIC REGRESSION ---
        str_detect(originalMethod, regex("logistic", ignore_case = TRUE)) &
          str_detect(originalMethod, regex("regression", ignore_case = TRUE)) &
          str_detect(originalMethod, regex("conditional", ignore_case = TRUE)) ~ "Conditional logistic regression",

        # --- LOGISTIC REGRESSION ---
        str_detect(originalMethod, regex("logistic", ignore_case = TRUE)) &
          str_detect(originalMethod, regex("regression", ignore_case = TRUE)) ~ "Logistic regression",

        # --- COX ---
        str_detect(originalMethod, regex("cox", ignore_case = TRUE)) &
          !str_detect(originalMethod, regex("wilcoxon|mantel|signed|rank", ignore_case = TRUE)) ~ "Cox",

        # --- Logrank / Mantel-Cox ---
        str_detect(originalMethod, regex("log ?rank|mantel[- ]?cox", ignore_case = TRUE)) &
          !str_detect(originalMethod, regex("wilcoxon|signed", ignore_case = TRUE)) ~ "Logrank",

        # --- MANN-WHITNEY ---
        str_detect(originalMethod, regex("mann[- ]?whitney", ignore_case = TRUE)) |
          (str_detect(originalMethod, regex("wilcoxon", ignore_case = TRUE)) &
             str_detect(originalMethod, regex("rank", ignore_case = TRUE)) &
             !str_detect(originalMethod, regex("signed", ignore_case = TRUE))) ~ "Mann-Whitney",

        # --- WILCOXON ---
        str_detect(originalMethod, regex("wilcoxon", ignore_case = TRUE)) &
          str_detect(originalMethod, regex("signed", ignore_case = TRUE)) &
          str_detect(originalMethod, regex("rank", ignore_case = TRUE)) ~ "Wilcoxon",

        # --- KRUSKAL-WALLIS ---
        str_detect(originalMethod, regex("kruskal|wallis", ignore_case = TRUE)) ~ "Kruskal-Wallis",

        # --- CHI-SQUARE ---
        str_detect(originalMethod, regex("chi[- ]?square|χ2|chi2|chisq", ignore_case = TRUE)) ~ "Chi-square test",

        # --- REPEATED MEASURES ---
        str_detect(originalMethod, regex("repeated[- ]?measures?", ignore_case = TRUE)) ~ "Repeated measures analysis",

        # --- LINEAR REGRESSION ---
        str_detect(originalMethod, regex("linear", ignore_case = TRUE)) &
          str_detect(originalMethod, regex("regression", ignore_case = TRUE)) ~ "Linear regression",

        # --- ANOVA ---
        str_detect(originalMethod, regex("anova|analysis of variance", ignore_case = TRUE)) ~ "ANOVA",

        # --- DEFAULT: drop if no pattern matched ---
        TRUE ~ NA_character_
      )
    ) %>%
    # Drop observations where we failed to detect a test
    filter(!is.na(statisticalMethod)) %>%
    mutate(
      pValue = if_else(
        str_detect(originalMethod, regex("one[- ]?sided|1[- ]?sided", ignore_case = TRUE)),
        pmin(pValue * 2, 1),
        pValue
      )
    ) %>%
    group_by(analysisId) %>%
    filter(!(first(statisticalMethod) %in% c("Paired t-test", "One-sample t-test", "Wilcoxon") & n() != 1)) %>%
    ungroup()

}

study3_cleaned <- categorize_study3(study3_clean)

study3_summary <- study3_cleaned %>%
  group_by(analysisId) %>%
  summarise(
    # pValue and nctId are constant within each analysisId, so we take the first value.
    nctId = first(nctId),
    pValue = first(pValue),
    statisticalMethod = first(statisticalMethod),
    studyType = first(studyType),
    EARLY_PHASE1 = first(EARLY_PHASE1),
    PHASE1 = first(PHASE1),
    PHASE2 = first(PHASE2),
    PHASE3 = first(PHASE3),
    PHASE4 = first(PHASE4),
    outcomeTitle = first(outcomeTitle),
    outcomeType = first(outcomeType),
    groupDescription = first(groupDescription),
    conditionIds = first(condition_ids),
    briefTitle = first(briefTitle),
    has_less_than = first(has_less_than),
    # Sample size:
    # For one-sample t-tests and wilcoxon tests (i.e. Wilcoxon signed rank), take the first value,
    # For all other tests (e.g. two-sample t-tests, paired t-tests, etc.), sum the values.
    N = case_when(
      statisticalMethod %in% c("One-sample t-test", "Wilcoxon") ~ first(value),
      statisticalMethod == "Cox" ~ sum(value, na.rm = TRUE) * 0.5,
      statisticalMethod == "Logrank" ~ sum(value, na.rm = TRUE) * 0.5,
      statisticalMethod == "Conditional logistic regression" ~ max(value, na.rm = TRUE),
      TRUE ~ sum(value, na.rm = TRUE)
    ),

    # Map statisticalMethod to the expected model string for JAB01.
    model = case_when(
      statisticalMethod %in% c("Two-sample t-test", "Paired t-test", "One-sample t-test") ~ "t-test",
      statisticalMethod == "Linear regression" ~ "linear_regression",
      statisticalMethod %in% c("Logistic regression","Conditional logistic regression") ~ "logistic_regression",
      #statisticalMethod == "Cox" ~ "cox",
      #statisticalMethod == "Logrank" ~ "cox",
      statisticalMethod == "ANOVA" ~ "anova",
      statisticalMethod == "Kruskal-Wallis" ~ "kruskal_wallis",   # kruskal wallis follows ANOVA logic for parameters.
      statisticalMethod == "Mann-Whitney" ~ "mann_whitney",       # mann-whitney follows two-sample t-test logic.
      statisticalMethod == "Wilcoxon" ~ "wilcoxon",      # wilcoxon (signed rank) uses first value only.
      #statisticalMethod == "Repeated measures analysis" ~ "repeated_measures",
      #statisticalMethod == "Chi-square test" ~ "chi_squared",
      TRUE ~ NA_character_
    ),

    # For ANOVA, kruskal wallis, and repeated measures, use the number of rows in the group.
    I = if_else(statisticalMethod %in% c("ANOVA", "Repeated measures analysis", "Kruskal-Wallis"),
                n(),
                NA_integer_),

    # For chi-square tests: R is the number of outcome groups (rows) and C is fixed at 2.
    R = if_else(statisticalMethod == "Chi-square test", n(), NA_integer_),
    C = if_else(statisticalMethod == "Chi-square test", 2L, NA_integer_)
  ) %>%
  # Drop rows where the model mapping did not apply.
  filter(!is.na(model)) %>%
  rowwise() %>%  # Ensure rowwise operation for computing JAB.
  mutate(
    JAB = tryCatch(
      JAB01(N, pValue, model, R = R, C = C, I = I),
      error = function(e) NA_real_
    )
  ) %>%
  ungroup() %>% filter(!is.na(JAB), !is.na(pValue), N != 1)
