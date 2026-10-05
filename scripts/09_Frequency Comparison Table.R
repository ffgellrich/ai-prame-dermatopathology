# ==========================================================================
# Project: AI PRAME Study
# Script:  09_Supplementary_Frequency_Deep_Dive
# Author:  Frank Friedrich Gellrich
# Date:    2026-04-10
# Description: Supplementary tables with calibrated AI feature frequencies.
#              Binary AI feature calls are recalibrated via logistic regression 
#              (GT ~ AI_Score) on the training set.
#              Outputs: Formatted Word document (landscape) and CSV data.
# ==========================================================================

# Part 1 - Packages ----
needed_packages <- c("tidyverse", "here", "flextable", "officer")
invisible(lapply(needed_packages, library, character.only = TRUE))

# Part 2 - Configuration & Data ----
path_cv80     <- here("data", "processed", "Master_Table_AI_Features_CV80.csv")
path_deploy20 <- here("data", "processed", "Master_Table_AI_Features_Deploy20.csv")
table_dir     <- here("results", "tables")
csv_dir       <- here("results", "csv")
dir.create(table_dir, showWarnings = FALSE, recursive = TRUE)

df_train <- read_csv(path_cv80, show_col_types = FALSE)
df_test  <- read_csv(path_deploy20, show_col_types = FALSE)

# English feature name mapping
morph_map <- c(
  "KI_maligne" = "Overall Malignancy",
  "Symmetrie" = "Symmetry", "Aufsteigende_Melanozyten" = "Upward Melanocytes",
  "Pagetoide_Durchsetzung" = "Pagetoid Spread", "Konfluenz_junktionaler_Nester" = "Nest Confluence",
  "Ausgepraegte_junktionale_Atypien" = "Junctional Atypia", "Atypische_junktionale_Mitosen" = "Junctional Mitoses",
  "Fehlende_Ausreifung_in_Nestern" = "Lack of Maturation", "Ausgepraegte_dermale_Atypien" = "Dermal Atypia",
  "Dermale_Mitosen" = "Dermal Mitoses", "Geordnete_Vereinzelung_Benignes_Muster" = "Ordered Single Cells",
  "Atypische_Vereinzelung_Dissoziation_Malignes_Muster" = "Atypical Dissociation",
  "Regression" = "Regression", "Ulzeration" = "Ulceration",
  "Junktionale_Beteiligung" = "Junctional Component", "Dermale_Beteiligung" = "Dermal Component",
  "Dichtes_lymphozytaeres_Infiltrat" = "Dense Lymphocytic Infiltrate",
  "Intratumorales_lymphozytaeres_Infiltrat_TILs" = "Tumor-Infiltrating Lymphocytes (TILs)",
  "Infiltration_der_Subkutis" = "Subcutis Infiltration"
)

# Part 3 - Calibration & Table Generation Function ----

get_calibrated_table <- function(train_data, test_data, target_set = "train") {
  
  # Select target dataset
  data <- if(target_set == "train") train_data else test_data
  
  # Fixed denominators based on Ground Truth malignancy status
  n_mal <- sum(data$KI_maligne_GroundTruth %in% c("ja", 1), na.rm = TRUE)
  n_ben <- sum(!(data$KI_maligne_GroundTruth %in% c("ja", 1)), na.rm = TRUE)
  
  results <- list()
  
  # A. Morphology & Malignancy (Calibrated binary prediction)
  for(orig in names(morph_map)) {
    eng <- morph_map[orig]
    gt_col <- paste0(orig, "_GroundTruth")
    score_col <- paste0(orig, "_ja") # Continuous prediction score
    
    # 1. Fit logistic regression calibration model on training set
    train_temp <- train_data %>%
      mutate(y = ifelse(!!sym(gt_col) %in% c("ja", 1), 1, 0)) %>%
      select(y, score = !!sym(score_col)) %>% drop_na()
    
    fit <- glm(y ~ score, data = train_temp, family = binomial)
    
    # 2. Predict probabilities on target dataset (training or test)
    data_temp <- data %>% 
      mutate(gt_mal = ifelse(KI_maligne_GroundTruth %in% c("ja", 1), "mal", "ben")) %>%
      select(gt_mal, gt_feat = !!sym(gt_col), score = !!sym(score_col))
    
    # Calibrated binary classification (threshold > 0.5)
    data_temp$ai_ja <- predict(fit, newdata = data_temp, type = "response") > 0.5
    
    # Count frequencies
    c_gt_mal <- sum(data_temp$gt_feat %in% c("ja", 1) & data_temp$gt_mal == "mal", na.rm = TRUE)
    c_gt_ben <- sum(data_temp$gt_feat %in% c("ja", 1) & data_temp$gt_mal == "ben", na.rm = TRUE)
    c_ai_mal <- sum(data_temp$ai_ja == TRUE & data_temp$gt_mal == "mal", na.rm = TRUE)
    c_ai_ben <- sum(data_temp$ai_ja == TRUE & data_temp$gt_mal == "ben", na.rm = TRUE)
    
    results[[eng]] <- tibble(
      Criterion = eng, Attribute = "Present",
      `GT Malignant` = sprintf("%d/%d (%.1f%%)", c_gt_mal, n_mal, (c_gt_mal/n_mal)*100),
      `GT Benign`    = sprintf("%d/%d (%.1f%%)", c_gt_ben, n_ben, (c_gt_ben/n_ben)*100),
      `AI Malignant` = sprintf("%d/%d (%.1f%%)", c_ai_mal, n_mal, (c_ai_mal/n_mal)*100),
      `AI Benign`    = sprintf("%d/%d (%.1f%%)", c_ai_ben, n_ben, (c_ai_ben/n_ben)*100)
    )
  }
  
  # B. PRAME Status (Using existing predicted classes 0 to 4)
  for(v in 0:4) {
    gt_c <- "PRAME_Status_GroundTruth"
    ai_c <- "PRAME_Status_pred"
    c_gt_mal <- sum(data[[gt_c]] == v & data$KI_maligne_GroundTruth %in% c("ja", 1), na.rm = TRUE)
    c_gt_ben <- sum(data[[gt_c]] == v & !(data$KI_maligne_GroundTruth %in% c("ja", 1)), na.rm = TRUE)
    c_ai_mal <- sum(data[[ai_c]] == v & data$KI_maligne_GroundTruth %in% c("ja", 1), na.rm = TRUE)
    c_ai_ben <- sum(data[[ai_c]] == v & !(data$KI_maligne_GroundTruth %in% c("ja", 1)), na.rm = TRUE)
    
    results[[paste0("PRAME_", v)]] <- tibble(
      Criterion = "PRAME Status", Attribute = paste0("Score ", v),
      `GT Malignant` = sprintf("%d/%d (%.1f%%)", c_gt_mal, n_mal, (c_gt_mal/n_mal)*100),
      `GT Benign`    = sprintf("%d/%d (%.1f%%)", c_gt_ben, n_ben, (c_gt_ben/n_ben)*100),
      `AI Malignant` = sprintf("%d/%d (%.1f%%)", c_ai_mal, n_mal, (c_ai_mal/n_mal)*100),
      `AI Benign`    = sprintf("%d/%d (%.1f%%)", c_ai_ben, n_ben, (c_ai_ben/n_ben)*100)
    )
  }
  
  # C. Ki-67 (Proliferation index levels)
  mib_map <- c("Niedrig" = "Low (<5%)", "Intermediaer" = "Intermediate (5-20%)", "Hoch" = "High (>20%)")
  for(orig_m in names(mib_map)) {
    gt_c <- "Mib1_GroundTruth"
    ai_c <- "Mib1_pred"
    n_mal_mib <- sum(!is.na(data[[gt_c]]) & data$KI_maligne_GroundTruth %in% c("ja", 1))
    n_ben_mib <- sum(!is.na(data[[gt_c]]) & !(data$KI_maligne_GroundTruth %in% c("ja", 1)))
    
    c_gt_mal <- sum(str_detect(data[[gt_c]], orig_m) & data$KI_maligne_GroundTruth %in% c("ja", 1), na.rm = TRUE)
    c_gt_ben <- sum(str_detect(data[[gt_c]], orig_m) & !(data$KI_maligne_GroundTruth %in% c("ja", 1)), na.rm = TRUE)
    c_ai_mal <- sum(str_detect(data[[ai_c]], orig_m) & data$KI_maligne_GroundTruth %in% c("ja", 1), na.rm = TRUE)
    c_ai_ben <- sum(str_detect(data[[ai_c]], orig_m) & !(data$KI_maligne_GroundTruth %in% c("ja", 1)), na.rm = TRUE)
    
    results[[paste0("Mib1_", orig_m)]] <- tibble(
      Criterion = "Ki-67 Index (Mib1)", Attribute = mib_map[orig_m],
      `GT Malignant` = sprintf("%d/%d (%.1f%%)", c_gt_mal, n_mal_mib, (c_gt_mal/n_mal_mib)*100),
      `GT Benign`    = sprintf("%d/%d (%.1f%%)", c_gt_ben, n_ben_mib, (c_gt_ben/n_ben_mib)*100),
      `AI Malignant` = sprintf("%d/%d (%.1f%%)", c_ai_mal, n_mal_mib, (c_ai_mal/n_mal_mib)*100),
      `AI Benign`    = sprintf("%d/%d (%.1f%%)", c_ai_ben, n_ben_mib, (c_ai_ben/n_ben_mib)*100)
    )
  }
  return(bind_rows(results))
}

# Part 4 - Execution & Flextable Formatting ----

tab_cv80     <- get_calibrated_table(df_train, df_test, "train")
tab_deploy20 <- get_calibrated_table(df_train, df_test, "test")

ft_formatted <- function(df, title) {
  flextable(df) %>%
    set_caption(title) %>% theme_vanilla() %>% merge_v(j = 1) %>% fix_border_issues() %>% autofit() %>%
    bold(part = "header") %>% align(j = 3:6, align = "center", part = "all") %>%
    add_footer_lines("Note: AI binary predictions for morphological criteria were recalculated using logistic regression models calibrated on the training set (GT_Feature ~ AI_Feature_Score) to provide optimized frequency estimates beyond the standard 0.5 threshold. Malignant/Benign stratification is based on the Ground Truth diagnosis.")
}

ft_cv80     <- ft_formatted(tab_cv80, "Supp. Table S1: Calibrated Feature Frequencies (Training Set, CV80)")
ft_deploy20 <- ft_formatted(tab_deploy20, "Supp. Table S2: Calibrated Feature Frequencies (Independent Test Set, Deploy20)")

# Part 5 - Word Export (Landscape) ----
doc <- read_docx() %>%
  body_add_par("Supplementary Material: Optimized Feature Frequency Analysis", style = "heading 1") %>%
  body_add_flextable(ft_cv80) %>% body_add_break() %>%
  body_add_flextable(ft_deploy20) %>%
  body_end_section_landscape()

print(doc, target = file.path(table_dir, "09_Supplementary_Calibrated_Frequencies.docx"))

# Part 6 - CSV Export (Data for downstream Python visualizations) ----
dir.create(csv_dir, recursive = TRUE, showWarnings = FALSE)

write_csv(tab_cv80, file.path(csv_dir, "09_Calibrated_Frequencies_Training_Set.csv"))
write_csv(tab_deploy20, file.path(csv_dir, "09_Calibrated_Frequencies_Test_Set.csv"))

cat("\n--- Script 09 (Calibrated Frequency Table) successfully executed. ---\n")