# ==========================================================================
# Project: AI PRAME Study
# Script:  10_Supplementary_Four_Field_Analysis
# Author:  Frank Friedrich Gellrich
# Date:    2026-04-11
# Description: Human-AI Agreement Analysis (Confusion Matrix).
#              - Specific denominator for Ki-67 (only evaluated cases).
#              - Comprehensive legend and consistent formatting.
#              Outputs: Formatted Word document (landscape) and CSV data.
# ==========================================================================

# Part 1 - Packages ----
needed_packages <- c("tidyverse", "here", "flextable", "officer", "broom")
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
  "Intratumorales_lymphozytaeres_Infiltrat_TILs" = "Tumor-Infiltrating Lymphocytes", 
  "Infiltration_der_Subkutis" = "Subcutis Infiltration"
)

# Part 3 - Calculation Function ----
get_agreement_table <- function(train_data, test_data, target_set = "train") {
  data <- if(target_set == "train") train_data else test_data
  n_total <- nrow(data)
  results <- list()
  
  # A. Morphological Criteria (Calibrated Prediction)
  for(orig in names(morph_map)) {
    eng <- morph_map[orig]
    gt_col <- paste0(orig, "_GroundTruth"); score_col <- paste0(orig, "_ja")
    
    train_temp <- train_data %>%
      mutate(y = ifelse(!!sym(gt_col) %in% c("ja", 1), 1, 0)) %>%
      select(y, score = !!sym(score_col)) %>% drop_na()
    fit <- glm(y ~ score, data = train_temp, family = binomial)
    
    data_temp <- data %>%
      mutate(gt_pos = !!sym(gt_col) %in% c("ja", 1)) %>%
      select(gt_pos, score = !!sym(score_col))
    data_temp$ai_pos <- predict(fit, newdata = data_temp, type = "response") > 0.5
    
    tp <- sum(data_temp$gt_pos & data_temp$ai_pos, na.rm = TRUE)
    fn <- sum(data_temp$gt_pos & !data_temp$ai_pos, na.rm = TRUE)
    fp <- sum(!data_temp$gt_pos & data_temp$ai_pos, na.rm = TRUE)
    tn <- sum(!data_temp$gt_pos & !data_temp$ai_pos, na.rm = TRUE)
    
    results[[eng]] <- tibble(
      Criterion = eng,
      col_tp = sprintf("%d (%.1f%%)", tp, (tp/n_total)*100),
      col_fn = sprintf("%d (%.1f%%)", fn, (fn/n_total)*100),
      col_fp = sprintf("%d (%.1f%%)", fp, (fp/n_total)*100),
      col_tn = sprintf("%d (%.1f%%)", tn, (tn/n_total)*100),
      `Total Human GT` = sprintf("%d (%.1f%%)", (tp+fn), ((tp+fn)/n_total)*100),
      `Total AI Pred`  = sprintf("%d (%.1f%%)", (tp+fp), ((tp+fp)/n_total)*100)
    )
  }
  
  # B. PRAME Status
  for(v in 0:4) {
    gt_vec <- data$PRAME_Status_GroundTruth == v
    ai_vec <- data$PRAME_Status_pred == v
    tp <- sum(gt_vec & ai_vec, na.rm = TRUE); fn <- sum(gt_vec & !ai_vec, na.rm = TRUE)
    fp <- sum(!gt_vec & ai_vec, na.rm = TRUE); tn <- sum(!gt_vec & !ai_vec, na.rm = TRUE)
    name <- paste0("PRAME Status: ", v)
    results[[name]] <- tibble(
      Criterion = name, col_tp = sprintf("%d (%.1f%%)", tp, (tp/n_total)*100),
      col_fn = sprintf("%d (%.1f%%)", fn, (fn/n_total)*100), col_fp = sprintf("%d (%.1f%%)", fp, (fp/n_total)*100),
      col_tn = sprintf("%d (%.1f%%)", tn, (tn/n_total)*100),
      `Total Human GT` = sprintf("%d (%.1f%%)", (tp+fn), ((tp+fn)/n_total)*100),
      `Total AI Pred`  = sprintf("%d (%.1f%%)", (tp+fp), ((tp+fp)/n_total)*100)
    )
  }
  
  # C. Ki-67 (Specific denominator based on available data)
  mib_map <- c("Niedrig" = "Low", "Intermediaer" = "Intermediate", "Hoch" = "High")
  n_mib_available <- sum(!is.na(data$Mib1_GroundTruth))
  
  for(orig_m in names(mib_map)) {
    eng_m <- mib_map[orig_m]
    gt_vec <- str_detect(as.character(data$Mib1_GroundTruth), orig_m)
    ai_vec <- str_detect(as.character(data$Mib1_pred), orig_m)
    tp <- sum(gt_vec & ai_vec, na.rm = TRUE); fn <- sum(gt_vec & !ai_vec, na.rm = TRUE)
    fp <- sum(!gt_vec & ai_vec, na.rm = TRUE); tn <- sum(!gt_vec & !ai_vec, na.rm = TRUE)
    
    name <- paste0("Ki-67 Index: ", eng_m)
    results[[name]] <- tibble(
      Criterion = name,
      col_tp = sprintf("%d (%.1f%%)", tp, (tp/n_mib_available)*100),
      col_fn = sprintf("%d (%.1f%%)", fn, (fn/n_mib_available)*100),
      col_fp = sprintf("%d (%.1f%%)", fp, (fp/n_mib_available)*100),
      col_tn = sprintf("%d (%.1f%%)", tn, (tn/n_mib_available)*100),
      `Total Human GT` = sprintf("%d (%.1f%%)", (tp+fn), ((tp+fn)/n_mib_available)*100),
      `Total AI Pred`  = sprintf("%d (%.1f%%)", (tp+fp), ((tp+fp)/n_mib_available)*100)
    )
  }
  return(bind_rows(results))
}

# Part 4 - Flextable Formatting ----
ft_formatted <- function(df, title) {
  legend_text <- "Note: Confusion matrix of human-AI agreement. TP (True Positive) and TN (True Negative) indicate consensus; FN (False Negative) represents features identified by humans only; FP (False Positive) represents features identified by AI only. Marginal totals (Total Human GT / Total AI Pred) provide the overall frequency of feature detection. AI classifications for morphological criteria were optimized via probability-calibrated logistic regression models. Percentages for Ki-67 reflect only evaluable cases."
  
  flextable(df) %>%
    set_caption(title) %>% theme_vanilla() %>% autofit() %>%
    bold(part = "header") %>% bold(j = 1, part = "body") %>%
    align(j = 2:7, align = "center", part = "all") %>%
    set_header_labels(col_tp = "Agreement Both (TP)", col_fn = "Human only (FN)", col_fp = "AI only (FP)", col_tn = "Agreement Both (TN)",
                      `Total Human GT` = "Total Human GT (+)", `Total AI Pred` = "Total AI Pred (+)") %>%
    add_footer_lines(legend_text) %>%
    bg(j = 6:7, bg = "#F2F2F2", part = "all")
}

tab_cv80     <- get_agreement_table(df_train, df_test, "train")
tab_deploy20 <- get_agreement_table(df_train, df_test, "test")

# Part 5 - Word Export (Landscape) ----
doc <- read_docx() %>%
  body_add_par("Supplementary Material: Human-AI Agreement Analysis", style = "heading 1") %>%
  body_add_flextable(ft_formatted(tab_cv80, "Supplementary Table S3: Agreement Matrix (Training Set, N=800)")) %>% body_add_break() %>%
  body_add_flextable(ft_formatted(tab_deploy20, "Supplementary Table S4: Agreement Matrix (Independent Test Set, N=200)")) %>%
  body_end_section_landscape()

print(doc, target = file.path(table_dir, "10_Supplementary_Four_Field_Agreement.docx"))

# Part 6 - CSV Export (Data for downstream Python visualizations) ----
dir.create(csv_dir, recursive = TRUE, showWarnings = FALSE)

write_csv(tab_cv80, file.path(csv_dir, "10_Agreement_Analysis_Training_Set.csv"))
write_csv(tab_deploy20, file.path(csv_dir, "10_Agreement_Analysis_Test_Set.csv"))

cat("\n--- Script 10 (Supplementary Four Field Analysis) successfully executed. ---\n")