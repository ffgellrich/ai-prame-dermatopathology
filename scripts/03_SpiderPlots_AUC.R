# ==========================================================================
# Project: AI PRAME Study - Morphological Features & Malignancy
# Script:  03_Analysis_Morphology_Performance_Advanced
# Author:  Frank Friedrich Gellrich
# Date:    2026-04-11
# Description: Spider plots (AUC) including One-vs-Rest analyses.
#              Extended to include PRAME Status 0 and 4.
#              Outputs: Formatted Word tables, Spider plots, and CSV data.
# ==========================================================================

# Part 1 - Packages ----
needed_packages <- c("tidyverse", "here", "pROC", "fmsb", "scales", "flextable", "officer", "magrittr")
invisible(lapply(needed_packages, library, character.only = TRUE))

# Part 2 - Configuration & Paths ----
path_cv80     <- here("data", "processed", "Master_Table_AI_Features_CV80.csv")
path_deploy20 <- here("data", "processed", "Master_Table_AI_Features_Deploy20.csv")
table_dir     <- here("results", "tables")
figure_dir    <- here("results", "figures")
csv_dir       <- here("results", "csv")

base_features <- c(
  "KI_maligne" = "Overall Malignancy", "Symmetrie" = "Symmetry",
  "Aufsteigende_Melanozyten" = "Upward Melanocytes", "Pagetoide_Durchsetzung" = "Pagetoid Spread",
  "Konfluenz_junktionaler_Nester" = "Nest Confluence", "Ausgepraegte_junktionale_Atypien" = "Junctional Atypia",
  "Atypische_junktionale_Mitosen" = "Junctional Mitoses", "Fehlende_Ausreifung_in_Nestern" = "Lack of Maturation",
  "Ausgepraegte_dermale_Atypien" = "Dermal Atypia", "Dermale_Mitosen" = "Dermal Mitoses",
  "Geordnete_Vereinzelung_Benignes_Muster" = "Ordered Single Cells",
  "Atypische_Vereinzelung_Dissoziation_Malignes_Muster" = "Atypical Dissociation",
  "Regression" = "Regression", "Ulzeration" = "Ulceration"
)

# Part 3 - Calculation Functions ----
df_cv <- read_csv(path_cv80, show_col_types = FALSE)
df_dp <- read_csv(path_deploy20, show_col_types = FALSE)

get_auc_ovr <- function(df, target_col, score_col, target_value) {
  clean <- df %>% select(all_of(c(target_col, score_col))) %>% drop_na()
  if(nrow(clean) < 10) return(NULL)
  
  # One-vs-Rest Labeling
  clean$binary_label <- ifelse(clean[[target_col]] == target_value, 1, 0)
  if(length(unique(clean$binary_label)) < 2) return(NULL)
  
  roc_obj <- pROC::roc(response = clean$binary_label, predictor = clean[[score_col]], 
                       levels = c(0, 1), direction = "<", quiet = TRUE)
  ci_val  <- pROC::ci.auc(roc_obj, method = "delong")
  
  tibble(N = nrow(clean), AUC = as.numeric(roc_obj$auc), Low = ci_val[1], High = ci_val[3])
}

results_list <- list()

# 1. Baseline features
for(f in names(base_features)) {
  results_list[[f]] <- tibble(
    English_Name = base_features[f],
    stats_cv = list(get_auc_ovr(df_cv, paste0(f, "_GroundTruth"), paste0(f, "_ja"), "ja")),
    stats_dp = list(get_auc_ovr(df_dp, paste0(f, "_GroundTruth"), paste0(f, "_ja"), "ja"))
  )
}

# 2. PRAME Status 0 & 4 (One-vs-Rest)
prame_scores <- c("0" = "PRAME Status 0", "4" = "PRAME Status 4")
for(s in names(prame_scores)) {
  id <- paste0("PRAME_", s)
  results_list[[id]] <- tibble(
    English_Name = prame_scores[s],
    stats_cv = list(get_auc_ovr(df_cv, "PRAME_Status_GroundTruth", paste0("PRAME_Status_", s), as.numeric(s))),
    stats_dp = list(get_auc_ovr(df_dp, "PRAME_Status_GroundTruth", paste0("PRAME_Status_", s), as.numeric(s)))
  )
}

# 3. Diagnostic categories
diag_types <- c("NZN" = "Diagnosis: Nevus", "Melanom_inSitu" = "Diagnosis: MiS", "Melanom_invasiv" = "Diagnosis: Inv. Melanoma")
for(d in names(diag_types)) {
  results_list[[d]] <- tibble(
    English_Name = diag_types[d],
    stats_cv = list(get_auc_ovr(df_cv, "KI_korr_Diagnose_GroundTruth", paste0("KI_korr_Diagnose_", d), d)),
    stats_dp = list(get_auc_ovr(df_dp, "KI_korr_Diagnose_GroundTruth", paste0("KI_korr_Diagnose_", d), d))
  )
}

# 4. Ki-67 categories
mib_types <- c("Niedrig" = "Ki-67: Low", "Intermediaer" = "Ki-67: Intermediate", "Hoch" = "Ki-67: High")
for(m in names(mib_types)) {
  results_list[[m]] <- tibble(
    English_Name = mib_types[m],
    stats_cv = list(get_auc_ovr(df_cv, "Mib1_GroundTruth", paste0("Mib1_", m), m)),
    stats_dp = list(get_auc_ovr(df_dp, "Mib1_GroundTruth", paste0("Mib1_", m), m))
  )
}

results_final <- bind_rows(results_list) %>%
  unnest_wider(stats_cv, names_sep = "_") %>%
  unnest_wider(stats_dp, names_sep = "_") %>%
  mutate(English_Label_N = paste0(English_Name, " (", stats_cv_N, "/", stats_dp_N, ")"))

# Part 4 - Spider Plots ----
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)

draw_spider_final <- function(data_subset, title, file_name) {
  plot_data <- data_subset %>%
    select(English_Label_N, stats_cv_AUC, stats_dp_AUC) %>%
    pivot_longer(cols = starts_with("stats"), names_to = "Set", values_to = "AUC") %>%
    pivot_wider(names_from = English_Label_N, values_from = AUC) %>%
    select(-Set)
  
  final_df <- as.data.frame(rbind(rep(1.0, ncol(plot_data)), rep(0.5, ncol(plot_data)), plot_data))
  
  png(file.path(figure_dir, file_name), width = 1600, height = 1400, res = 150)
  radarchart(final_df, axistype = 1, seg = 5,
             pcol = c("#2C3E50", "#E74C3C"), pfcol = c(alpha("#2C3E50", 0.1), alpha("#E74C3C", 0.1)), 
             plwd = 4, plty = 1, cglcol = "grey", cglty = 1, caxislabels = seq(0.5, 1, 0.1), 
             vlcex = 0.6, title = title)
  
  legend(x = "bottom", inset = c(0, -0.1), legend = c("Training Set (80%)", "Independent Test Set (20%)"), horiz = TRUE,
         bty = "n", pch = 20, col = c("#2C3E50", "#E74C3C"), pt.cex = 2, cex = 1, xpd = TRUE)
  dev.off()
}

draw_spider_final(results_final, "AUC Fingerprint: Extended Feature Analysis", "03_Spider_Extended_Analysis.png")
draw_spider_final(results_final %>% filter(English_Name %in% base_features), 
                  "AUC Fingerprint: Key Diagnostic Features", "03_Spider_Selected_Features.png")

# Part 5 - Flextables (Consistent Legend) ----
dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)

ft_data <- results_final %>%
  mutate(
    `AUC CV [95% CI]` = sprintf("%.3f [%.3f-%.3f]", stats_cv_AUC, stats_cv_Low, stats_cv_High),
    `AUC Deploy [95% CI]` = sprintf("%.3f [%.3f-%.3f]", stats_dp_AUC, stats_dp_Low, stats_dp_High)
  ) %>%
  select(Feature = English_Name, `N (CV)` = stats_cv_N, `AUC CV [95% CI]`, 
         `N (Deploy)` = stats_dp_N, `AUC Deploy [95% CI]`)

# Legend for Table 03
full_legend <- "Note: Performance metrics (AUC) are stratified by Training Set (80%) and Independent Test Set (20%). For multiclass variables such as PRAME status, pathological diagnosis, and Ki-67 index, category-specific AUCs were calculated using a One-vs-Rest (OvR) approach against all other categories. Confidence Intervals (CI) were determined using the DeLong method. Abbreviations: AUC, Area Under the Receiver Operating Characteristic Curve; CI, Confidence Interval."

ft <- flextable(ft_data) %>%
  set_caption(caption = "Table 03: Comprehensive AUC Performance Metrics across all diagnostic features.") %>%
  autofit() %>% theme_vanilla() %>% bold(part = "header") %>%
  add_footer_lines(full_legend) %>%
  fontsize(size = 9, part = "footer")

read_docx() %>% body_add_flextable(value = ft) %>% body_end_section_landscape() %>%
  print(target = file.path(table_dir, "03_Table_Performance_Comprehensive.docx"))

# Part 6 - CSV Export (Data for downstream Python visualizations) ----
dir.create(csv_dir, recursive = TRUE, showWarnings = FALSE)

csv_export_df <- results_final %>%
  select(
    Feature_Name = English_Name,
    Label_With_N = English_Label_N,
    CV_N = stats_cv_N,
    CV_AUC = stats_cv_AUC,
    CV_CI_Lower = stats_cv_Low,
    CV_CI_Upper = stats_cv_High,
    Deploy_N = stats_dp_N,
    Deploy_AUC = stats_dp_AUC,
    Deploy_CI_Lower = stats_dp_Low,
    Deploy_CI_Upper = stats_dp_High
  )

write_csv(csv_export_df, file.path(csv_dir, "03_Comprehensive_Performance_Metrics.csv"))

cat("\n--- Script 03 (Analysis Morphology Performance) successfully executed. ---\n")