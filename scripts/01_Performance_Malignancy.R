# ==========================================================================
# Project: AI PRAME Study - Morphological Features & Malignancy
# Script:  01_Performance_Malignancy
# Author:  Frank Friedrich Gellrich
# Date:    2026-04-11
# Description: Performance of the "Malignancy Prediction" (CV80 vs. Deploy20).
#              Generates ROC/PR plots, formatted tables, and CSV exports.
# ==========================================================================

# Part 1 - Packages ----
needed_packages <- c("tidyverse", "here", "pROC", "PRROC", "patchwork", "flextable", "officer")

load_packages <- function(pkgs) {
  new_pkgs <- pkgs[!(pkgs %in% installed.packages()[, "Package"])]
  if (length(new_pkgs)) install.packages(new_pkgs)
  invisible(lapply(pkgs, library, character.only = TRUE))
}
load_packages(needed_packages)

# Part 2 - Configuration & Paths ----
path_cv80     <- here("data", "processed", "Master_Table_AI_Features_CV80.csv")
path_deploy20 <- here("data", "processed", "Master_Table_AI_Features_Deploy20.csv")
table_dir     <- here("results", "tables")
figure_dir    <- here("results", "figures")
csv_dir       <- here("results", "csv")

# Part 3 - Data Preparation ----
df_cv80   <- read_csv(path_cv80, col_types = cols(Schnitt_ID = col_character()), show_col_types = FALSE)
df_deploy <- read_csv(path_deploy20, col_types = cols(Schnitt_ID = col_character()), show_col_types = FALSE)

prepare_metrics_data <- function(df, feature_name) {
  gt_col <- paste0(feature_name, "_GroundTruth")
  prob_col <- paste0(feature_name, "_ja")
  
  df %>%
    select(Schnitt_ID, !!sym(gt_col), !!sym(prob_col)) %>%
    drop_na() %>%
    mutate(gt_numeric = ifelse(!!sym(gt_col) == "ja" | !!sym(gt_col) == 1, 1, 0))
}

data_cv <- prepare_metrics_data(df_cv80, "KI_maligne")
data_dp <- prepare_metrics_data(df_deploy, "KI_maligne")

# Part 4 - Calculation of Metrics & CIs ----

# 1. ROC-AUC (pROC) incl. 95% DeLong CI
roc_cv <- pROC::roc(data_cv$gt_numeric, data_cv$KI_maligne_ja, quiet=TRUE)
roc_dp <- pROC::roc(data_dp$gt_numeric, data_dp$KI_maligne_ja, quiet=TRUE)

ci_cv <- pROC::ci.auc(roc_cv, method = "delong")
ci_dp <- pROC::ci.auc(roc_dp, method = "delong")

# 2. PR-AUC (PRROC)
pr_cv <- PRROC::pr.curve(scores.class0 = data_cv$KI_maligne_ja[data_cv$gt_numeric == 1],
                         scores.class1 = data_cv$KI_maligne_ja[data_cv$gt_numeric == 0], curve = TRUE)
pr_dp <- PRROC::pr.curve(scores.class0 = data_dp$KI_maligne_ja[data_dp$gt_numeric == 1],
                         scores.class1 = data_dp$KI_maligne_ja[data_dp$gt_numeric == 0], curve = TRUE)

# Part 5 - Creation of Visualization Labels ----

label_roc_cv <- sprintf("CV (80%%) (AUC: %.3f [%.3f-%.3f])", roc_cv$auc, ci_cv[1], ci_cv[3])
label_roc_dp <- sprintf("Deployment (20%%) (AUC: %.3f [%.3f-%.3f])", roc_dp$auc, ci_dp[1], ci_dp[3])

label_pr_cv  <- sprintf("CV (80%%) (AUC: %.3f)", pr_cv$auc.integral)
label_pr_dp  <- sprintf("Deployment (20%%) (AUC: %.3f)", pr_dp$auc.integral)

get_clean_curve_df <- function(roc_obj, pr_obj, label_roc, label_pr) {
  df_roc <- data.frame(x = 1 - roc_obj$specificities, y = roc_obj$sensitivities, type = "ROC", set = label_roc) %>% 
    distinct() %>% arrange(x, y)
  df_pr  <- data.frame(x = pr_obj$curve[,1], y = pr_obj$curve[,2], type = "PR", set = label_pr) %>% 
    distinct() %>% arrange(x, desc(y))
  bind_rows(df_roc, df_pr)
}

plot_data_roc <- bind_rows(
  get_clean_curve_df(roc_cv, pr_cv, label_roc_cv, label_pr_cv) %>% filter(type == "ROC"),
  get_clean_curve_df(roc_dp, pr_dp, label_roc_dp, label_pr_dp) %>% filter(type == "ROC")
)

plot_data_pr <- bind_rows(
  get_clean_curve_df(roc_cv, pr_cv, label_roc_cv, label_pr_cv) %>% filter(type == "PR"),
  get_clean_curve_df(roc_dp, pr_dp, label_roc_dp, label_pr_dp) %>% filter(type == "PR")
)

# Part 6 - Creation and Export of Plots ----

p_roc <- ggplot(plot_data_roc, aes(x, y, color = set)) +
  geom_step(size = 1.1) + geom_abline(linetype = "dashed", color = "grey60") +
  scale_color_manual(values = setNames(c("#2C3E50", "#E74C3C"), c(label_roc_cv, label_roc_dp))) +
  theme_minimal(base_size = 14) + 
  labs(title = "ROC Curves: Malignancy Prediction", x = "1 - Specificity", y = "Sensitivity", color = "Dataset")

p_pr <- ggplot(plot_data_pr, aes(x, y, color = set)) +
  geom_step(size = 1.1) +
  scale_color_manual(values = setNames(c("#2C3E50", "#E74C3C"), c(label_pr_cv, label_pr_dp))) +
  theme_minimal(base_size = 14) + 
  labs(title = "PR Curves: Malignancy Prediction", x = "Recall", y = "Precision", color = "Dataset")

final_plot <- (p_roc | p_pr) + plot_layout(guides = "collect") & 
  theme(legend.position = "bottom", legend.direction = "vertical")

dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)
ggsave(file.path(figure_dir, "01_Malignancy_Performance_Final.png"), plot = final_plot, width = 12, height = 8, dpi = 300)

# Part 7 - Export of Tables (flextable & Word) ----

performance_base_df <- tibble(
  Dataset = c("Cross-Validation (80%)", "Independent Deployment (20%)"),
  N = c(nrow(data_cv), nrow(data_dp)),
  `ROC-AUC [95% CI]` = c(
    sprintf("%.4f [%.4f-%.4f]", roc_cv$auc, ci_cv[1], ci_cv[3]),
    sprintf("%.4f [%.4f-%.4f]", roc_dp$auc, ci_dp[1], ci_dp[3])
  ),
  `PR-AUC` = c(round(pr_cv$auc.integral, 4), round(pr_dp$auc.integral, 4))
)

# Legend for Table 01
full_legend <- "Note: Diagnostic performance metrics for the holistic AI malignancy prediction score (baseline). Performance was evaluated on the training set (80%) via cross-validation and on an independent test set (20%) consisting of the most recent cases (historical split). AUC values are reported with 95% Confidence Intervals calculated using the DeLong method. Abbreviations: AUC, Area Under the Curve; CI, Confidence Interval."

dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)

ft <- flextable(performance_base_df) %>%
  set_caption(caption = "Table 01: Performance metrics for Malignancy Prediction (Baseline).") %>%
  autofit() %>% theme_vanilla() %>% bold(part = "header") %>%
  add_footer_lines(full_legend) %>%
  fontsize(size = 9, part = "footer")

read_docx() %>% 
  body_add_par("Baseline Performance: Malignancy Prediction", style = "heading 1") %>%
  body_add_flextable(value = ft) %>% 
  body_end_section_landscape() %>%
  print(target = file.path(table_dir, "01_Table_Malignancy_Performance_Baseline.docx"))

# Part 8 - CSV Export (Data for downstream Python visualizations) ----
dir.create(csv_dir, recursive = TRUE, showWarnings = FALSE)

# 1. Export metrics (including raw values for CI)
metrics_csv_df <- performance_base_df %>%
  mutate(
    ROC_AUC_value = c(roc_cv$auc, roc_dp$auc),
    ROC_AUC_CI_lower = c(ci_cv[1], ci_dp[1]),
    ROC_AUC_CI_upper = c(ci_cv[3], ci_dp[3])
  )
write_csv(metrics_csv_df, file.path(csv_dir, "01_Malignancy_Performance_Metrics.csv"))

# 2. Export complete curve data (ROC and PR)
curve_data_csv <- bind_rows(
  plot_data_roc %>% mutate(Dataset_Group = ifelse(grepl("CV", set), "CV_80", "Deploy_20")),
  plot_data_pr %>% mutate(Dataset_Group = ifelse(grepl("CV", set), "CV_80", "Deploy_20"))
)
write_csv(curve_data_csv, file.path(csv_dir, "01_Malignancy_Curve_Data.csv"))

cat("\n--- Script 01 (Performance Malignancy) successfully executed. ---\n")