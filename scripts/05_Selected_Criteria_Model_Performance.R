# ==========================================================================
# Project: AI PRAME Study - Morphological Features & Malignancy
# Script:  05_Pure_Morphology_Model_Performance
# Author:  Frank Friedrich Gellrich
# Date:    2026-04-11
# Description: 1. Modeling based SOLELY on 13 morphological features
#              2. Comparison against Baseline (Malignancy Prediction)
#              3. ROC & PR curves including AUC values in the legend
#              4. Forest Plot & Word export (including detailed legends)
#              Outputs: Formatted Word tables, plots, and CSV data
# ==========================================================================

# Part 1 - Packages ----
needed_packages <- c("tidyverse", "here", "pROC", "PRROC", "patchwork", "flextable", "officer", "broom")
load_packages <- function(pkgs) {
  new_pkgs <- pkgs[!(pkgs %in% installed.packages()[, "Package"])]
  if (length(new_pkgs)) install.packages(new_pkgs)
  invisible(lapply(pkgs, library, character.only = TRUE))
}
load_packages(needed_packages)

# Part 2 - Configuration & Data ----
path_cv80     <- here("data", "processed", "Master_Table_AI_Features_CV80.csv")
path_deploy20 <- here("data", "processed", "Master_Table_AI_Features_Deploy20.csv")
table_dir     <- here("results", "tables")
figure_dir    <- here("results", "figures")
csv_dir       <- here("results", "csv")

features_map <- c(
  "Symmetrie" = "Symmetry", "Aufsteigende_Melanozyten" = "Upward Melanocytes",
  "Pagetoide_Durchsetzung" = "Pagetoid Spread", "Konfluenz_junktionaler_Nester" = "Nest Confluence",
  "Ausgepraegte_junktionale_Atypien" = "Junctional Atypia", "Atypische_junktionale_Mitosen" = "Junctional Mitoses",
  "Fehlende_Ausreifung_in_Nestern" = "Lack of Maturation", "Ausgepraegte_dermale_Atypien" = "Dermal Atypia",
  "Dermale_Mitosen" = "Dermal Mitoses", "Geordnete_Vereinzelung_Benignes_Muster" = "Ordered Single Cells",
  "Atypische_Vereinzelung_Dissoziation_Malignes_Muster" = "Atypical Dissociation",
  "Regression" = "Regression", "Ulzeration" = "Ulceration"
)

predictor_cols <- paste0(names(features_map), "_ja")

df_train <- read_csv(path_cv80, col_types = cols(Schnitt_ID = col_character()), show_col_types = FALSE)
df_test  <- read_csv(path_deploy20, col_types = cols(Schnitt_ID = col_character()), show_col_types = FALSE)

prep_data <- function(df, cols) {
  df %>%
    mutate(target = ifelse(KI_maligne_GroundTruth %in% c("ja", 1), 1, 0)) %>%
    select(target, all_of(cols), KI_maligne_ja) %>%
    drop_na()
}

train_clean <- prep_data(df_train, predictor_cols)
test_clean  <- prep_data(df_test, predictor_cols)

# Part 3 - Modeling & Odds Ratios ----
formula_str    <- paste("target ~", paste(paste0("`", predictor_cols, "`"), collapse = " + "))
model_morph    <- glm(as.formula(formula_str), data = train_clean, family = binomial)

model_stats <- tidy(model_morph, exponentiate = TRUE, conf.int = TRUE) %>%
  filter(term != "(Intercept)") %>%
  mutate(
    term_clean = str_remove_all(term, "`|_ja"),
    Feature = features_map[term_clean],
    label_text = sprintf("%.2f [%.2f - %.2f]", estimate, conf.low, conf.high)
  ) %>%
  arrange(estimate) %>%
  mutate(Feature = factor(Feature, levels = Feature))

# Part 4 - Visualization: Forest Plot ----
p_main <- ggplot(model_stats, aes(x = estimate, y = Feature)) +
  geom_vline(xintercept = 1, linetype = "dashed", color = "grey30") +
  geom_errorbarh(aes(xmin = conf.low, xmax = conf.high), height = 0.2, color = "black") +
  geom_point(size = 3) +
  scale_x_log10(breaks = c(0.01, 0.1, 1, 10, 100), labels = c("0.01", "0.1", "1", "10", "100"), limits = c(0.001, 250)) +
  theme_bw() + theme(panel.grid.minor = element_blank(), axis.text = element_text(color = "black")) +
  labs(title = "Forest Plot: Pure Morphological Features", x = "Odds Ratio (log scale)", y = "")

p_table <- ggplot(model_stats, aes(y = Feature)) +
  geom_text(aes(x = 0, label = label_text), hjust = 0, size = 3.5) +
  theme_void() + theme(plot.margin = margin(l = 10, r = 10)) + labs(title = "OR [95% CI]")

final_forest <- (p_main | p_table) + plot_layout(widths = c(3, 1))

# Part 5 - Visualization: ROC/PR Curves ----
test_clean$prob_morph    <- predict(model_morph, newdata = test_clean, type = "response")
test_clean$prob_baseline <- test_clean$KI_maligne_ja

roc_morph <- pROC::roc(test_clean$target, test_clean$prob_morph, quiet=TRUE)
roc_base  <- pROC::roc(test_clean$target, test_clean$prob_baseline, quiet=TRUE)
ci_morph  <- pROC::ci.auc(roc_morph)
ci_base   <- pROC::ci.auc(roc_base)

pr_morph  <- PRROC::pr.curve(scores.class0 = test_clean$prob_morph[test_clean$target == 1], scores.class1 = test_clean$prob_morph[test_clean$target == 0], curve = TRUE)
pr_base   <- PRROC::pr.curve(scores.class0 = test_clean$prob_baseline[test_clean$target == 1], scores.class1 = test_clean$prob_baseline[test_clean$target == 0], curve = TRUE)

label_roc_morph <- sprintf("Pure Morph. Model (AUC: %.3f [%.3f-%.3f])", roc_morph$auc, ci_morph[1], ci_morph[3])
label_roc_base  <- sprintf("Malignancy Prediction (AUC: %.3f [%.3f-%.3f])", roc_base$auc, ci_base[1], ci_base[3])
label_pr_morph  <- sprintf("Pure Morph. Model (AUC: %.3f)", pr_morph$auc.integral)
label_pr_base   <- sprintf("Malignancy Prediction (AUC: %.3f)", pr_base$auc.integral)

get_curve_df <- function(roc_obj, pr_obj, label) {
  df_roc <- data.frame(x = 1 - roc_obj$specificities, y = roc_obj$sensitivities, type = "ROC", model = label) %>% distinct() %>% arrange(x, y)
  df_pr  <- data.frame(x = pr_obj$curve[,1], y = pr_obj$curve[,2], type = "PR", model = label) %>% distinct() %>% arrange(x, desc(y))
  bind_rows(df_roc, df_pr)
}

plot_data_roc <- bind_rows(get_curve_df(roc_morph, pr_morph, label_roc_morph), get_curve_df(roc_base, pr_base, label_roc_base))
plot_data_pr  <- bind_rows(get_curve_df(roc_morph, pr_morph, label_pr_morph), get_curve_df(roc_base, pr_base, label_pr_base))

p_roc <- ggplot(plot_data_roc %>% filter(type == "ROC"), aes(x, y, color = model)) +
  geom_step(size = 1) + geom_abline(linetype = "dashed", color = "grey") +
  scale_color_manual(values = setNames(c("#E74C3C", "#2C3E50"), c(label_roc_morph, label_roc_base))) +
  theme_minimal() + labs(title = "ROC Comparison", x = "1 - Specificity", y = "Sensitivity", color = "Model")

p_pr <- ggplot(plot_data_pr %>% filter(type == "PR"), aes(x, y, color = model)) +
  geom_step(size = 1) + 
  scale_color_manual(values = setNames(c("#E74C3C", "#2C3E50"), c(label_pr_morph, label_pr_base))) +
  theme_minimal() + labs(title = "PR-Curve Comparison", x = "Recall", y = "Precision", color = "Model")

final_curves <- (p_roc | p_pr) + plot_layout(guides = "collect") & theme(legend.position = "bottom", legend.direction = "vertical")

# Part 6 - Export Figures ----
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)
ggsave(file.path(figure_dir, "05_Forest_Plot_Pure_Morphology.png"), plot = final_forest, width = 12, height = 7, dpi = 300)
ggsave(file.path(figure_dir, "05_Model_Comparison_Pure_Morphology.png"), plot = final_curves, width = 12, height = 8, dpi = 300)

# Part 7 - Export Tables (incl. Legends) ----
dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)

# Legend for Performance Table
perf_legend <- "Note: Performance comparison between the holistic Malignancy Prediction baseline and the 'Pure Morphological Model', which utilizes only the 13 listed AI-derived morphological features. AUC values include 95% Confidence Intervals (DeLong method)."

performance_df_raw <- tibble(
  Model = c("Baseline (Malignancy Prediction)", "Pure Morphological Model (13)"),
  `ROC-AUC [95% CI]` = c(sprintf("%.4f [%.4f-%.4f]", roc_base$auc, ci_base[1], ci_base[3]), sprintf("%.4f [%.4f-%.4f]", roc_morph$auc, ci_morph[1], ci_morph[3])),
  `PR-AUC` = c(round(pr_base$auc.integral, 4), round(pr_morph$auc.integral, 4))
)

ft_perf <- performance_df_raw %>% flextable() %>% 
  set_caption(caption = "Table 5: Performance metrics of the pure morphological model.") %>% 
  autofit() %>% theme_vanilla() %>% bold(part = "header") %>%
  add_footer_lines(perf_legend)

# Legend for OR Table
or_legend <- "Note: Multivariate logistic regression analysis of AI-derived morphological features as predictors of malignancy. Bold p-values indicate statistical significance (p < 0.05)."

or_table_data <- model_stats %>% arrange(desc(estimate)) %>%
  mutate(`OR [95% CI]` = label_text, p_num = p.value, `p-value` = ifelse(p.value < 0.001, "< .001", sprintf("%.3f", p.value))) %>%
  select(Feature, `OR [95% CI]`, `p-value`, p_num)

ft_or <- flextable(or_table_data %>% select(-p_num)) %>% 
  set_caption(caption = "Table 6: Odds Ratios for pure morphological features.") %>% 
  autofit() %>% theme_vanilla() %>% bold(part = "header") %>% 
  bold(i = ~ or_table_data$p_num < 0.05, j = "p-value") %>%
  add_footer_lines(or_legend)

doc <- read_docx() %>% 
  body_add_par("Study Results: Pure Morphological Model Performance", style = "heading 1") %>%
  body_add_flextable(value = ft_perf) %>% body_add_break() %>% 
  body_add_flextable(value = ft_or) %>% body_end_section_landscape()

print(doc, target = file.path(table_dir, "05_Pure_Morphology_Model_Results.docx"))

# Part 8 - CSV Export (Data for downstream Python visualizations) ----
dir.create(csv_dir, recursive = TRUE, showWarnings = FALSE)

# 1. Export performance comparison
performance_csv_df <- performance_df_raw %>%
  mutate(
    ROC_AUC_value = c(roc_base$auc, roc_morph$auc),
    ROC_AUC_CI_lower = c(ci_base[1], ci_morph[1]),
    ROC_AUC_CI_upper = c(ci_base[3], ci_morph[3])
  )
write_csv(performance_csv_df, file.path(csv_dir, "05_Pure_Morphology_Performance_Metrics.csv"))

# 2. Export model statistics (Odds Ratios)
or_csv_df <- model_stats %>%
  select(
    Feature_Name = Feature,
    Variable_Clean = term_clean,
    Odds_Ratio = estimate,
    OR_CI_Lower = conf.low,
    OR_CI_Upper = conf.high,
    p_value = p.value,
    Formatted_Text = label_text
  )
write_csv(or_csv_df, file.path(csv_dir, "05_Pure_Morphology_Odds_Ratios.csv"))

# 3. Export curve coordinates
curve_csv_df <- bind_rows(
  plot_data_roc %>% mutate(Model_Group = ifelse(grepl("Pure", model), "Pure_Morphology_Model", "Baseline_Malignancy")),
  plot_data_pr %>% mutate(Model_Group = ifelse(grepl("Pure", model), "Pure_Morphology_Model", "Baseline_Malignancy"))
)
write_csv(curve_csv_df, file.path(csv_dir, "05_Pure_Morphology_Curve_Data.csv"))

cat("\n--- Script 05 (Pure Morphology Model Performance) successfully executed. ---\n")