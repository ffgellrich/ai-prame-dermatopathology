# ==========================================================================
# Project: AI PRAME Study
# Script:  08_Ground_Truth_Morphology_Benchmark
# Author:  Frank Friedrich Gellrich
# Date:    2026-04-11
# Description: Benchmarking: Human GT vs. AI-Morphology vs. Baseline.
#              Outputs: Detailed legends, comparative tables, Forest Plots,
#                       and CSV data for further visualization.
# ==========================================================================

# Part 1 - Packages ----
needed_packages <- c("tidyverse", "here", "pROC", "PRROC", "patchwork", "flextable", "officer", "broom", "scales")
invisible(lapply(needed_packages, library, character.only = TRUE))

# Part 2 - Configuration & Data ----
path_cv80     <- here("data", "processed", "Master_Table_AI_Features_CV80.csv")
path_deploy20 <- here("data", "processed", "Master_Table_AI_Features_Deploy20.csv")
figure_dir    <- here("results", "figures")
table_dir     <- here("results", "tables")
csv_dir       <- here("results", "csv")

dir.create(figure_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(table_dir, showWarnings = FALSE, recursive = TRUE)

morph_map <- c(
  "Symmetrie" = "Symmetry", "Aufsteigende_Melanozyten" = "Upward Melanocytes",
  "Pagetoide_Durchsetzung" = "Pagetoid Spread", "Konfluenz_junktionaler_Nester" = "Nest Confluence",
  "Ausgepraegte_junktionale_Atypien" = "Junctional Atypia", "Atypische_junktionale_Mitosen" = "Junctional Mitoses",
  "Fehlende_Ausreifung_in_Nestern" = "Lack of Maturation", "Ausgepraegte_dermale_Atypien" = "Dermal Atypia",
  "Dermale_Mitosen" = "Dermal Mitoses", "Geordnete_Vereinzelung_Benignes_Muster" = "Ordered Single Cells",
  "Atypische_Vereinzelung_Dissoziation_Malignes_Muster" = "Atypical Dissociation",
  "Regression" = "Regression", "Ulzeration" = "Ulceration"
)

features_13  <- names(morph_map)
cols_gt      <- paste0(features_13, "_GroundTruth")
cols_pred    <- paste0(features_13, "_ja")

df_train <- read_csv(path_cv80, show_col_types = FALSE)
df_test  <- read_csv(path_deploy20, show_col_types = FALSE)

prep_data <- function(df) {
  df %>%
    mutate(target = ifelse(KI_maligne_GroundTruth %in% c("ja", 1), 1, 0)) %>%
    mutate(across(all_of(cols_gt), ~ ifelse(. %in% c("ja", 1), 1, 0))) %>%
    select(target, all_of(cols_gt), all_of(cols_pred), KI_maligne_ja) %>%
    drop_na()
}

train_clean <- prep_data(df_train)
test_clean  <- prep_data(df_test)

# Part 3 - Modeling ----
mod_gt <- glm(as.formula(paste("target ~", paste(paste0("`", cols_gt, "`"), collapse = " + "))), 
              data = train_clean, family = binomial)
mod_ai <- glm(as.formula(paste("target ~", paste(paste0("`", cols_pred, "`"), collapse = " + "))), 
              data = train_clean, family = binomial)

# Part 4 - Performance Metrics ----
test_clean <- test_clean %>%
  mutate(prob_gt = predict(mod_gt, newdata = test_clean, type = "response"),
         prob_ai = predict(mod_ai, newdata = test_clean, type = "response"),
         prob_ho = KI_maligne_ja)

roc_gt <- pROC::roc(test_clean$target, test_clean$prob_gt, quiet = TRUE)
roc_ai <- pROC::roc(test_clean$target, test_clean$prob_ai, quiet = TRUE)
roc_ho <- pROC::roc(test_clean$target, test_clean$prob_ho, quiet = TRUE)
ci_gt <- pROC::ci.auc(roc_gt); ci_ai <- pROC::ci.auc(roc_ai); ci_ho <- pROC::ci.auc(roc_ho)

pr_gt <- PRROC::pr.curve(scores.class0 = test_clean$prob_gt[test_clean$target == 1], scores.class1 = test_clean$prob_gt[test_clean$target == 0], curve = TRUE)
pr_ai <- PRROC::pr.curve(scores.class0 = test_clean$prob_ai[test_clean$target == 1], scores.class1 = test_clean$prob_ai[test_clean$target == 0], curve = TRUE)
pr_ho <- PRROC::pr.curve(scores.class0 = test_clean$prob_ho[test_clean$target == 1], scores.class1 = test_clean$prob_ho[test_clean$target == 0], curve = TRUE)

l_roc_gt <- sprintf("Human Ground Truth (AUC: %.3f [%.3f-%.3f])", roc_gt$auc, ci_gt[1], ci_gt[3])
l_roc_ai <- sprintf("AI-Predicted Morphology (AUC: %.3f [%.3f-%.3f])", roc_ai$auc, ci_ai[1], ci_ai[3])
l_roc_ho <- sprintf("Malignancy Prediction (AUC: %.3f [%.3f-%.3f])", roc_ho$auc, ci_ho[1], ci_ho[3])

l_pr_gt <- sprintf("Human Ground Truth (AUC: %.3f)", pr_gt$auc.integral)
l_pr_ai <- sprintf("AI-Predicted Morphology (AUC: %.3f)", pr_ai$auc.integral)
l_pr_ho <- sprintf("Malignancy Prediction (AUC: %.3f)", pr_ho$auc.integral)

get_curve_df <- function(roc_obj, pr_obj, label_roc, label_pr) {
  df_roc <- data.frame(x = 1 - roc_obj$specificities, y = roc_obj$sensitivities, type = "ROC", model = label_roc) %>% distinct() %>% arrange(x, y)
  df_pr  <- data.frame(x = pr_obj$curve[,1], y = pr_obj$curve[,2], type = "PR", model = label_pr) %>% distinct() %>% arrange(x, desc(y))
  bind_rows(df_roc, df_pr)
}

plot_roc <- bind_rows(get_curve_df(roc_gt, pr_gt, l_roc_gt, l_pr_gt), get_curve_df(roc_ai, pr_ai, l_roc_ai, l_pr_ai), get_curve_df(roc_ho, pr_ho, l_roc_ho, l_pr_ho)) %>% filter(type == "ROC")
plot_pr  <- bind_rows(get_curve_df(roc_gt, pr_gt, l_roc_gt, l_pr_gt), get_curve_df(roc_ai, pr_ai, l_roc_ai, l_pr_ai), get_curve_df(roc_ho, pr_ho, l_roc_ho, l_pr_ho)) %>% filter(type == "PR")

p_roc <- ggplot(plot_roc, aes(x, y, color = model)) + geom_step(size = 1) + geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey") +
  scale_color_manual(values = setNames(c("#27AE60", "#E67E22", "#2C3E50"), c(l_roc_gt, l_roc_ai, l_roc_ho))) + theme_minimal() + labs(title = "ROC Comparison", x = "1 - Specificity", y = "Sensitivity", color = "Model")

p_pr <- ggplot(plot_pr, aes(x, y, color = model)) + geom_step(size = 1) + 
  scale_color_manual(values = setNames(c("#27AE60", "#E67E22", "#2C3E50"), c(l_pr_gt, l_pr_ai, l_pr_ho))) + theme_minimal() + labs(title = "PR Comparison", x = "Recall", y = "Precision", color = "Model")

final_curves <- (p_roc | p_pr) + plot_layout(guides = "collect") & theme(legend.position = "bottom", legend.direction = "vertical")
ggsave(file.path(figure_dir, "08_ROC_PR_Benchmark.png"), final_curves, width = 12, height = 8, dpi = 300, bg = "white")

# Part 5 - Forest Plot ----
get_stats <- function(model, label, suffix) {
  tidy(model, exponentiate = TRUE, conf.int = TRUE) %>%
    filter(term != "(Intercept)") %>%
    mutate(term_clean = str_remove_all(term, paste0("`|_ja|", suffix)),
           Feature = morph_map[term_clean],
           model_source = label)
}

stats_gt <- get_stats(mod_gt, "Human Ground Truth", "_GroundTruth")
stats_ai <- get_stats(mod_ai, "AI Prediction", "")

feature_order  <- stats_gt %>% arrange(estimate) %>% pull(Feature)
stats_combined <- bind_rows(stats_gt, stats_ai) %>% mutate(Feature = factor(Feature, levels = feature_order))

stats_plot <- stats_combined %>%
  mutate(is_sep = (Feature %in% c("Dermal Mitoses", "Ulceration") & model_source == "Human Ground Truth"),
         display_est = ifelse(is_sep | estimate > 500, 500, estimate),
         display_low = ifelse(is_sep, NA, ifelse(conf.low < 0.001, 0.001, conf.low)),
         display_high = ifelse(is_sep, NA, ifelse(conf.high > 500, 500, conf.high)),
         label_star = ifelse(is_sep, "*", ""))

p_forest <- ggplot(stats_plot, aes(x = display_est, y = Feature, color = model_source)) +
  geom_vline(xintercept = 1, linetype = "dashed", color = "grey30") +
  geom_errorbarh(aes(xmin = display_low, xmax = display_high), height = 0.3, position = position_dodge(width = 0.6), na.rm = TRUE) +
  geom_point(size = 2.5, position = position_dodge(width = 0.6)) +
  geom_text(aes(label = label_star), position = position_dodge(width = 0.6), hjust = -0.5, vjust = 0.4, size = 6, show.legend = FALSE) +
  scale_x_log10(limits = c(0.0008, 800), breaks = c(0.001, 0.01, 0.1, 1, 10, 100, 500), labels = c("0.001", "0.01", "0.1", "1", "10", "100", "500+")) +
  scale_color_manual(values = c("#E67E22", "#27AE60")) + theme_bw() + 
  labs(title = "OR Comparison: Human GT vs. AI-Morphology", subtitle = "* indicates perfect separation (OR -> Inf, CIs omitted)", x = "Odds Ratio (log scale)", y = "", color = "Data Source")

ggsave(file.path(figure_dir, "08_Forest_Benchmark.png"), p_forest, width = 11, height = 8, dpi = 300, bg = "white")

# Part 6 - Word Export with Legends ----
format_num <- function(x) { ifelse(abs(x) >= 1000 | (abs(x) < 0.001 & abs(x) > 0), formatC(x, format = "e", digits = 2), sprintf("%.2f", x)) }

final_table_df <- stats_gt %>% 
  select(Feature, est_gt=estimate, l_gt=conf.low, h_gt=conf.high, p_gt=p.value) %>%
  inner_join(stats_ai %>% select(Feature, est_ai=estimate, l_ai=conf.low, h_ai=conf.high, p_ai=p.value), by = "Feature") %>%
  mutate(GT_OR_CI = sprintf("%s [%s - %s]", format_num(est_gt), format_num(l_gt), format_num(h_gt)),
         GT_p = ifelse(p_gt < 0.001, "< .001", sprintf("%.3f", p_gt)),
         AI_OR_CI = sprintf("%s [%s - %s]", format_num(est_ai), format_num(l_ai), format_num(h_ai)),
         AI_p = ifelse(p_ai < 0.001, "< .001", sprintf("%.3f", p_ai))) %>%
  arrange(desc(est_gt))

# Legend for OR Table
or_legend <- "Note: Multivariate logistic regression comparing the predictive impact of individual features on malignancy. 'Human GT' refers to the model based on pathologist-defined ground truth; 'AI Pred' refers to the model based on AI-derived feature probabilities. Features with perfect separation in the human ground truth (where presence always correlates with malignancy) are marked with an asterisk (*), and their Confidence Intervals (CI) are omitted as the Odds Ratio (OR) approaches infinity. Bold p-values indicate statistical significance (p < 0.05)."

ft_or_wide <- flextable(final_table_df %>% select(Feature, GT_OR_CI, GT_p, AI_OR_CI, AI_p)) %>%
  set_header_labels(Feature = "Criterion", GT_OR_CI = "Human GT: OR [95% CI]", GT_p = "Human GT: p", AI_OR_CI = "AI Pred: OR [95% CI]", AI_p = "AI Pred: p") %>%
  theme_vanilla() %>% autofit() %>%
  bold(i = ~ final_table_df$p_gt < 0.05, j = "GT_p") %>%
  bold(i = ~ final_table_df$p_ai < 0.05, j = "AI_p") %>%
  add_footer_lines(or_legend)

# Legend for Performance Table
perf_legend <- "Note: Comparison of diagnostic performance metrics across different models. Human GT and AI-Predicted Morphology models each utilize 13 morphological features. Malignancy Prediction refers to the holistic baseline AI score. AUC values are presented with 95% Confidence Intervals (DeLong method)."

perf_df <- tibble(Model = c("Human Ground Truth (13 Features)", "AI-Predicted Morphology (13 Features)", "Malignancy Prediction (Baseline)"),
                  `ROC-AUC [95% CI]` = c(sprintf("%.3f [%.3f-%.3f]", roc_gt$auc, ci_gt[1], ci_gt[3]), sprintf("%.3f [%.3f-%.3f]", roc_ai$auc, ci_ai[1], ci_ai[3]), sprintf("%.3f [%.3f-%.3f]", roc_ho$auc, ci_ho[1], ci_ho[3])),
                  `PR-AUC` = c(round(pr_gt$auc.integral, 4), round(pr_ai$auc.integral, 4), round(pr_ho$auc.integral, 4)))

ft_perf <- flextable(perf_df) %>% 
  theme_vanilla() %>% autofit() %>% 
  bold(part = "header") %>%
  add_footer_lines(perf_legend)

# Word Export
doc <- read_docx() %>% body_add_par("Benchmark Results: Human vs. AI", style = "heading 1") %>%
  body_add_flextable(ft_perf) %>%
  body_add_break() %>% body_add_par("Comparison of Odds Ratios", style = "heading 2") %>%
  body_add_flextable(ft_or_wide) %>% body_end_section_landscape()

print(doc, target = file.path(table_dir, "08_Benchmark_Results_Final.docx"))

# Part 7 - CSV Export (Data for downstream Python visualizations) ----
dir.create(csv_dir, recursive = TRUE, showWarnings = FALSE)

# 1. Export performance metrics
write_csv(
  perf_df %>% mutate(
    ROC_AUC_value = c(roc_gt$auc, roc_ai$auc, roc_ho$auc),
    ROC_AUC_CI_lower = c(ci_gt[1], ci_ai[1], ci_ho[1]),
    ROC_AUC_CI_upper = c(ci_gt[3], ci_ai[3], ci_ho[3])
  ),
  file.path(csv_dir, "08_Benchmark_Performance_Metrics.csv")
)

# 2. Export Odds Ratios comparison (Human vs. AI)
write_csv(
  final_table_df %>% select(
    Criterion = Feature,
    Human_GT_OR = est_gt, Human_GT_CI_Lower = l_gt, Human_GT_CI_Upper = h_gt, Human_GT_p = p_gt, Human_GT_Formatted = GT_OR_CI,
    AI_Pred_OR = est_ai,  AI_Pred_CI_Lower = l_ai,  AI_Pred_CI_Upper = h_ai,  AI_Pred_p = p_ai,  AI_Pred_Formatted = AI_OR_CI
  ),
  file.path(csv_dir, "08_Benchmark_Odds_Ratios_Comparison.csv")
)

# 3. Export curve coordinates (ROC & PR)
curve_csv_df <- bind_rows(
  plot_roc %>% mutate(Model_Group = case_when(
    grepl("Human", model) ~ "Human_Ground_Truth",
    grepl("AI-Predicted", model) ~ "AI_Morphology",
    TRUE ~ "Baseline_Malignancy"
  )),
  plot_pr %>% mutate(Model_Group = case_when(
    grepl("Human", model) ~ "Human_Ground_Truth",
    grepl("AI-Predicted", model) ~ "AI_Morphology",
    TRUE ~ "Baseline_Malignancy"
  ))
)
write_csv(curve_csv_df, file.path(csv_dir, "08_Benchmark_Curve_Data.csv"))

cat("\n--- Script 08 (Ground Truth Morphology Benchmark) successfully executed. ---\n")