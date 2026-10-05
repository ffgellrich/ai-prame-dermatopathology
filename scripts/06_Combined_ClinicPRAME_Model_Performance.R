# ==========================================================================
# Project: AI PRAME Study
# Script:  06_Combined_ClinicPRAME_Model_Performance
# Author:  Frank Friedrich Gellrich
# Date:    2026-04-11
# Description: Triple Model including synthetic row for "PRAME 4 vs 0".
#              Visualizes the maximum clinical impact of PRAME.
#              Outputs: Formatted Word tables, plots, and CSV data.
# ==========================================================================

# Part 1 - Packages ----
needed_packages <- c("tidyverse", "here", "pROC", "PRROC", "patchwork", "flextable", "officer", "broom")
invisible(lapply(needed_packages, library, character.only = TRUE))

# Part 2 - Configuration & Data ----
path_cv80     <- here("data", "processed", "Master_Table_AI_Features_CV80.csv")
path_deploy20 <- here("data", "processed", "Master_Table_AI_Features_Deploy20.csv")
table_dir     <- here("results", "tables")
figure_dir    <- here("results", "figures")
csv_dir       <- here("results", "csv")

features_map <- c(
  "KI_maligne" = "Malignancy Prediction",
  "PRAME_Status_GroundTruth" = "PRAME Status (per point)",
  "Symmetrie" = "Symmetry (AI)",
  "Aufsteigende_Melanozyten" = "Upward Melanocytes (AI)",
  "Pagetoide_Durchsetzung" = "Pagetoid Spread (AI)",
  "Konfluenz_junktionaler_Nester" = "Nest Confluence (AI)",
  "Ausgepraegte_junktionale_Atypien" = "Junctional Atypia (AI)",
  "Atypische_junktionale_Mitosen" = "Junctional Mitoses (AI)",
  "Fehlende_Ausreifung_in_Nestern" = "Lack of Maturation (AI)",
  "Ausgepraegte_dermale_Atypien" = "Dermal Atypia (AI)",
  "Dermale_Mitosen" = "Dermal Mitoses (AI)",
  "Geordnete_Vereinzelung_Benignes_Muster" = "Ordered Single Cells (AI)",
  "Atypische_Vereinzelung_Dissoziation_Malignes_Muster" = "Atypical Dissociation (AI)",
  "Regression" = "Regression (AI)",
  "Ulzeration" = "Ulceration (AI)"
)

morph_features <- c("Symmetrie", "Aufsteigende_Melanozyten", "Pagetoide_Durchsetzung",
                    "Konfluenz_junktionaler_Nester", "Ausgepraegte_junktionale_Atypien", 
                    "Atypische_junktionale_Mitosen", "Fehlende_Ausreifung_in_Nestern", 
                    "Ausgepraegte_dermale_Atypien", "Dermale_Mitosen", 
                    "Geordnete_Vereinzelung_Benignes_Muster", 
                    "Atypische_Vereinzelung_Dissoziation_Malignes_Muster", 
                    "Regression", "Ulzeration")

predictor_cols <- c("KI_maligne_ja", paste0(morph_features, "_ja"), "PRAME_Status_GroundTruth")

df_train <- read_csv(path_cv80, col_types = cols(Schnitt_ID = col_character()), show_col_types = FALSE)
df_test  <- read_csv(path_deploy20, col_types = cols(Schnitt_ID = col_character()), show_col_types = FALSE)

prep_data <- function(df, cols) {
  df %>%
    mutate(target = ifelse(KI_maligne_GroundTruth %in% c("ja", 1), 1, 0)) %>%
    mutate(PRAME_Status_GroundTruth = as.numeric(as.character(PRAME_Status_GroundTruth))) %>%
    select(target, all_of(cols), KI_maligne_ja) %>%
    drop_na()
}

train_clean <- prep_data(df_train, predictor_cols)
test_clean  <- prep_data(df_test, predictor_cols)

# Part 3 - Modeling & Odds Ratios ----
formula_str   <- paste("target ~", paste(paste0("`", predictor_cols, "`"), collapse = " + "))
model_triple  <- glm(as.formula(formula_str), data = train_clean, family = binomial)

model_stats <- tidy(model_triple, exponentiate = TRUE, conf.int = TRUE) %>%
  filter(term != "(Intercept)") %>%
  mutate(
    term_clean = str_remove_all(term, "`|_ja"),
    Feature = features_map[term_clean]
  )

# Calculate synthetic row for PRAME 4 vs 0
prame_row <- model_stats %>% 
  filter(term_clean == "PRAME_Status_GroundTruth") %>%
  mutate(
    Feature = "PRAME Status (Score 4 vs 0)*",
    estimate = estimate^4,
    conf.low = conf.low^4,
    conf.high = conf.high^4,
    term_clean = "PRAME_4vs0"
  )

model_stats_final <- bind_rows(model_stats, prame_row) %>%
  mutate(label_text = sprintf("%.2f [%.2f - %.2f]", estimate, conf.low, conf.high)) %>%
  arrange(estimate) %>%
  mutate(Feature = factor(Feature, levels = Feature))

# Part 4 - Visualization: Forest Plot ----
p_main <- ggplot(model_stats_final, aes(x = estimate, y = Feature)) +
  geom_vline(xintercept = 1, linetype = "dashed", color = "grey30") +
  geom_errorbarh(aes(xmin = conf.low, xmax = conf.high), height = 0.2, color = "black") +
  geom_point(aes(size = ifelse(term_clean == "PRAME_4vs0", 4, 3))) +
  scale_x_log10(breaks = c(0.1, 1, 10, 100), labels = c("0.1", "1", "10", "100"), limits = c(0.001, 300)) +
  theme_bw() + theme(panel.grid.minor = element_blank(), axis.text = element_text(color = "black"), legend.position = "none") +
  labs(title = "Forest Plot: Predictors of Malignancy (Triple Model)", 
       subtitle = "*Synthetic row: Calculated effect for a PRAME score jump from 0 to 4",
       x = "Odds Ratio (log scale)", y = "")

p_table <- ggplot(model_stats_final, aes(y = Feature)) +
  geom_text(aes(x = 0, label = label_text), hjust = 0, size = 3.5, 
            fontface = ifelse(model_stats_final$term_clean == "PRAME_4vs0", "bold", "plain")) +
  theme_void() + theme(plot.margin = margin(l = 10, r = 10)) + labs(title = "OR [95% CI]")

final_forest <- (p_main | p_table) + plot_layout(widths = c(3, 1))

# Part 5 - Visualization: ROC/PR Curves ----
test_clean$prob_triple <- predict(model_triple, newdata = test_clean, type = "response")
test_clean$prob_base   <- test_clean$KI_maligne_ja

roc_triple <- pROC::roc(test_clean$target, test_clean$prob_triple, quiet=TRUE)
roc_base   <- pROC::roc(test_clean$target, test_clean$prob_base, quiet=TRUE)
ci_triple  <- pROC::ci.auc(roc_triple); ci_base <- pROC::ci.auc(roc_base)

pr_triple  <- PRROC::pr.curve(scores.class0 = test_clean$prob_triple[test_clean$target == 1], scores.class1 = test_clean$prob_triple[test_clean$target == 0], curve = TRUE)
pr_base    <- PRROC::pr.curve(scores.class0 = test_clean$prob_base[test_clean$target == 1], scores.class1 = test_clean$prob_base[test_clean$target == 0], curve = TRUE)

label_roc_triple <- sprintf("Triple Model (AUC: %.3f [%.3f-%.3f])", roc_triple$auc, ci_triple[1], ci_triple[3])
label_roc_base   <- sprintf("Malignancy Prediction (AUC: %.3f [%.3f-%.3f])", roc_base$auc, ci_base[1], ci_base[3])
label_pr_triple  <- sprintf("Triple Model (AUC: %.3f)", pr_triple$auc.integral)
label_pr_base    <- sprintf("Malignancy Prediction (AUC: %.3f)", pr_base$auc.integral)

get_curve_df <- function(roc_obj, pr_obj, label) {
  df_roc <- data.frame(x = 1 - roc_obj$specificities, y = roc_obj$sensitivities, type = "ROC", model = label) %>% distinct() %>% arrange(x, y)
  df_pr  <- data.frame(x = pr_obj$curve[,1], y = pr_obj$curve[,2], type = "PR", model = label) %>% distinct() %>% arrange(x, desc(y))
  bind_rows(df_roc, df_pr)
}

plot_data_roc <- bind_rows(get_curve_df(roc_triple, pr_triple, label_roc_triple), get_curve_df(roc_base, pr_base, label_roc_base))
plot_data_pr  <- bind_rows(get_curve_df(roc_triple, pr_triple, label_pr_triple), get_curve_df(roc_base, pr_base, label_pr_base))

p_roc <- ggplot(plot_data_roc %>% filter(type == "ROC"), aes(x, y, color = model)) +
  geom_step(size = 1) + geom_abline(linetype = "dashed", color = "grey") +
  scale_color_manual(values = setNames(c("#E74C3C", "#2C3E50"), c(label_roc_triple, label_roc_base))) +
  theme_minimal() + labs(title = "ROC Comparison", x = "1 - Specificity", y = "Sensitivity", color = "Model")

p_pr <- ggplot(plot_data_pr %>% filter(type == "PR"), aes(x, y, color = model)) +
  geom_step(size = 1) + scale_color_manual(values = setNames(c("#E74C3C", "#2C3E50"), c(label_pr_triple, label_pr_base))) +
  theme_minimal() + labs(title = "PR Comparison", x = "Recall", y = "Precision", color = "Model")

final_curves <- (p_roc | p_pr) + plot_layout(guides = "collect") & theme(legend.position = "bottom", legend.direction = "vertical")

# Part 6 - Export Figures ----
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)
ggsave(file.path(figure_dir, "06_Forest_Plot_Triple_Impact.png"), plot = final_forest, width = 12, height = 8, dpi = 300)
ggsave(file.path(figure_dir, "06_Model_Curves_Triple.png"), plot = final_curves, width = 12, height = 8, dpi = 300)

# Part 7 - Export Tables (Flextable & Word) ----
dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)

perf_legend <- "Note: Performance comparison between the holistic Malignancy Prediction baseline and the integrated Triple Model. AUC values include 95% Confidence Intervals (DeLong method)."
or_legend   <- "Note: Multivariate logistic regression of predictors in the Triple Model. The synthetic 'PRAME Score 4 vs 0' row estimates the cumulative effect size for a maximum PRAME score increase. Bold p-values indicate statistical significance (p < 0.05)."

performance_df_raw <- tibble(
  Model = c("Baseline (Malignancy Prediction)", "Triple Model (AI+Morph+PRAME)"),
  `ROC-AUC [95% CI]` = c(sprintf("%.4f [%.4f-%.4f]", roc_base$auc, ci_base[1], ci_base[3]), sprintf("%.4f [%.4f-%.4f]", roc_triple$auc, ci_triple[1], ci_triple[3])),
  `PR-AUC` = c(round(pr_base$auc.integral, 4), round(pr_triple$auc.integral, 4))
)

ft_perf <- performance_df_raw %>% flextable() %>% set_caption("Table 7: Triple Model Performance Metrics.") %>% 
  autofit() %>% theme_vanilla() %>% bold(part = "header") %>% add_footer_lines(perf_legend)

or_table_data <- model_stats_final %>% arrange(desc(estimate)) %>%
  mutate(`OR [95% CI]` = label_text, p_formatted = ifelse(p.value < 0.001, "< .001", sprintf("%.3f", p.value))) %>%
  select(Feature, `OR [95% CI]`, `p-value` = p_formatted, p_num = p.value, term_clean)

ft_or <- flextable(or_table_data %>% select(-p_num, -term_clean)) %>% 
  set_caption("Table 8: Odds Ratios for the Triple Model (including PRAME Score 4 vs 0 comparison).") %>% 
  autofit() %>% theme_vanilla() %>% bold(part = "header") %>%
  bold(i = ~ or_table_data$p_num < 0.05, j = "p-value") %>%
  bold(i = ~ or_table_data$term_clean == "PRAME_4vs0") %>% add_footer_lines(or_legend)

doc <- read_docx() %>% body_add_par("Supplementary Material: Triple Model Impact Analysis", style = "heading 1") %>%
  body_add_flextable(value = ft_perf) %>% body_add_break() %>%
  body_add_flextable(value = ft_or) %>% body_end_section_landscape()

print(doc, target = file.path(table_dir, "06_Triple_Model_Full_Results.docx"))

# Part 8 - CSV Export (Data for downstream Python visualizations) ----
dir.create(csv_dir, recursive = TRUE, showWarnings = FALSE)

# 1. Export performance comparison
performance_csv_df <- performance_df_raw %>%
  mutate(
    ROC_AUC_value = c(roc_base$auc, roc_triple$auc),
    ROC_AUC_CI_lower = c(ci_base[1], ci_triple[1]),
    ROC_AUC_CI_upper = c(ci_base[3], ci_triple[3])
  )
write_csv(performance_csv_df, file.path(csv_dir, "06_Triple_Model_Performance_Metrics.csv"))

# 2. Export model statistics (Odds Ratios) including the synthetic PRAME row
or_csv_df <- model_stats_final %>%
  select(
    Feature_Name = Feature,
    Variable_Clean = term_clean,
    Odds_Ratio = estimate,
    OR_CI_Lower = conf.low,
    OR_CI_Upper = conf.high,
    p_value = p.value, # For the synthetic row, this automatically reflects the baseline PRAME effect p-value
    Formatted_Text = label_text
  )
write_csv(or_csv_df, file.path(csv_dir, "06_Triple_Model_Odds_Ratios.csv"))

# 3. Export curve coordinates
curve_csv_df <- bind_rows(
  plot_data_roc %>% 
    filter(type == "ROC") %>% 
    mutate(Model_Group = ifelse(grepl("Triple", model), "Triple_Model", "Baseline_Malignancy")),
  plot_data_pr %>% 
    filter(type == "PR") %>% 
    mutate(Model_Group = ifelse(grepl("Triple", model), "Triple_Model", "Baseline_Malignancy"))
)
write_csv(curve_csv_df, file.path(csv_dir, "06_Triple_Model_Curve_Data.csv"))

cat("\n--- Script 06 (Combined Clinic/PRAME Model Performance) successfully executed. ---\n")