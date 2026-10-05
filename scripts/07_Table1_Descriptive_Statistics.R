# ==========================================================================
# Project: AI PRAME Study
# Script:  07_Table1_Descriptive_Statistics
# Author:  Frank Friedrich Gellrich
# Date:    2026-05-31
# Description: Publication-ready Table 1. 
#              Integrates demographic data, clinicopathological subtypes, 
#              and morphological baseline criteria.
#              Applies Benjamini-Hochberg (FDR) correction strictly 
#              separated within each hypothesis family (class-wise).
#              Outputs: Formatted Word table and CSV data.
# ==========================================================================

# Part 1 - Packages ----
needed_packages <- c("tidyverse", "readxl", "here", "flextable", "officer")
load_packages <- function(pkgs) {
  new_pkgs <- pkgs[!(pkgs %in% installed.packages()[, "Package"])]
  if (length(new_pkgs)) install.packages(new_pkgs)
  invisible(lapply(pkgs, library, character.only = TRUE))
}
load_packages(needed_packages)

# Part 2 - Paths & Mapping ----
path_cv80     <- here("data", "processed", "Master_Table_AI_Features_CV80.csv")
path_deploy20 <- here("data", "processed", "Master_Table_AI_Features_Deploy20.csv")
path_excel    <- here("data", "processed", "PRAME_Daten_v121_reduced.xlsx")
path_subtypes <- here("data", "processed", "PRAME_lesion_subtypes_n955.xlsx")
table_dir     <- here("results", "tables")
csv_dir       <- here("results", "csv")

recode_map <- list(
  Geschlecht = c("männlich" = "male", "weiblich" = "female"),
  diagnosis = c("NZN" = "Common Nevus", "Melanom_inSitu" = "Melanoma in situ", "Melanom_invasiv" = "Invasive melanoma"),
  KI_korr_Diagnose = c("NZN" = "Common Nevus", "Melanom_inSitu" = "Melanoma in situ", "Melanom_invasiv" = "Invasive melanoma"),
  KI_maligne = c("nein" = "benign", "ja" = "malignant"),
  Mib1 = c("Niedrig (<5%)" = "Low (<5%)", "Intermediär (5-20%)" = "Intermediate (5-20%)", "Hoch (>20%)" = "High (>20%)"),
  Lokalisation = c("Torso" = "Torso", "Kopf" = "Head/Neck", "UntereExtremitaet" = "Lower extremity", "ObereExtremitaet" = "Upper extremity", "Sonstige" = "Other locations")
)

# Part 3 - Data Import & Recoding ----
df_cv80_ids   <- read_csv(path_cv80, col_types = cols(Schnitt_ID = col_character()), show_col_types = FALSE) %>% select(Schnitt_ID) %>% mutate(Set = "Training Set (Historical 80%)")
df_deploy_ids <- read_csv(path_deploy20, col_types = cols(Schnitt_ID = col_character()), show_col_types = FALSE) %>% select(Schnitt_ID) %>% mutate(Set = "Independent Test Set (Most Recent 20%)")
id_mapping    <- bind_rows(df_cv80_ids, df_deploy_ids)

df_raw <- read_excel(path_excel) %>% mutate(Schnitt_ID = as.character(Schnitt_ID))
colnames(df_raw) <- trimws(colnames(df_raw))

df_subtypes_raw <- read_excel(path_subtypes) %>% mutate(Schnitt_ID = as.character(Schnitt_ID))
colnames(df_subtypes_raw) <- trimws(colnames(df_subtypes_raw))

df_subtypes_clean <- df_subtypes_raw %>%
  select(
    Schnitt_ID, NZN_blue, NZN_spitz, 
    NZN_irritiert, NZN_combined, NZN_reed, NZN_special, NZN_dysplastic, 
    Lentigo_maligna
  )

df_clean <- df_raw %>%
  inner_join(id_mapping, by = "Schnitt_ID") %>%
  left_join(df_subtypes_clean, by = "Schnitt_ID") %>%
  mutate(
    Set = factor(Set, levels = c("Training Set (Historical 80%)", "Independent Test Set (Most Recent 20%)")),
    Alter = as.numeric(Alter),
    Gender = recode(Geschlecht, !!!recode_map$Geschlecht),
    Localization = recode(`Lokalisation grob`, !!!recode_map$Lokalisation, .default = "Other locations"),
    `Original Diagnosis` = recode(diagnosis, !!!recode_map$diagnosis),
    `Corrected Diagnosis (Slide-level)` = recode(KI_korr_Diagnose, !!!recode_map$KI_korr_Diagnose),
    Malignancy = recode(KI_maligne, !!!recode_map$KI_maligne),
    `Ki-67 Index (Mib1)` = factor(recode(Mib1, !!!recode_map$Mib1), levels = c("Low (<5%)", "Intermediate (5-20%)", "High (>20%)")),
    `PRAME Status` = factor(PRAME_Status, levels = 0:4),
    `Patient Frequency` = factor(ifelse(unique_pat == 1, "Primary lesion", "Multiple specimens"), levels = c("Primary lesion", "Multiple specimens")),
    `Slide Density` = factor(ifelse(unique_case == 1, "Index level", "Additional levels"), levels = c("Index level", "Additional levels")),
    
    `Tumor Thickness (mm)` = as.numeric(`Tumordicke (mm)`),
    `Clark Level` = factor(Clark_Level, levels = 1:5, labels = c("I", "II", "III", "IV", "V")),
    
    `NZN Unspecified` = ifelse(
      `Corrected Diagnosis (Slide-level)` == "Common Nevus" &
        !(NZN_blue %in% 1 | NZN_spitz %in% 1 | NZN_irritiert %in% 1 | 
            NZN_combined %in% 1 | NZN_reed %in% 1 | NZN_special %in% 1 | 
            NZN_dysplastic %in% 1),
      1, 0
    ),
    
    `MIS Classic` = ifelse(
      `Corrected Diagnosis (Slide-level)` == "Melanoma in situ" & !(Lentigo_maligna %in% 1),
      1, 0
    )
  )

# Part 4 - Helper Functions for Table Generation ----
get_cat_row <- function(data, var_name, label, show_missing = TRUE) {
  sub <- data %>% select(Set, val = !!sym(var_name)) %>% mutate(val = as.character(val))
  if(show_missing) sub$val[is.na(sub$val)] <- "Missing / Unknown" else sub <- sub %>% filter(!is.na(val))
  if(nrow(sub) == 0) return(tibble(Characteristic = label, `Training Set` = "0 (0.0%)", `Test Set` = "0 (0.0%)", p_val_raw = NA))
  tab <- table(sub$Set, sub$val)
  p_val <- if(ncol(tab) > 1 && nrow(tab) > 1) suppressWarnings(chisq.test(tab)$p.value) else NA
  lvls <- sort(unique(sub$val))
  if("Missing / Unknown" %in% lvls) lvls <- c(lvls[lvls != "Missing / Unknown"], "Missing / Unknown")
  summary_stats <- sub %>% group_by(Set, val) %>% tally() %>% group_by(Set) %>% 
    mutate(perc = (n / sum(n)) * 100, text = sprintf("%d (%.1f%%)", n, perc))
  res_list <- lapply(lvls, function(lvl) {
    tibble(Characteristic = paste0("  ", lvl),
           `Training Set` = summary_stats$text[summary_stats$Set == "Training Set (Historical 80%)" & summary_stats$val == lvl] %||% "0 (0.0%)",
           `Test Set` = summary_stats$text[summary_stats$Set == "Independent Test Set (Most Recent 20%)" & summary_stats$val == lvl] %||% "0 (0.0%)",
           p_val_raw = NA)
  })
  bind_rows(tibble(Characteristic = label, `Training Set` = "", `Test Set` = "", p_val_raw = p_val), bind_rows(res_list))
}

get_binary_row <- function(data, var_name, label) {
  sub <- data %>% select(Set, val = !!sym(var_name)) %>% mutate(val = tolower(as.character(val))) %>% filter(!is.na(val))
  tab <- table(sub$Set, sub$val)
  p_val <- if(ncol(tab) > 1 && nrow(tab) > 1) suppressWarnings(chisq.test(tab)$p.value) else NA
  stats <- sub %>% group_by(Set) %>% summarise(n_total = n(), n_ja = sum(val == "ja", na.rm = TRUE), perc = (n_ja / n_total) * 100, text = sprintf("%d/%d (%.1f%%)", n_ja, n_total, perc))
  tibble(Characteristic = label, 
         `Training Set` = stats$text[stats$Set == "Training Set (Historical 80%)"] %||% "0/0 (0.0%)",
         `Test Set` = stats$text[stats$Set == "Independent Test Set (Most Recent 20%)"] %||% "0/0 (0.0%)",
         p_val_raw = p_val)
}

get_binary_subcohort_row <- function(data, var_name, label) {
  sub <- data %>% select(Set, val = !!sym(var_name)) %>% mutate(val = ifelse(val %in% 1, "yes", "no"))
  tab <- table(sub$Set, sub$val)
  p_val <- if("yes" %in% sub$val && ncol(tab) > 1 && nrow(tab) > 1) suppressWarnings(chisq.test(tab)$p.value) else NA
  stats <- sub %>% group_by(Set) %>% 
    summarise(
      n_total = n(), 
      n_ja = sum(val == "yes", na.rm = TRUE), 
      perc = (n_ja / n_total) * 100, 
      text = sprintf("%d (%.1f%%)", n_ja, perc),
      .groups = "drop"
    )
  tibble(Characteristic = paste0("    ", label), 
         `Training Set` = stats$text[stats$Set == "Training Set (Historical 80%)"] %||% "0 (0.0%)",
         `Test Set` = stats$text[stats$Set == "Independent Test Set (Most Recent 20%)"] %||% "0 (0.0%)",
         p_val_raw = p_val)
}

get_na_row <- function(data, var_name, label) {
  stats <- data %>% group_by(Set) %>% 
    summarise(
      n_total = n(),
      n_na = sum(is.na(!!sym(var_name))),
      perc = (n_na / n_total) * 100,
      text = sprintf("%d (%.1f%%)", n_na, perc),
      .groups = "drop"
    )
  tibble(Characteristic = label, 
         `Training Set` = stats$text[stats$Set == "Training Set (Historical 80%)"] %||% "0 (0.0%)",
         `Test Set` = stats$text[stats$Set == "Independent Test Set (Most Recent 20%)"] %||% "0 (0.0%)",
         p_val_raw = NA)
}

# Part 5 - Assembling the Table ----
df_nevi <- df_clean %>% filter(`Corrected Diagnosis (Slide-level)` == "Common Nevus")
df_mis  <- df_clean %>% filter(`Corrected Diagnosis (Slide-level)` == "Melanoma in situ")

table_list <- list(
  tibble(Characteristic = "Demographics & Clinical Background", `Training Set` = "", `Test Set` = "", p_val_raw = NA),
  tibble(Characteristic = "Age (years), Median [IQR]", 
         `Training Set` = sprintf("%.1f [%.1f-%.1f]", median(df_clean$Alter[df_clean$Set=="Training Set (Historical 80%)"], na.rm=T), quantile(df_clean$Alter[df_clean$Set=="Training Set (Historical 80%)"], 0.25, na.rm=T), quantile(df_clean$Alter[df_clean$Set=="Training Set (Historical 80%)"], 0.75, na.rm=T)), 
         `Test Set` = sprintf("%.1f [%.1f-%.1f]", median(df_clean$Alter[df_clean$Set=="Independent Test Set (Most Recent 20%)"], na.rm=T), quantile(df_clean$Alter[df_clean$Set=="Independent Test Set (Most Recent 20%)"], 0.25, na.rm=T), quantile(df_clean$Alter[df_clean$Set=="Independent Test Set (Most Recent 20%)"], 0.75, na.rm=T)), 
         p_val_raw = wilcox.test(Alter ~ Set, data=df_clean)$p.value),
  get_cat_row(df_clean, "Gender", "Gender"),
  get_cat_row(df_clean, "Localization", "Anatomical Localization"),
  get_cat_row(df_clean, "Patient Frequency", "Patient Presentation Status"),
  get_cat_row(df_clean, "Slide Density", "Specimen Sampling Density"),
  
  tibble(Characteristic = "Histopathology & Clinicopathological Subtypes", `Training Set` = "", `Test Set` = "", p_val_raw = NA),
  get_cat_row(df_clean, "Original Diagnosis", "Original Pathological Diagnosis"),
  get_cat_row(df_clean, "Corrected Diagnosis (Slide-level)", "Corrected Diagnosis (Slide-level)"),
  
  tibble(Characteristic = "  ↳ Nevus Subclassification (overlapping)", `Training Set` = "", `Test Set` = "", p_val_raw = NA),
  get_binary_subcohort_row(df_nevi, "NZN_blue", "Blue Nevus"),
  get_binary_subcohort_row(df_nevi, "NZN_spitz", "Spitz Nevus"),
  get_binary_subcohort_row(df_nevi, "NZN_irritiert", "Irritated Nevus"),
  get_binary_subcohort_row(df_nevi, "NZN_combined", "Combined Nevus"),
  get_binary_subcohort_row(df_nevi, "NZN_reed", "Reed Nevus"),
  get_binary_subcohort_row(df_nevi, "NZN_special", "Special Nevus"),
  get_binary_subcohort_row(df_nevi, "NZN_dysplastic", "Dysplastic Nevus"),
  get_binary_subcohort_row(df_nevi, "NZN Unspecified", "Common Nevus (Unspecified)"),
  
  tibble(Characteristic = "  ↳ Melanoma in situ Subclassification", `Training Set` = "", `Test Set` = "", p_val_raw = NA),
  get_binary_subcohort_row(df_mis, "Lentigo_maligna", "Lentigo Maligna"),
  get_binary_subcohort_row(df_mis, "MIS Classic", "Melanoma in situ (Classic)"),
  
  tibble(Characteristic = "  ↳ Invasive Melanoma Parameters", `Training Set` = "", `Test Set` = "", p_val_raw = NA),
  
  tibble(Characteristic = "    ↳ Tumor Thickness (mm), Median [IQR]", 
         `Training Set` = sprintf("%.2f [%.2f-%.2f]", 
                                  median(df_clean$`Tumor Thickness (mm)`[df_clean$Set=="Training Set (Historical 80%)" & df_clean$`Corrected Diagnosis (Slide-level)` == "Invasive melanoma"], na.rm=T), 
                                  quantile(df_clean$`Tumor Thickness (mm)`[df_clean$Set=="Training Set (Historical 80%)" & df_clean$`Corrected Diagnosis (Slide-level)` == "Invasive melanoma"], 0.25, na.rm=T), 
                                  quantile(df_clean$`Tumor Thickness (mm)`[df_clean$Set=="Training Set (Historical 80%)" & df_clean$`Corrected Diagnosis (Slide-level)` == "Invasive melanoma"], 0.75, na.rm=T)), 
         `Test Set` = sprintf("%.2f [%.2f-%.2f]", 
                              median(df_clean$`Tumor Thickness (mm)`[df_clean$Set=="Independent Test Set (Most Recent 20%)" & df_clean$`Corrected Diagnosis (Slide-level)` == "Invasive melanoma"], na.rm=T), 
                              quantile(df_clean$`Tumor Thickness (mm)`[df_clean$Set=="Independent Test Set (Most Recent 20%)" & df_clean$`Corrected Diagnosis (Slide-level)` == "Invasive melanoma"], 0.25, na.rm=T), 
                              quantile(df_clean$`Tumor Thickness (mm)`[df_clean$Set=="Independent Test Set (Most Recent 20%)" & df_clean$`Corrected Diagnosis (Slide-level)` == "Invasive melanoma"], 0.75, na.rm=T)), 
         p_val_raw = wilcox.test(`Tumor Thickness (mm)` ~ Set, data = df_clean %>% filter(`Corrected Diagnosis (Slide-level)` == "Invasive melanoma"))$p.value),
  
  get_cat_row(df_clean %>% filter(`Corrected Diagnosis (Slide-level)` == "Invasive melanoma"), "Clark Level", "    ↳ Clark Level (I-V)", show_missing = FALSE),
  get_na_row(df_clean %>% filter(`Corrected Diagnosis (Slide-level)` == "Invasive melanoma"), "Clark Level", "    ↳ Not documented (NA)"),
  
  get_cat_row(df_clean, "Malignancy", "Overall Malignancy Status"),
  get_cat_row(df_clean, "PRAME Status", "PRAME Immunohistochemistry (Score 0-4)"),
  get_cat_row(df_clean, "Ki-67 Index (Mib1)", "Proliferation Index (Ki-67/Mib1)"),
  
  tibble(Characteristic = "Morphological Criteria", `Training Set` = "", `Test Set` = "", p_val_raw = NA)
)

morph_map <- c("Junktionale_Beteiligung" = "Junctional Component", "Dermale_Beteiligung" = "Dermal Component", "Symmetrie" = "Symmetry", "Aufsteigende_Melanozyten" = "Upward Melanocytes", "Pagetoide_Durchsetzung" = "Pagetoid Spread", "Konfluenz_junktionaler_Nester" = "Nest Confluence", "Ausgepraegte_junktionale_Atypien" = "Junctional Atypia", "Atypische_junktionale_Mitosen" = "Junctional Mitoses", "Fehlende_Ausreifung_in_Nestern" = "Lack of Maturation", "Ausgepraegte_dermale_Atypien" = "Dermal Atypia", "Dermale_Mitosen" = "Dermal Mitoses", "Geordnete_Vereinzelung_Benignes_Muster" = "Ordered Single Cells", "Atypische_Vereinzelung_Dissoziation_Malignes_Muster" = "Atypical Dissociation", "Regression" = "Regression", "Dichtes_lymphozytaeres_Infiltrat" = "Dense Infiltrate", "Intratumorales_lymphozytaeres_Infiltrat_TILs" = "TILs", "Ulzeration" = "Ulceration", "Infiltration_der_Subkutis" = "Subcutis Infiltration")

for(i in seq_along(morph_map)) table_list[[length(table_list)+1]] <- get_binary_row(df_clean, names(morph_map)[i], morph_map[i])

# Part 6 - Multiple Testing Adjustment (Benjamini-Hochberg) ----
# Exact definition of hypothesis classes (vectors for separated adjustment)
nevi_rows  <- paste0("    ", c("Blue Nevus", "Spitz Nevus", "Irritated Nevus", "Combined Nevus", "Reed Nevus", "Special Nevus", "Dysplastic Nevus", "Common Nevus (Unspecified)"))
mis_rows   <- paste0("    ", c("Lentigo Maligna", "Melanoma in situ (Classic)"))
morph_rows <- unname(morph_map)

final_df <- bind_rows(table_list) %>% mutate(p_val_final = p_val_raw)

# Strictly isolated Benjamini-Hochberg correction (FDR) within each class
idx_nevi  <- final_df$Characteristic %in% nevi_rows
idx_mis   <- final_df$Characteristic %in% mis_rows
idx_morph <- final_df$Characteristic %in% morph_rows

final_df$p_val_final[idx_nevi]  <- p.adjust(final_df$p_val_raw[idx_nevi], method = "BH")
final_df$p_val_final[idx_mis]   <- p.adjust(final_df$p_val_raw[idx_mis], method = "BH")
final_df$p_val_final[idx_morph] <- p.adjust(final_df$p_val_raw[idx_morph], method = "BH")

# Final string formatting for Word export
final_df <- final_df %>%
  mutate(p_value = case_when(
    is.na(p_val_final) ~ "", 
    p_val_final < 0.001 ~ "< .001", 
    TRUE ~ sprintf("%.3f", p_val_final)
  ))

# Part 7 - Flextable Design & Legend ----
full_legend <- "Note: The study cohort was partitioned using a chronological (historical) split, with the most recent 20% of cases assigned to an independent test set. Continuous variables are reported as Median [Interquartile Range (IQR)]; categorical variables are reported as N (%). P-values were determined using the Wilcoxon rank-sum test for continuous data and Chi-square or Fisher's exact tests for categorical data. To control for type I error accumulation, multiple testing adjustments were performed using the Benjamini-Hochberg (False Discovery Rate, FDR) procedure, applied strictly and independently within each individual class of hypotheses (the nevus subclassification family, the melanoma in situ subclassification family, and the group of 18 morphological criteria, respectively). 'Corrected Diagnosis (Slide-level)' refers to a retrospective histopathological re-evaluation of the specific slide used for AI analysis, noting that in some cases only in situ components were present despite a case-level diagnosis of invasive melanoma. Subclassification percentages for nevi, melanomas in situ, and invasive parameters are calculated relative to their respective main diagnostic categories and may exceed 100% due to overlapping clinicopathological features. Bold p-values indicate statistical significance (p < 0.05 after class-wise FDR adjustment)."

ft <- flextable(final_df %>% select(-p_val_raw, -p_val_final)) %>%
  set_caption("Table 1: Baseline characteristics, clinicopathological subtypes, and morphological criteria of the study cohort.") %>%
  theme_vanilla() %>% autofit() %>%
  bold(i = ~ `Training Set` == "", j = "Characteristic", bold = TRUE) %>% 
  bold(i = ~ !is.na(as.numeric(str_replace(p_value, "< ", ""))) & as.numeric(str_replace(p_value, "< ", "")) < 0.05, 
       j = "p_value", bold = TRUE) %>%
  italic(i = ~ str_starts(Characteristic, "  ") & !str_detect(Characteristic, "↳"), italic = TRUE) %>%
  align(j = 2:4, align = "center", part = "all") %>%
  add_footer_lines(full_legend) %>%
  fontsize(size = 9, part = "footer")

# Word Export
dir.create(table_dir, showWarnings = FALSE)
doc <- read_docx() %>% body_add_flextable(value = ft) %>% body_end_section_landscape() %>%
  print(target = file.path(table_dir, "07_Table1_Descriptive_Statistics.docx"))

# Part 8 - CSV Export (Data for downstream Python visualizations) ----
dir.create(csv_dir, recursive = TRUE, showWarnings = FALSE)

# Saving the computed table including raw and adjusted p-values
write_csv(
  final_df %>% select(Characteristic, Training_Set = `Training Set`, Test_Set = `Test Set`, p_val_raw, p_val_FDR_adjusted = p_val_final, p_value_formatted = p_value),
  file.path(csv_dir, "07_Table1_Descriptive_Statistics.csv")
)

cat("\n--- Script 07 (Table 1 Descriptive Statistics) successfully executed. ---\n")