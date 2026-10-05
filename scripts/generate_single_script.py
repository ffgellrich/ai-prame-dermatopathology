# =============================================================================
# Project: AI PRAME Study - Morphological Features & Malignancy
# Script:  11_generate_publication_figure2.py
# Author:  Frank Friedrich Gellrich
# Date:    2026-04-11
# Description: Generates the publication 9-panel master grid (Figure 2)
#              integrating AUROC/AUPRC curves, radar plot, and forest plots.
#              Outputs: High-resolution PNG and vector PDF in results/figures/.
# =============================================================================

from pathlib import Path
import textwrap
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import seaborn as sns

# =============================================================================
# 1. PATH CONFIGURATION (RELATIVE TO PROJECT ROOT)
# =============================================================================
# Resolves root directory assuming this script resides in 'scripts/'
script_dir = Path(__file__).parent.resolve()
project_root = script_dir.parent

csv_dir = project_root / "results" / "csv"
output_dir = project_root / "results" / "figures"

output_path_png = output_dir / "Figure_2_Master_Grid.png"
output_path_pdf = output_dir / "Figure_2_Master_Grid.pdf"

# H&E (Hematoxylin-Eosin) inspired publication palette
ACTIVE_PALETTE = ["#C95C85", "#2C467A", "#338981", "#E08F44", "#785185", "#708090"]

# Contrast assignment for forest and curve plots
COLOR_CV = "#2C467A"      # Hematoxylin Blue (Human GT / CV 80% / Baseline)
COLOR_DEPLOY = "#C95C85"  # Eosin Pink (AI Prediction / Deploy 20% / Novel Method)


# =============================================================================
# 2. HELPER FUNCTIONS
# =============================================================================
def complete_curve_points(df_subset, curve_type, group_col):
    """Ensures boundary anchor points (0,0) and (1,1) for ROC or (1,0) for PR."""
    completed_groups = []
    for _, grp in df_subset.groupby(group_col, sort=False):
        grp = grp.sort_values(by=["x", "y"]).reset_index(drop=True)

        if curve_type == "ROC":
            if grp.iloc[0]["x"] > 0.0 or grp.iloc[0]["y"] > 0.0:
                start_row = grp.iloc[[0]].copy()
                start_row["x"], start_row["y"] = 0.0, 0.0
                grp = pd.concat([start_row, grp], ignore_index=True)
            if grp.iloc[-1]["x"] < 1.0 or grp.iloc[-1]["y"] < 1.0:
                end_row = grp.iloc[[-1]].copy()
                end_row["x"], end_row["y"] = 1.0, 1.0
                grp = pd.concat([grp, end_row], ignore_index=True)

        elif curve_type == "PR":
            if grp.iloc[-1]["x"] < 1.0 or grp.iloc[-1]["y"] > 0.0:
                end_row = grp.iloc[[-1]].copy()
                end_row["x"], end_row["y"] = 1.0, 0.0
                grp = pd.concat([grp, end_row], ignore_index=True)

        completed_groups.append(grp)
    return pd.concat(completed_groups, ignore_index=True) if completed_groups else df_subset


def format_feature_label(label, max_width=14):
    """Wraps feature labels for the polar radar plot."""
    clean_label = str(label).replace("_", " ")
    return "\n".join(textwrap.wrap(clean_label, width=max_width))


def ensure_synthetic_prame_row(df):
    """Ensures presence of the synthetic PRAME Score 4 vs 0 contrast row for Panel I."""
    df = df.copy()
    df["Feature_Name"] = df["Feature_Name"].astype(str).str.strip()
    is_synth_mask = df["Feature_Name"].str.contains(r"Score 4|4 vs 0", case=False, na=False)

    if not is_synth_mask.any():
        prame_mask = df["Feature_Name"].str.contains(r"PRAME", case=False, na=False)
        if prame_mask.any():
            p_row = df[prame_mask].iloc[0]
            or_4 = float(p_row["Odds_Ratio"]) ** 4
            ci_l_4 = float(p_row["OR_CI_Lower"]) ** 4
            ci_u_4 = float(p_row["OR_CI_Upper"]) ** 4

            synth_row = pd.DataFrame([{
                "Feature_Name": "PRAME Status (Score 4 vs 0)*",
                "Variable_Clean": "PRAME_Score_4_vs_0",
                "Odds_Ratio": or_4,
                "OR_CI_Lower": ci_l_4,
                "OR_CI_Upper": ci_u_4,
                "p_value": p_row.get("p_value", np.nan),
                "Formatted_Text": f"{or_4:.2f} [{ci_l_4:.2f} - {ci_u_4:.2f}]"
            }])
            df = pd.concat([df, synth_row], ignore_index=True)

    df["is_synthetic"] = df["Feature_Name"].str.contains(r"Score 4|4 vs 0|\*", case=False, na=False)
    return df


# =============================================================================
# 3. PANEL PLOTTING FUNCTIONS
# =============================================================================
def plot_panels_A_B(ax_roc, ax_pr):
    csv_path = csv_dir / "01_Malignancy_Curve_Data.csv"
    if not csv_path.exists():
        return
    df = pd.read_csv(csv_path)
    df["set"] = df["set"].astype(str).str.replace(r"\s*\[.*?\]", "", regex=True)

    df_roc = complete_curve_points(df[df["type"] == "ROC"].copy(), "ROC", "Dataset_Group")
    df_pr = complete_curve_points(df[df["type"] == "PR"].copy(), "PR", "Dataset_Group")

    palette = {
        df["Dataset_Group"].unique()[0]: COLOR_CV,
        df["Dataset_Group"].unique()[1]: COLOR_DEPLOY,
    }

    def clean_legend_label(lbl, curve_kind):
        lbl = lbl.replace("AUC", curve_kind).replace("Deployment", "Test Set (20%)").replace("Deploy", "Test Set (20%)")
        if "CV" in lbl and "Training" not in lbl:
            lbl = lbl.replace("CV", "Training Set (CV; 80%)")
        return lbl

    # Panel A: ROC
    ax_roc.set_box_aspect(1)
    ax_roc.plot([0, 1], [0, 1], color="#8F9499", linestyle="--", linewidth=1.0, label="Chance")
    sns.lineplot(data=df_roc, x="x", y="y", hue="Dataset_Group", palette=palette, linewidth=2.0, estimator=None, ax=ax_roc)
    ax_roc.set(xlabel="1 - Specificity (False Positive Rate)", ylabel="Sensitivity (True Positive Rate)", xlim=[-0.02, 1.02], ylim=[-0.02, 1.02])
    ax_roc.grid(True, linestyle=":", alpha=0.55, color="#c0c0c0")
    sns.despine(ax=ax_roc, top=True, right=True)

    handles, labels = ax_roc.get_legend_handles_labels()
    label_map = df_roc.drop_duplicates("Dataset_Group").set_index("Dataset_Group")["set"].to_dict()
    new_labels_roc = [clean_legend_label(label_map.get(l, l), "AUROC") for l in labels]
    ax_roc.legend(handles=handles, labels=new_labels_roc, loc="lower right", frameon=True, fontsize=8.0)

    # Panel B: PR
    ax_pr.set_box_aspect(1)
    sns.lineplot(data=df_pr, x="x", y="y", hue="Dataset_Group", palette=palette, linewidth=2.0, estimator=None, ax=ax_pr)
    ax_pr.set(xlabel="Recall (Sensitivity)", ylabel="Precision (Positive Predictive Value)", xlim=[-0.02, 1.02], ylim=[-0.02, 1.02])
    ax_pr.grid(True, linestyle=":", alpha=0.55, color="#c0c0c0")
    sns.despine(ax=ax_pr, top=True, right=True)

    handles_pr, labels_pr = ax_pr.get_legend_handles_labels()
    label_map_pr = df_pr.drop_duplicates("Dataset_Group").set_index("Dataset_Group")["set"].to_dict()
    new_labels_pr = [clean_legend_label(label_map_pr.get(l, l), "AUPRC") for l in labels_pr]
    ax_pr.legend(handles=handles_pr, labels=new_labels_pr, loc="lower left", frameon=True, fontsize=8.0)


def plot_panel_C(ax):
    csv_path = csv_dir / "03_Comprehensive_Performance_Metrics.csv"
    if not csv_path.exists():
        return
    df = pd.read_csv(csv_path)

    feature_labels = [format_feature_label(name) for name in df["Feature_Name"]]
    angles = np.linspace(0, 2 * np.pi, len(feature_labels), endpoint=False).tolist()
    angles_closed = angles + angles[:1]
    cv_closed = df["CV_AUC"].tolist() + [df["CV_AUC"].iloc[0]]
    deploy_closed = df["Deploy_AUC"].tolist() + [df["Deploy_AUC"].iloc[0]]

    ax.set_theta_offset(np.pi / 2)
    ax.set_theta_direction(-1)
    ax.set_xticks(angles)
    ax.set_xticklabels(feature_labels, fontsize=7.0, color="#222222")
    ax.tick_params(axis="x", pad=10)

    for angle, label in zip(angles, ax.get_xticklabels()):
        if np.isclose(angle, 0, atol=0.08):
            label.set_horizontalalignment("center")
            label.set_verticalalignment("bottom")
        elif np.isclose(angle, np.pi, atol=0.08):
            label.set_horizontalalignment("center")
            label.set_verticalalignment("top")
        elif 0 < angle < np.pi:
            label.set_horizontalalignment("left")
            label.set_verticalalignment("center")
        else:
            label.set_horizontalalignment("right")
            label.set_verticalalignment("center")

    ax.set_ylim(0.50, 1.00)
    ax.set_yticks([0.60, 0.70, 0.80, 0.90, 1.00])
    ax.set_yticklabels(["0.60", "0.70", "0.80", "0.90", "1.00"], fontsize=6.5, color="#666666")
    ax.set_rlabel_position(38)
    ax.grid(True, linestyle=":", alpha=0.6, color="#b0b0b0")
    ax.spines["polar"].set_color("#cccccc")

    ax.plot(np.linspace(0, 2 * np.pi, 100), [0.50] * 100, color="#8F9499", linestyle="--", linewidth=0.9, label="Chance (0.50)")
    ax.plot(angles_closed, cv_closed, color=COLOR_CV, linewidth=1.8, marker="o", markersize=3.5, label="Training Set (CV; 80%)")
    ax.fill(angles_closed, cv_closed, color=COLOR_CV, alpha=0.14)
    ax.plot(angles_closed, deploy_closed, color=COLOR_DEPLOY, linewidth=1.8, marker="s", markersize=3.5, label="Test Set (20%)")
    ax.fill(angles_closed, deploy_closed, color=COLOR_DEPLOY, alpha=0.14)

    ax.legend(loc="upper center", bbox_to_anchor=(0.5, -0.22), ncol=3, frameon=True, fontsize=8.0, columnspacing=1.5)


def plot_panels_D_E(ax_roc, ax_pr):
    csv_path = csv_dir / "08_Benchmark_Curve_Data.csv"
    if not csv_path.exists():
        return
    df = pd.read_csv(csv_path)
    df["model"] = df["model"].astype(str).str.replace(r"\s*\[.*?\]", "", regex=True)

    # Filter out slide-level baseline curves to focus strictly on morphological models
    mask_morphology = ~df["model"].str.contains(r"Malignancy|Slide|Score", case=False, na=False)
    df = df[mask_morphology].copy()

    df_roc = complete_curve_points(df[df["type"] == "ROC"].copy(), "ROC", "Model_Group")
    df_pr = complete_curve_points(df[df["type"] == "PR"].copy(), "PR", "Model_Group")

    palette = {}
    for mg in df["Model_Group"].unique():
        mg_str = str(mg).lower()
        if "human" in mg_str or "expert" in mg_str or "gt" in mg_str or "ground" in mg_str:
            palette[mg] = COLOR_CV
        else:
            palette[mg] = COLOR_DEPLOY

    # Panel D: ROC
    ax_roc.set_box_aspect(1)
    ax_roc.plot([0, 1], [0, 1], color="#8F9499", linestyle="--", linewidth=1.0, label="Chance")
    sns.lineplot(data=df_roc, x="x", y="y", hue="Model_Group", palette=palette, linewidth=2.0, estimator=None, ax=ax_roc)
    ax_roc.set(xlabel="1 - Specificity (False Positive Rate)", ylabel="Sensitivity (True Positive Rate)", xlim=[-0.02, 1.02], ylim=[-0.02, 1.02])
    ax_roc.grid(True, linestyle=":", alpha=0.55, color="#c0c0c0")
    sns.despine(ax=ax_roc, top=True, right=True)

    handles, labels = ax_roc.get_legend_handles_labels()
    label_map = df_roc.drop_duplicates("Model_Group").set_index("Model_Group")["model"].to_dict()
    new_labels_roc = [
        label_map.get(l, l).replace("AUC", "AUROC")
        .replace("Human Ground Truth", "Expert-Defined Morphology")
        .replace("AI Prediction", "AI-Predicted Morphology")
        for l in labels
    ]
    ax_roc.legend(handles=handles, labels=new_labels_roc, loc="lower right", frameon=True, fontsize=8.0)

    # Panel E: PR
    ax_pr.set_box_aspect(1)
    sns.lineplot(data=df_pr, x="x", y="y", hue="Model_Group", palette=palette, linewidth=2.0, estimator=None, ax=ax_pr)
    ax_pr.set(xlabel="Recall (Sensitivity)", ylabel="Precision (Positive Predictive Value)", xlim=[-0.02, 1.02], ylim=[-0.02, 1.02])
    ax_pr.grid(True, linestyle=":", alpha=0.55, color="#c0c0c0")
    sns.despine(ax=ax_pr, top=True, right=True)

    handles_pr, labels_pr = ax_pr.get_legend_handles_labels()
    label_map_pr = df_pr.drop_duplicates("Model_Group").set_index("Model_Group")["model"].to_dict()
    new_labels_pr = [
        label_map_pr.get(l, l).replace("AUC", "AUPRC")
        .replace("Human Ground Truth", "Expert-Defined Morphology")
        .replace("AI Prediction", "AI-Predicted Morphology")
        for l in labels_pr
    ]
    ax_pr.legend(handles=handles_pr, labels=new_labels_pr, loc="lower left", frameon=True, fontsize=8.0)


def plot_panel_F(ax):
    csv_path = csv_dir / "08_Benchmark_Odds_Ratios_Comparison.csv"
    if not csv_path.exists():
        return
    df = pd.read_csv(csv_path)
    df["Criterion"] = df["Criterion"].astype(str).str.strip()

    feature_order = [
        "Ordered Single Cells", "Symmetry", "Nest Confluence", "Atypical Dissociation",
        "Lack of Maturation", "Regression", "Junctional Mitoses", "Dermal Atypia",
        "Pagetoid Spread", "Junctional Atypia", "Upward Melanocytes", "Ulceration", "Dermal Mitoses"
    ]
    present_features = [f for f in feature_order if f in df["Criterion"].values]
    y_indices = {feat: idx for idx, feat in enumerate(present_features)}

    ax.set_xscale("log")
    ax.axvline(1.0, color="#777777", linestyle="--", linewidth=1.0, zorder=1)
    offset = 0.16

    for is_human, color, y_shift in [(True, COLOR_CV, offset), (False, COLOR_DEPLOY, -offset)]:
        first_entry = True
        for _, row in df.iterrows():
            if row["Criterion"] not in y_indices:
                continue
            y_pos = y_indices[row["Criterion"]] + y_shift
            or_val = row["Human_GT_OR"] if is_human else row["AI_Pred_OR"]
            ci_low = row["Human_GT_CI_Lower"] if is_human else row["AI_Pred_CI_Lower"]
            ci_high = row["Human_GT_CI_Upper"] if is_human else row["AI_Pred_CI_Upper"]
            label_text = "Expert" if is_human else "AI Prediction"

            if is_human and (or_val > 500 or pd.isna(ci_high) or ci_high > 1e4):
                ax.plot(500.0, y_pos, marker="o", markersize=5.0, color=color, label=label_text if first_entry else "", zorder=3)
                ax.text(620.0, y_pos - 0.05, "*", color=color, fontsize=11.5, fontweight="bold", va="center", zorder=3)
            else:
                err_low = max(0.0, or_val - max(1e-4, ci_low))
                err_high = max(0.0, ci_high - or_val)
                ax.errorbar(
                    x=or_val, y=y_pos, xerr=[[err_low], [err_high]], fmt="o", color=color,
                    markersize=5.0, capsize=2.8, elinewidth=1.1, label=label_text if first_entry else "", zorder=3
                )
            first_entry = False

    ax.set(xlim=[0.0005, 1200], xlabel="Odds Ratio (log scale)")
    ax.set_xticks([0.001, 0.01, 0.1, 1, 10, 100, 500])
    ax.set_xticklabels(["0.001", "0.01", "0.1", "1", "10", "100", "500+"], fontsize=8.5)
    ax.set_yticks(range(len(present_features)))
    ax.set_yticklabels(present_features, fontsize=8.2)
    ax.set_ylim([-0.6, len(present_features) - 0.4])
    ax.grid(True, which="major", linestyle=":", alpha=0.55, color="#c0c0c0")
    sns.despine(ax=ax, top=True, right=True)
    ax.legend(loc="lower right", bbox_to_anchor=(0.98, 0.02), frameon=True, fontsize=8.0, title="Data Source", title_fontsize=8.5)


def plot_panels_G_H(ax_roc, ax_pr):
    csv_path = csv_dir / "06_Triple_Model_Curve_Data.csv"
    if not csv_path.exists():
        return
    df = pd.read_csv(csv_path)
    df["model"] = df["model"].astype(str).str.replace(r"\s*\[.*?\]", "", regex=True)

    df_roc = complete_curve_points(df[df["type"] == "ROC"].copy(), "ROC", "Model_Group")
    df_pr = complete_curve_points(df[df["type"] == "PR"].copy(), "PR", "Model_Group")
    palette = dict(zip(df["Model_Group"].unique(), ACTIVE_PALETTE))

    # Panel G: ROC
    ax_roc.set_box_aspect(1)
    ax_roc.plot([0, 1], [0, 1], color="#8F9499", linestyle="--", linewidth=1.0, label="Chance")
    sns.lineplot(data=df_roc, x="x", y="y", hue="Model_Group", palette=palette, linewidth=2.0, estimator=None, ax=ax_roc)
    ax_roc.set(xlabel="1 - Specificity (False Positive Rate)", ylabel="Sensitivity (True Positive Rate)", xlim=[-0.02, 1.02], ylim=[-0.02, 1.02])
    ax_roc.grid(True, linestyle=":", alpha=0.55, color="#c0c0c0")
    sns.despine(ax=ax_roc, top=True, right=True)

    handles, labels = ax_roc.get_legend_handles_labels()
    label_map = df_roc.drop_duplicates("Model_Group").set_index("Model_Group")["model"].to_dict()
    new_labels_roc = [label_map.get(l, l).replace("AUC", "AUROC") for l in labels]
    ax_roc.legend(handles=handles, labels=new_labels_roc, loc="lower right", frameon=True, fontsize=8.0)

    # Panel H: PR
    ax_pr.set_box_aspect(1)
    sns.lineplot(data=df_pr, x="x", y="y", hue="Model_Group", palette=palette, linewidth=2.0, estimator=None, ax=ax_pr)
    ax_pr.set(xlabel="Recall (Sensitivity)", ylabel="Precision (Positive Predictive Value)", xlim=[-0.02, 1.02], ylim=[-0.02, 1.02])
    ax_pr.grid(True, linestyle=":", alpha=0.55, color="#c0c0c0")
    sns.despine(ax=ax_pr, top=True, right=True)

    handles_pr, labels_pr = ax_pr.get_legend_handles_labels()
    label_map_pr = df_pr.drop_duplicates("Model_Group").set_index("Model_Group")["model"].to_dict()
    new_labels_pr = [label_map_pr.get(l, l).replace("AUC", "AUPRC") for l in labels_pr]
    ax_pr.legend(handles=handles_pr, labels=new_labels_pr, loc="lower left", frameon=True, fontsize=8.0)


def plot_panel_I(ax):
    csv_path = csv_dir / "06_Triple_Model_Odds_Ratios.csv"
    if not csv_path.exists():
        return
    df = ensure_synthetic_prame_row(pd.read_csv(csv_path))

    df_reg = df[~df["is_synthetic"]].sort_values(by="Odds_Ratio", ascending=True).copy()
    df_synth = df[df["is_synthetic"]].copy()
    df_sorted = pd.concat([df_reg, df_synth], ignore_index=True)
    n_rows, n_regular = len(df_sorted), len(df_reg)

    ax.set_xscale("log")
    ax.axvline(1.0, color="#777777", linestyle="--", linewidth=1.0, zorder=1)
    if n_rows > n_regular:
        ax.axhline(n_regular - 0.5, color="#d0d0d0", linestyle=":", linewidth=0.9, zorder=1)

    for idx, row in df_sorted.iterrows():
        is_syn = row["is_synthetic"]
        point_color, point_marker, marker_size, line_width = (COLOR_DEPLOY, "D", 6.0, 1.4) if is_syn else (COLOR_CV, "o", 5.0, 1.1)
        err_low = max(0.0, row["Odds_Ratio"] - row["OR_CI_Lower"])
        err_high = max(0.0, row["OR_CI_Upper"] - row["Odds_Ratio"])
        ax.errorbar(
            x=row["Odds_Ratio"], y=idx, xerr=[[err_low], [err_high]], fmt=point_marker,
            color=point_color, markersize=marker_size, capsize=2.8, elinewidth=line_width, zorder=3
        )

    ax.set(xlim=[0.05, 80], xlabel="Odds Ratio (log scale)")
    ax.set_xticks([0.1, 1, 10])
    ax.set_xticklabels(["0.1", "1", "10"], fontsize=8.5)
    ax.set_yticks(range(n_rows))
    ax.set_yticklabels(df_sorted["Feature_Name"], fontsize=8.2)
    ax.set_ylim([-0.6, n_rows - 0.4])
    ax.grid(True, which="major", linestyle=":", alpha=0.55, color="#c0c0c0")
    sns.despine(ax=ax, top=True, right=True)


# =============================================================================
# 4. MASTER GRID ORCHESTRATION
# =============================================================================
def main():
    output_dir.mkdir(parents=True, exist_ok=True)
    sns.set_theme(style="ticks", context="paper")

    # Dimensions: 13.8 x 14.0 inches
    fig = plt.figure(figsize=(13.8, 14.0))

    # Outer layout specification
    outer_gs = gridspec.GridSpec(
        nrows=1, ncols=2, figure=fig,
        width_ratios=[2.0, 1.25],
        wspace=0.14, left=0.055, right=0.975, top=0.925, bottom=0.045
    )

    # Left grid (ROC & PR curves)
    gs_left = gridspec.GridSpecFromSubplotSpec(
        nrows=3, ncols=2, subplot_spec=outer_gs[0],
        wspace=0.24, hspace=0.42
    )

    # Right grid (Radar plot & Forest plots)
    gs_right = gridspec.GridSpecFromSubplotSpec(
        nrows=3, ncols=2, subplot_spec=outer_gs[1],
        width_ratios=[0.55, 1.0], wspace=0.0,
        height_ratios=[1.0, 1.0, 1.0], hspace=0.55
    )

    # ==================== ROW 1 ====================
    ax_a = fig.add_subplot(gs_left[0, 0])
    ax_b = fig.add_subplot(gs_left[0, 1])
    ax_c = fig.add_subplot(gs_right[0, :], projection="polar")

    plot_panels_A_B(ax_a, ax_b)
    ax_a.set_title("A   Slide-Level Malignancy (AUROC)", loc="left", fontsize=11.5, fontweight="bold", pad=8)
    ax_b.set_title("B   Slide-Level Malignancy (AUPRC)", loc="left", fontsize=11.5, fontweight="bold", pad=8)

    plot_panel_C(ax_c)
    ax_c.set_title("C   Morphological Feature AUROCs", loc="left", fontsize=11.5, fontweight="bold", pad=14)

    # ==================== ROW 2 ====================
    ax_d = fig.add_subplot(gs_left[1, 0])
    ax_e = fig.add_subplot(gs_left[1, 1])

    ax_f_title = fig.add_subplot(gs_right[1, :])
    ax_f_title.axis("off")
    ax_f_title.set_title("F   Feature Odds Ratios (Expert vs. AI)", loc="left", fontsize=11.5, fontweight="bold", pad=22)
    ax_f_title.text(
        0.0, 1.022, "* indicates near-perfect class separation (OR > 10⁶, CIs omitted)",
        transform=ax_f_title.transAxes, fontsize=7.5, fontstyle="italic", color="#555555", va="bottom"
    )

    ax_f = fig.add_subplot(gs_right[1, 1])

    plot_panels_D_E(ax_d, ax_e)
    ax_d.set_title("D   Diagnostic Benchmark (AUROC)", loc="left", fontsize=11.5, fontweight="bold", pad=8)
    ax_e.set_title("E   Diagnostic Benchmark (AUPRC)", loc="left", fontsize=11.5, fontweight="bold", pad=8)
    plot_panel_F(ax_f)

    # ==================== ROW 3 ====================
    ax_g = fig.add_subplot(gs_left[2, 0])
    ax_h = fig.add_subplot(gs_left[2, 1])

    ax_i_title = fig.add_subplot(gs_right[2, :])
    ax_i_title.axis("off")
    ax_i_title.set_title("I   Predictors of Malignancy (Triple Model)", loc="left", fontsize=11.5, fontweight="bold", pad=22)
    ax_i_title.text(
        0.0, 1.022, "* Pre-specified contrast: PRAME score 4 vs 0",
        transform=ax_i_title.transAxes, fontsize=7.5, fontstyle="italic", color="#555555", va="bottom"
    )

    ax_i = fig.add_subplot(gs_right[2, 1])

    plot_panels_G_H(ax_g, ax_h)
    ax_g.set_title("G   Integrated 'Triple Model' (AUROC)", loc="left", fontsize=11.5, fontweight="bold", pad=8)
    ax_h.set_title("H   Integrated 'Triple Model' (AUPRC)", loc="left", fontsize=11.5, fontweight="bold", pad=8)
    plot_panel_I(ax_i)

    # Export figure
    print(f"\nExporting Master Figure 2:\n- {output_path_png}\n- {output_path_pdf}")
    fig.savefig(output_path_png, dpi=300, bbox_inches="tight")
    fig.savefig(output_path_pdf, bbox_inches="tight")
    plt.close(fig)
    print("Execution completed successfully.\n")


if __name__ == "__main__":
    main()