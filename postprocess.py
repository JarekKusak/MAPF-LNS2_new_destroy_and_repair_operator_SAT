#!/usr/bin/env python3
import pandas as pd, seaborn as sns, matplotlib.pyplot as plt, numpy as np
from pathlib import Path
from scipy.stats import wilcoxon

sns.set_theme()
RES = Path("results")
df = pd.read_csv(RES/"results_all.csv")

for col in ("sat_ok", "sat_repairs"):
    if col not in df.columns:
        df[col] = 0

# jen doběhnuté případy
df = df[df.state == "OK"].copy()

# ---------------- doplňkové sloupce ----------------
# NaN pokud SAT nevolal –- potom se ve skupinovém průměru ignoruje
df["sat_succ_rate_run"] = df.apply(
    lambda r: r.sat_ok / r.sat_calls if r.sat_calls > 0 else np.nan,
    axis=1)

df["iter_succ_rate_run"] = df.apply(
    lambda r: r.succesful_iterations / r.outer_iterations
    if r.outer_iterations else np.nan,
    axis=1)

df["sat_success_with_repair_pct"] = df.apply(
    lambda r: 100.0 * r.sat_repairs / r.sat_ok if r.sat_ok else np.nan,
    axis=1)

# ---------------- agregace ----------------
grp = ["map", "k", "algo", "satProb"]
agg = (
    df.groupby(grp)
      .agg(mean_runtime      = ("runtime", "mean"),
           std_runtime       = ("runtime", "std"),
           mean_impr_pct     = ("soc_improvement_pct", "mean"),
           mean_sat_ratio    = ("sat_ratio_ops", "mean"),
           mean_repairs_pct  = ("sat_success_with_repair_pct", "mean"),
           mean_sat_succ_rate = ("sat_succ_rate_run", "mean"),
           mean_iter_succ_rate = ("iter_succ_rate_run", "mean"),
           n_runs            = ("run_id", "count"))
      .reset_index()
)

agg.to_csv(RES / "summary.csv", index=False)

# ---------------- grafy -------------------
# ------------------------------------------------------------------
# Pomocná funkce pro ukládání obrázků zároveň do PNG i PDF ----------
OUT_PNG = True        # přepni na False, když chceš jen PDF

def _save(fig: plt.Figure, stem: str):
    """Uloží figuru do results/{category}/png a results/{category}/pdf (bez okrajů)."""
    category = stem.split('_')[0]
    png_dir = RES / category / "png"
    pdf_dir = RES / category / "pdf"
    png_dir.mkdir(parents=True, exist_ok=True)
    pdf_dir.mkdir(parents=True, exist_ok=True)

    if OUT_PNG:
        fig.savefig(png_dir / f"{stem}.png", dpi=300, bbox_inches="tight")
    fig.savefig(pdf_dir / f"{stem}.pdf", bbox_inches="tight")
    plt.close(fig)
# ------------------------------------------------------------------

# === Počet agentů k × vybraná metrika ==========================
METRIC = "mean_impr_pct"          # lze snadno změnit
fig, ax = plt.subplots(figsize=(6, 4))
sns.lineplot(data=agg, x="k", y=METRIC,
             hue="algo", style="satProb",
             marker="o", ax=ax)
ax.set_ylabel(METRIC.replace("_", " "))
ax.set_title(f"{METRIC.replace('_', ' ')} v závislosti na počtu agentů k")
ax.legend(title="algo / satProb", bbox_to_anchor=(1.02, 1))
_save(fig, "fig_k_vs_metric")

# === Relace podílu SAT‑runtime a zlepšení SoC =========================
fig, ax = plt.subplots(figsize=(6, 4))

# Globální trend (obyčejná regrese bez LOWESS, ať není závislost na statsmodels)
sns.regplot(
    data=agg,
    x="mean_sat_ratio",
    y="mean_impr_pct",
    scatter=False,
    color="grey",
    line_kws={"linewidth": 1},
    ax=ax,
)

# Body: barvy → satProb, tvary → algo, velikost → k
sns.scatterplot(
    data=agg,
    x="mean_sat_ratio",
    y="mean_impr_pct",
    hue="satProb",
    palette="rocket_r",
    style="algo",
    size="k",
    sizes=(40, 200),
    alpha=.9,
    ax=ax,
)

ax.set_xlabel("podíl SAT v runtime (%)")
ax.set_ylabel("% zlepšení SoC")
ax.set_title("Zlepšení SoC vs. využití SATu")
ax.legend(title="satProb / algo", bbox_to_anchor=(1.02, 1))
_save(fig, "fig_satshare_vs_impr")

# ---------------- významnost PP vs. čistý SAT -------------
def sig_test(sub):
    pp  = sub[(sub.algo == "PP")  & (sub.satProb == 0)]["final_soc"]
    sat = sub[(sub.algo == "SAT") & (sub.satProb == 100)]["final_soc"]
    # Wilcoxon požaduje shodný počet pozorování –- jinak vrátíme NaN
    return wilcoxon(pp, sat).pvalue if len(pp) == len(sat) and len(pp) >= 3 else np.nan

sig = (df.groupby(["map", "k"])
         .apply(sig_test)
         .rename("p")
         .reset_index())
sig.to_csv(RES / "significance.csv", index=False)


print("Post-processing dokončeno ✔")

# ------------------------------------------------------------------
# ============ Instance‑profile křivky pro vybrané metriky =========
#
# Pro každou kombinaci (map, k) vykreslíme křivky instancí
#   – každá křivka = konkrétní konfigurace (algo + satProb)
#   – osa X = instance seřazené podle hodnoty metriky
#   – osa Y = hodnota metriky
#
# Přidáváme čtyři metriky (sloupce se vytvoří, pokud chybí):
#   1) delta_soc_run            = initial_soc - final_soc
#   2) soc_improvement_pct      = už existuje
#   3) sat_ratio_ops            = už existuje
#   4) sat_repairs              = už existuje
# ------------------------------------------------------------------
if "delta_soc_run" not in df.columns:
    df["delta_soc_run"] = df["initial_soc"] - df["final_soc"]

PROFILE_METRICS = {
    "delta_soc_run":       {"label": "Δ SoC (initial - final)"},
    "soc_improvement_pct": {"label": "%-zlepšení SoC"}
}

def _profile_plot(sub_df: pd.DataFrame, metric: str, map_name: str, k_val: int):
    fig, ax = plt.subplots(figsize=(6, 4))

    # Vykreslíme křivku pro každou (algo, satProb) konfiguraci,
    # kromě „SAT 0 %“ (je to duplicitní k PP).
    filtered = sub_df[~((sub_df["algo"] == "SAT") & (sub_df["satProb"] == 0))]
    grouped = (
        filtered.groupby(["algo", "satProb"])
        .apply(lambda g: g.sort_values(metric))
        .reset_index(drop=True)
        .groupby(["algo", "satProb"])
    )
    for (algo, p), grp in grouped:
        y = grp[metric].values
        x = np.arange(1, len(y) + 1)
        label = f"{algo} / {p}%" if algo == "SAT" else "PP"
        linestyle = "--" if algo == "SAT" else "-"
        ax.plot(x, y, linestyle=linestyle, marker="", label=label)

    ax.set_xlabel("Instance (seřazeno)")
    ax.set_ylabel(PROFILE_METRICS[metric]["label"])
    ax.set_title(f"{map_name}, k={k_val}: profil {PROFILE_METRICS[metric]['label']}")
    ax.legend(title="konfigurace", bbox_to_anchor=(1.02, 1))
    stem = f"profile_{metric}_{map_name}_k{k_val}"
    _save(fig, stem)

# vytvoříme profily pro všechny mapy a k
for metric in PROFILE_METRICS:
    for (m, k_val), sub in df.groupby(["map", "k"]):
        _profile_plot(sub, metric, m, k_val)

# =========================
# Agregovaný průměr soc_improvement_pct vs k
# =========================
def _average_improvement_plot(df: pd.DataFrame, map_name: str):
    """Vykreslí průměr soc_improvement_pct pro všechny instance a konfigurace
    na dané mapě v závislosti na k."""
    fig, ax = plt.subplots(figsize=(6, 4))
    for (algo, sat_prob), grp in df.groupby(["algo", "satProb"]):
        ax.plot(grp["k"], grp["mean_impr_pct"],
                marker="o",
                linestyle="--" if algo == "SAT" else "-",
                label=f"{algo} / {sat_prob}%" if algo == "SAT" else "PP")

    ax.set_xlabel("Počet agentů k")
    ax.set_ylabel("Průměrné zlepšení SoC (%)")
    ax.set_title(f"{map_name}: Průměrné zlepšení SoC vs k")
    ax.legend(title="konfigurace", bbox_to_anchor=(1.02, 1))
    stem = f"avg_improvement_vs_k_{map_name}"
    _save(fig, stem)

# Vytvoříme průměry pro všechny mapy
for map_name, sub in agg.groupby("map"):
    _average_improvement_plot(sub, map_name)

# =========================
# Profilový styl trend mean_impr_pct vs k
# =========================
def _k_vs_metric_profile(df: pd.DataFrame, metric: str, map_name: str = "all_maps"):
    """Vykreslí trend mean_impr_pct vs k jako profile-style graf bez CI."""
    fig, ax = plt.subplots(figsize=(6, 4))
    grouped = df.groupby(["algo", "satProb"])
    for (algo, p), grp in grouped:
        grp_sorted = grp.sort_values("k")  # Seřadí k
        y = grp_sorted[metric].values
        x = grp_sorted["k"].values
        label = f"{algo} / {p}%" if algo == "SAT" else "PP"
        linestyle = "--" if algo == "SAT" else "-"
        ax.plot(x, y, linestyle=linestyle, marker="o", label=label)

    ax.set_xlabel("Počet agentů k")
    ax.set_ylabel(PROFILE_METRICS.get(metric, {"label": metric}).get("label", metric))
    ax.set_title(f"{map_name}: trend {metric} vs k (profilový styl)") 
    ax.legend(title="konfigurace", bbox_to_anchor=(1.02, 1))
    stem = f"profile_{metric}_vs_k_{map_name}"
    _save(fig, stem)

# Vytvoříme profilový trend (přes všechny mapy) pro vybranou metriku
all_maps_agg = agg.groupby(["k", "algo", "satProb"])["mean_impr_pct"].mean().reset_index()
for metric in ["mean_impr_pct"]:
    _k_vs_metric_profile(all_maps_agg, metric, "all_maps")