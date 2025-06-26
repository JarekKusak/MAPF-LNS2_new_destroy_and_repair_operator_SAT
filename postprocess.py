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

# --- nový ukazatel: zda běh s voláním SAT nic nezlepšil -------------
df["no_improvement_run"] = (
    (df["sat_ok"] > 0) & (df["soc_improvement_pct"] == 0)
)

# ---------------- agregace ----------------
# mean_runtime: průměr runtime
# mean_impr_pct: průměr %-zlepšení SoC
# mean_sat_ratio_ops: průměr podílu SAT v operacích
# mean_repairs_pct: Prům. % SAT oprav
# mean_no_improv: Podíl běhů, kdy SAT nic nezlepšil
grp = ["map", "k", "algo", "satProb"]
agg = (
    df.groupby(grp)
      .agg(mean_runtime      = ("runtime", "mean"),
           std_runtime       = ("runtime", "std"),
           mean_impr_pct     = ("soc_improvement_pct", "mean"),
           mean_sat_ratio_ops    = ("sat_ratio_ops", "mean"),
           mean_no_improv   = ("no_improvement_run", "mean"),
           mean_repairs_pct  = ("sat_success_with_repair_pct", "mean"),
           mean_sat_succ_rate = ("sat_succ_rate_run", "mean"),
           mean_iter_succ_rate = ("iter_succ_rate_run", "mean"),
           n_runs            = ("run_id", "count"))
      .reset_index()
)

agg.to_csv(RES / "results_compact.csv", index=False)   # hlavní tabulka do práce
agg.to_csv(RES / "summary.csv", index=False)           # zpětná kompatibilita

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
    if len(pp) == len(sat) and len(pp) >= 3:
        return pd.Series({
            "p": wilcoxon(pp, sat).pvalue,
            "n_pairs": len(pp)
        })
    else:
        return pd.Series({"p": np.nan, "n_pairs": len(pp)})

sig = (df.groupby(["map", "k"])
         .apply(sig_test)
         .reset_index())
# sig.to_csv(RES / "significance.csv", index=False)  # removed to avoid duplicate files
sig.to_csv(RES / "signif_pp_vs_sat.csv", index=False)

# ===========================
# Kontrola, proč chybí p-values
# ===========================
missing = sig[sig["p"].isna()]
if not missing.empty:
    print("\nKontrola map a k bez p-value (možná shodné final_soc):")
    for _, row in missing.iterrows():
        map_name, k_val = row["map"], row["k"]
        data_pp = df[(df["map"] == map_name) & (df["k"] == k_val) & (df["algo"] == "PP") & (df["satProb"] == 0)]["final_soc"]
        data_sat = df[(df["map"] == map_name) & (df["k"] == k_val) & (df["algo"] == "SAT") & (df["satProb"] == 100)]["final_soc"]
        print(f"  {map_name}, k={k_val}: PP instances={len(data_pp)}, SAT instances={len(data_sat)}. Final SOC PP: {data_pp.unique()}, SAT: {data_sat.unique()}")


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
    ax.set_title(f"{map_name}: trend {metric} vs k") 
    ax.legend(title="konfigurace", bbox_to_anchor=(1.02, 1))
    stem = f"profile_{metric}_vs_k_{map_name}"
    _save(fig, stem)

# Vytvoříme profilový trend (přes všechny mapy) pro vybranou metriku
all_maps_agg = agg.groupby(["k", "algo", "satProb"])["mean_impr_pct"].mean().reset_index()
for metric in ["mean_impr_pct"]:
    _k_vs_metric_profile(all_maps_agg, metric, "all_maps")



# ================================================================
# Agregovaný profil (průměr přes mapy) – pouze 10 bodů na ose X
# ================================================================

KS_TO_PLOT = [100, 200, 300, 400]

def _profile_mean_over_maps(k_value: int):
    """Průměrné %-zlepšení SoC přes všechny mapy.
    X-ová osa = index instance 1..10 (shodná numerace napříč mapami)."""
    sub = df[df["k"] == k_value]
    if sub.empty:
        return

    # vynechat duplicitní „SAT 0 %“
    sub = sub[~((sub["algo"] == "SAT") & (sub["satProb"] == 0))]

    # Vypočítat průměr přes mapy pro každou instanci a konfiguraci
    mean_per_inst = (
        sub.groupby(["inst", "algo", "satProb"])["soc_improvement_pct"]
           .mean()
           .reset_index()
    )

    fig, ax = plt.subplots(figsize=(6, 4))

    for (algo, p), grp in mean_per_inst.groupby(["algo", "satProb"]):
        # řadíme podle hodnoty metriky (vzestupně) – stejně jako u profilů jednotlivých map
        grp_sorted = grp.sort_values("soc_improvement_pct")
        x = np.arange(1, len(grp_sorted) + 1)
        y = grp_sorted["soc_improvement_pct"]
        label = f"{algo} / {p}%" if algo == "SAT" else "PP"
        linestyle = "--" if algo == "SAT" else "-"
        ax.plot(x, y, linestyle=linestyle, marker="o", label=label)

    ax.set_xlabel("Index instance (1-10)")
    ax.set_ylabel("%-zlepšení SoC (průměr přes mapy)")
    ax.set_title(f"Agregovaný profil všech map, k={k_value}")
    ax.legend(title="konfigurace", bbox_to_anchor=(1.02, 1))

    stem = f"profilekagg_soc_improvement_pct_allmaps_k{k_value}"
    _save(fig, stem)

for k_val in KS_TO_PLOT:
    _profile_mean_over_maps(k_val)