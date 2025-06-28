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
           mean_sat_ratio    = ("sat_ratio",      "mean"),
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
    category = stem.split('_')[0]           # prvni část slouží jako složka
    png_dir = RES / category / "png"
    pdf_dir = RES / category / "pdf"
    png_dir.mkdir(parents=True, exist_ok=True)
    pdf_dir.mkdir(parents=True, exist_ok=True)

    safe_stem = stem.replace("/", "_")      # nevkládat lomítka do názvu souboru
    if OUT_PNG:
        fig.savefig(png_dir / f"{safe_stem}.png", dpi=300, bbox_inches="tight")
    fig.savefig(pdf_dir / f"{safe_stem}.pdf", bbox_inches="tight")
    plt.close(fig)
# ------------------------------------------------------------------

# === Počet agentů k × vybraná metrika ==========================
METRIC = "mean_impr_pct"          # lze snadno změnit
fig, ax = plt.subplots(figsize=(6, 4))
sns.lineplot(data=agg, x="k", y=METRIC,
             hue="algo", style="satProb",
             marker="o", ax=ax)
ax.set_ylabel(METRIC.replace("_", " "))
ax.set_title(f"{METRIC.replace('_', ' ')} depending on the number of agents k")
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

ax.set_xlabel("SAT share in runtime (%)")
ax.set_ylabel("% SoC improvement")
ax.set_title("SoC improvement vs. SAT usage")
ax.legend(title="satProb / algo", bbox_to_anchor=(1.02, 1))

_save(fig, "fig_satshare_vs_impr")

# ================================================================
#  A) Agregovaný bar‑plot napříč všemi mapami – velikost sub‑mapy
# ================================================================
df_sat100 = df[(df["algo"] == "SAT") & (df["satProb"] == 100)]

agg_sub = (
    df_sat100
      .groupby("sub")["soc_improvement_pct"]
      .agg(mean="mean", std="std")
      .reset_index()
)

fig, ax = plt.subplots(figsize=(4, 3))
sns.barplot(
    data=agg_sub,
    x="sub",
    y="mean",
    errorbar=("sd"),
    ax=ax,
)
ax.set_xlabel("Sub-map size (side length)")
ax.set_ylabel("Average SoC improvement (%)")
ax.set_title("Influence of sub-map size – aggregation of all maps")
_save(fig, "bar_submap_allmaps")

# ================================================================
#  B) Agregovaný bar‑plot napříč všemi mapami – destroy heuristika
# ================================================================
agg_heur = (
    df_sat100
      .groupby("satHeur")["soc_improvement_pct"]
      .agg(mean="mean", std="std")
      .reset_index()
)

fig, ax = plt.subplots(figsize=(4, 3))
sns.barplot(
    data=agg_heur,
    x="satHeur",
    y="mean",
    errorbar=("sd"),
    ax=ax,
)
ax.set_xlabel("Destroy heuristic")
ax.set_ylabel("Average SoC improvement (%)")
ax.set_title("Influence of destroy heuristics – aggregation of all maps")
_save(fig, "bar_heur_allmaps")

# ================================================================
#  C) Párové Wilcoxonovy testy (9×9 vs ostatní; adaptive vs ostatní)
#     a export přehledové tabulky
# ================================================================
from scipy.stats import wilcoxon

wilcoxon_rows = []

# --- sub‑mapy ----------------------------------------------------
baseline_9 = df_sat100[df_sat100["sub"] == 9].sort_values(["map", "inst"]).reset_index(drop=True)

for sub_size in sorted(df_sat100["sub"].unique()):
    if sub_size == 9:
        continue
    test_set = df_sat100[df_sat100["sub"] == sub_size].sort_values(["map", "inst"]).reset_index(drop=True)
    if len(baseline_9) == len(test_set):
        stat, pval = wilcoxon(baseline_9["soc_improvement_pct"], test_set["soc_improvement_pct"])
        wilcoxon_rows.append({
            "category": "submap",
            "baseline": "9",
            "comparison": str(sub_size),
            "p_value": pval,
            "n_pairs": len(baseline_9)
        })

# --- heuristiky ---------------------------------------------------
baseline_adapt = df_sat100[df_sat100["satHeur"] == "adaptive"].sort_values(["map", "inst"]).reset_index(drop=True)

for heur in df_sat100["satHeur"].unique():
    if heur == "adaptive":
        continue
    test_set = df_sat100[df_sat100["satHeur"] == heur].sort_values(["map", "inst"]).reset_index(drop=True)
    if len(baseline_adapt) == len(test_set):
        stat, pval = wilcoxon(baseline_adapt["soc_improvement_pct"], test_set["soc_improvement_pct"])
        wilcoxon_rows.append({
            "category": "heuristic",
            "baseline": "adaptive",
            "comparison": heur,
            "p_value": pval,
            "n_pairs": len(baseline_adapt)
        })

wilcoxon_df = pd.DataFrame(wilcoxon_rows)
wilcoxon_df.to_csv(RES / "wilcoxon_summary.csv", index=False)

# Také uložíme agregovanou tabulku sub‑map & heuristiky
agg_sub.to_csv(RES / "submap_summary.csv", index=False)
agg_heur.to_csv(RES / "heuristic_summary.csv", index=False)

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
def _average_improvement_plot(map_name: str):
    """Generate three avg‑improvement‑vs‑k plots for *map_name*:
       1) plain means                → avg_plain_improvement_vs_k_<map>
       2) means ± SD error‑bars      → avg_sd_improvement_vs_k_<map>
       3) means ± SD + scatter points→ avg_sdpts_improvement_vs_k_<map>"""
    df_map = df[(df["map"] == map_name) & (df["state"] == "OK")]
    df_map = df_map[~((df_map["algo"] == "SAT") & (df_map["satProb"] == 0))]

    # jednotná legenda
    df_map = df_map.copy()
    df_map["config"] = df_map.apply(
        lambda r: "PP" if r.algo == "PP" else f"SAT {r.satProb}%", axis=1
    )

    def _make_plot(kind: str):
        """kind: 'plain', 'sd', 'sdpts'"""
        with_sd   = kind != "plain"
        add_pts   = kind == "sdpts"

        fig, ax = plt.subplots(figsize=(6, 4))
        sns.lineplot(
            data=df_map,
            x="k",
            y="soc_improvement_pct",
            hue="config",
            markers=True,
            estimator="mean",
            errorbar=("sd") if with_sd else None,
            err_style="bars",
            dashes=False,
            ax=ax,
        )
        if add_pts:
            # explicit scatter for the mean values
            means = (
                df_map.groupby(["config", "k"])["soc_improvement_pct"]
                      .mean()
                      .reset_index()
            )
            sns.scatterplot(
                data=means,
                x="k",
                y="soc_improvement_pct",
                hue="config",
                legend=False,
                ax=ax,
                zorder=5,
            )

        ax.set_xticks(sorted(df_map["k"].unique()))
        ax.set_xlabel("Number of agents k")
        ax.set_ylabel("Average %-improvement of SoC")
        title_extra = {"plain": "",
                       "sd":    " (± SD)",
                       "sdpts": " (± SD, points)"}[kind]
        ax.set_title(f"{map_name}: Average improvement vs k{title_extra}")
        ax.legend(title="configuration", bbox_to_anchor=(1.02, 1))

        stem = f"avg_{kind}_improvement_vs_k_{map_name}"
        _save(fig, stem)

    # create all three variants
    _make_plot("plain")
    _make_plot("sd")
    _make_plot("sdpts")

# Vytvoříme průměry pro všechny mapy
for map_name in sorted(df["map"].unique()):
    _average_improvement_plot(map_name)

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

    ax.set_xlabel("Number of agents k")
    ax.set_ylabel(PROFILE_METRICS.get(metric, {"label": metric}).get("label", metric))
    ax.set_title(f"{map_name}: trend {metric} vs k") 
    ax.legend(title="configuration", bbox_to_anchor=(1.02, 1))
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

    ax.set_xlabel("Index of instance (1-10)")
    ax.set_ylabel("%-improvement of SoC (average across maps)")
    ax.set_title(f"Aggregated profile of all maps, k={k_value}")
    ax.legend(title="configuration", bbox_to_anchor=(1.02, 1))

    stem = f"profilekagg_soc_improvement_pct_allmaps_k{k_value}"
    _save(fig, stem)

for k_val in KS_TO_PLOT:
    _profile_mean_over_maps(k_val)

# ================================================================
#  SAT‑100 % — průměrné zlepšení SoC podle velikosti sub‑mapy
#               a podle destroy‑heuristiky
# ================================================================

def _sat100_by_submap(df_map: pd.DataFrame, map_name: str):
    """Bar-plot: sub × mean_impr_pct   (hue = k)."""
    df_sat100 = df_map[(df_map["algo"] == "SAT") & (df_map["satProb"] == 100)]
    if df_sat100.empty or "sub" not in df_sat100.columns:
        return
    # agregace
    data = (df_sat100
            .groupby(["sub", "k"])["soc_improvement_pct"]
            .mean()
            .reset_index()
            .rename(columns={"soc_improvement_pct": "mean_impr_pct"}))

    fig, ax = plt.subplots(figsize=(6, 4))
    # Sloupcový graf: X = velikost sub‑mapy, Y = průměrné zlepšení, hue = k
    sns.barplot(
        data=data,
        x="sub",
        y="mean_impr_pct",
        hue="k",
        errorbar=("ci", 95),
        ax=ax,
    )
    ax.set_xlabel("Sub-map size (side length)")
    ax.set_ylabel("%-improvement of SoC (SAT 100%)")
    ax.set_title(f"{map_name}: SAT 100 % – influence of sub-map size")
    ax.legend(title="Number of agents k", bbox_to_anchor=(1.02, 1))
    _save(fig, f"bar_impr_sat_submap_{map_name}")

def _sat100_by_heur(df_map: pd.DataFrame, map_name: str):
    """Bar-plot: satHeur × mean_impr_pct   (hue = k)."""
    df_sat100 = df_map[(df_map["algo"] == "SAT") & (df_map["satProb"] == 100)]
    if df_sat100.empty or "satHeur" not in df_sat100.columns:
        return
    data = (df_sat100
            .groupby(["satHeur", "k"])["soc_improvement_pct"]
            .mean()
            .reset_index()
            .rename(columns={"soc_improvement_pct": "mean_impr_pct"}))

    fig, ax = plt.subplots(figsize=(6, 4))
    # Sloupcový graf: X = heuristika, Y = průměrné zlepšení, hue = k
    sns.barplot(
        data=data,
        x="satHeur",
        y="mean_impr_pct",
        hue="k",
        errorbar=("ci", 95),
        ax=ax,
    )
    ax.set_xlabel("Destroy heuristic")
    ax.set_ylabel("%-improvement of SoC (SAT 100%)")
    ax.set_title(f"{map_name}: SAT 100 % – influence of heuristics")
    ax.legend(title="Number of agents k", bbox_to_anchor=(1.02, 1))
    _save(fig, f"bar_impr_sat_heur_{map_name}")

# --- generování grafů pro každou mapu ---------------------------
for map_name, sub_df in df.groupby("map"):
    _sat100_by_submap(sub_df, map_name)
    _sat100_by_heur(sub_df, map_name)

# ================================================================
#  Export performance table (mean_impr_pct per satProb)
# ================================================================
print("Exporting summary table …")
tables_dir = RES / "tables"      # store inside results/
tables_dir.mkdir(exist_ok=True)

cols = ["map", "k", "satProb", "mean_impr_pct", "std_runtime"]
perf_df = (
    pd.read_csv(RES / "summary.csv")[cols]
      .pivot_table(index=["map", "k"],
                   columns="satProb",
                   values="mean_impr_pct")
      .round(2)
)

try:
    # pandas >=2.2 uses Styler which needs jinja2; basic export avoids it
    perf_df.to_latex(tables_dir / "sat_performance.tex",
                     index=True,
                     na_rep="",
                     # if available in your pandas:  method="basic"
                     )
    print("LaTeX table written to tables/sat_performance.tex")
except ImportError:
    # Fallback: write CSV instead and notify the user
    perf_df.to_csv(tables_dir / "sat_performance.csv")
    print("Jinja2 missing → wrote tables/sat_performance.csv instead.")
