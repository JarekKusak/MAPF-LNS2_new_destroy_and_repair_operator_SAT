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
    """Uloží figuru do results/ jako PNG a PDF (bez okrajů)."""
    if OUT_PNG:
        fig.savefig(RES / f"{stem}.png", dpi=300, bbox_inches="tight")
    fig.savefig(RES / f"{stem}.pdf", bbox_inches="tight")
    plt.close(fig)
# ------------------------------------------------------------------

for m, sub in agg.groupby("map"):
    ax = sns.lineplot(data=sub, x="k", y="mean_impr_pct",
                      hue="algo", style="satProb", marker="o")
    ax.set(title=f"{m}: %-zlepšení", ylabel="%-zlepšení (SoC)")
    ax.figure.savefig(RES / f"fig_impr_{m}.png", dpi=300)
    ax.figure.savefig(RES / f"fig_impr_{m}.pdf")
    plt.close()

# === a) Runtime × finální SoC =====================================
fig, ax = plt.subplots(figsize=(6, 5))
sns.kdeplot(data=df, x="runtime", y="final_soc",
            fill=True, thresh=.05, levels=60, cmap="mako", ax=ax)
sns.scatterplot(data=df, x="runtime", y="final_soc",
                hue="algo", style="satProb",
                alpha=.85, s=25, ax=ax)
ax.set_title("Runtime vs. finální SoC")
ax.set_xlabel("celkový runtime (s)")
ax.set_ylabel("finální SoC (↓ lepší)")
ax.legend(title="algo / satProb", bbox_to_anchor=(1.02, 1))
_save(fig, "fig_rt_vs_soc")

# === b) Počet agentů k × vybraná metrika ==========================
METRIC = "mean_impr_pct"          # lze snadno změnit
fig, ax = plt.subplots(figsize=(6, 4))
sns.lineplot(data=agg, x="k", y=METRIC,
             hue="algo", style="satProb",
             marker="o", ax=ax)
ax.set_ylabel(METRIC.replace("_", " "))
ax.set_title(f"{METRIC.replace('_', ' ')} v závislosti na počtu agentů k")
ax.legend(title="algo / satProb", bbox_to_anchor=(1.02, 1))
_save(fig, "fig_k_vs_metric")

# === c) Vztah podílu SAT‑runtime a zlepšení SoC ===================
# -- optional smooth LOWESS line needs statsmodels -----------------
try:
    import statsmodels.api  # noqa: F401
    _LOWESS = True
except ImportError:
    _LOWESS = False

fig, ax = plt.subplots(figsize=(6, 5))
sns.regplot(data=agg, x="mean_sat_ratio", y="mean_impr_pct",
            lowess=_LOWESS, scatter=False, color="grey", ax=ax)
sns.scatterplot(data=agg, x="mean_sat_ratio", y="mean_impr_pct",
                size="k", sizes=(20, 200),
                hue="map", style="algo", alpha=.9, ax=ax)
ax.set_xlabel("podíl SAT v runtime (%)")
ax.set_ylabel("%‑zlepšení SoC")
ax.set_title("Zlepšení SoC vs. využití SATu")
ax.legend(bbox_to_anchor=(1.02, 1), title="map / algo")
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