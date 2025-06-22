#!/usr/bin/env python3
import pandas as pd, seaborn as sns, matplotlib.pyplot as plt, numpy as np
from pathlib import Path
from scipy.stats import wilcoxon

sns.set_theme()
RES = Path("results")
df = pd.read_csv(RES/"results_all.csv")

# jen doběhnuté případy
df = df[df.state == "OK"].copy()

# ---------------- doplňkové sloupce ----------------
# NaN pokud SAT nevolal –- potom se ve skupinovém průměru ignoruje
df["succ_rate_run"] = df.apply(
    lambda r: r.sat_ok / r.sat_calls if r.sat_calls > 0 else np.nan,
    axis=1)

# ---------------- agregace ----------------
grp = ["map", "k", "algo", "satProb"]
agg = (
    df.groupby(grp)
      .agg(mean_runtime      = ("runtime", "mean"),
           mean_impr_pct     = ("soc_improvement_pct", "mean"),
           mean_sat_ratio    = ("sat_ratio_ops", "mean"),
           succ_rate         = ("succ_rate_run", "mean"),   # průměr jen z def. hodnot
           n_runs            = ("run_id", "count"))         # pro kontrolu počtu vzorků
      .reset_index()
)

agg.to_csv(RES / "summary.csv", index=False)

# ---------------- grafy -------------------
for m, sub in agg.groupby("map"):
    ax = sns.lineplot(data=sub, x="k", y="mean_impr_pct",
                      hue="algo", style="satProb", marker="o")
    ax.set(title=f"{m}: %-zlepšení", ylabel="%-zlepšení (SoC)")
    ax.figure.savefig(RES / f"fig_impr_{m}.png", dpi=300)
    plt.close()

sns.boxplot(data=df[df.satProb > 0],
            x="k", y="sat_ratio_ops", hue="algo")
plt.ylabel("podíl SAT v runtime (%)")
plt.savefig(RES / "fig_sat_share.png", dpi=300)
plt.close()

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

print("post-processing hotovo ✔")