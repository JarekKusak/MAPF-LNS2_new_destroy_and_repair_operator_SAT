import shutil
import pathlib as pl
import pandas as pd

# ────────────────────────────────────────────────────────────────
# 1) SMAZÁNÍ adresářů/souborů „It_*“
#    – projdeme několik kořenů (můžeš upravit podle potřeby)
# ────────────────────────────────────────────────────────────────
ROOTS_TO_CLEAN = [
    pl.Path("maps"),                       # mapové soubory
    pl.Path("instances"),                  # adresáře s .scen
    pl.Path("results"),                    # případné staré výsledky
]

for root in ROOTS_TO_CLEAN:
    if not root.exists():
        continue
    for item in root.glob("It_*"):         # jen názvy začínající „It_“
        if item.is_dir():
            print(f"[CLEAN] Removing directory {item}")
            shutil.rmtree(item)
        else:
            print(f"[CLEAN] Removing file      {item}")
            item.unlink()

# ────────────────────────────────────────────────────────────────
# 2) VYČIŠTĚNÍ CSV TABULEK (beze změny oproti tvému kódu)
# ────────────────────────────────────────────────────────────────
for name in (
    "results_all.csv",
    "results_T1_pure.csv",
    "results_T2_mix.csv",
    "results_T3_runtime.csv",
):
    p = pl.Path("results") / name
    if p.exists():
        df = pd.read_csv(p)
        # cut rows s mapou začínající „It_“
        df = df[~df["map"].str.startswith("It_")]
        df.to_csv(p, index=False)
        print(f"[CSV] Cleaned {name}")