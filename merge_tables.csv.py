#!/usr/bin/env python3
import pandas as pd
from pathlib import Path

RESULTS_DIR = Path("results")          # adresář s CSV

FILES = [
    "results_all",
    "results_T1_pure",
    "results_T2_mix",
    "results_T3_runtime",
]

for stem in FILES:
    new_csv  = RESULTS_DIR / f"{stem}.csv"
    old_csv  = RESULTS_DIR / f"{stem}_backup.csv"

    if not new_csv.exists():
        print(f"[SKIP] {new_csv} neexistuje – pravděpodobně ještě nebyl vygenerován.")
        continue
    if not old_csv.exists():
        print(f"[SKIP] {old_csv} neexistuje – není co slučovat.")
        continue

    print(f"[MERGE] {old_csv.name}  +  {new_csv.name}")

    df_old = pd.read_csv(old_csv)
    df_new = pd.read_csv(new_csv)

    merged = pd.concat([df_old, df_new], ignore_index=True)

    # ───────────────────────────────────────────────
    #  Odstranění duplicit – primárně podle sloupce `run_id`,
    #  který je ve všech „hlavních“ tabulkách (kromě T3).
    #  Pro T3 chybí, tak použijeme rozumný klíč složený z parametrů běhu.
    # ───────────────────────────────────────────────
    if "run_id" in merged.columns:
        merged = merged.drop_duplicates(subset=["run_id"], keep="first")
    else:
        key = ["map", "inst", "k", "algo", "satProb",
               "destFallback", "algoFallback"]
        merged = merged.drop_duplicates(subset=key, keep="first")

    merged.to_csv(new_csv, index=False)
    print(f"   → uloženo: {new_csv}")

print("\nHotovo – sloučené tabulky jsou uložené.")