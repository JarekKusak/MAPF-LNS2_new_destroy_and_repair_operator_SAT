import pandas as pd, pathlib, re, shutil
RES = pathlib.Path("results")

bad_cond = lambda df: (
        (df["map"].isin(["random-32-32-20", "room-64-64-16"])) &
        (df["k"] == 400) &
        (df["algo"] == "SAT") &
        (df["satProb"] == 100)
)

for fname in ["results_all-backup.csv", "results_T1_pure-backup.csv"]:
    df = pd.read_csv(RES / fname)
    df_clean = df[~bad_cond(df)]
    df_clean.to_csv(RES / fname, index=False)      # přepíšu zálohu očištěnou verzí