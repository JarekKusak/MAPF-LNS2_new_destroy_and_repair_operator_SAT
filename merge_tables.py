import pandas as pd, pathlib
RES = pathlib.Path("results")

def merge(new, old, out):
    df_new = pd.read_csv(RES / new)
    df_old = pd.read_csv(RES / old)
    df = pd.concat([df_old, df_new], ignore_index=True)
    df.to_csv(RES / out, index=False)

merge("results_all.csv",      "results_all-backup.csv",      "results_all.csv")
merge("results_T1_pure.csv",  "results_T1_pure-backup.csv",  "results_T1_pure.csv")