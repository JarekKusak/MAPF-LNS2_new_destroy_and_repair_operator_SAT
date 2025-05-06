#!/usr/bin/env python3
"""
batch_lns.py – hromadné spouštění MAPF‑LNS, parsování logu a tvorba statistik
"""

import argparse
import os
import re
import subprocess
import sys
import time
from pathlib import Path
from itertools import product

import pandas as pd
import matplotlib.pyplot as plt

# konfigurace
MAPS            = ["ost003d.map"]#, "random-32-32-20.map"]# soubory s mapou
INSTANCES_PER_MAP = 5 # kolik scén / mapu
AGENT_COUNTS    = [100, 300]# 500, 700, 900] # -k hodnoty
MAX_ITERS       = 2 # --maxIterations
LNS_BIN         = "./lns" # cesta k binárce
RESULTS_ROOT    = Path("results") # kam ukládat

# různé heuristiky
HEURISTIC_TAG   = "adaptive" # čistě do názvu
# HEURISTIC_TAG   = "roundRobin"
# HEURISTIC_TAG   = "mostDelayedAgent"

# regexy pro parsování logu
RE_SOC_POST   = re.compile(r"\[DEBUG\] sum_of_costs po opětovném přepočtu: (\d+)")
RE_SOC_PRE    = re.compile(r"\[DEBUG\] sum_of_costs před opětovným přepočtem: (\d+)")
RE_NEIGH_OLD  = re.compile(r"\[DEBUG\] neighbor\.old_sum_of_costs .*: (\d+)")
RE_NEIGH_NEW  = re.compile(r"\[DEBUG\] neighbor\.sum_of_costs .*: (\d+)")
RE_CONFLICT   = re.compile(r"\[WARNING\] Problem after SAT:")
RE_FINAL      = re.compile(
    r"LNS\([^)]+\): runtime = ([\d\.]+), iterations = (\d+), "
    r"solution cost = (\d+), initial solution cost = (\d+), failed iterations = (\d+)"
)

def run_single(map_file: str, scen_file: str, agent_num: int, run_id: str):
    """
    Spustí ./lns s danou mapou, instancí a počtem agentů.
    Vrátí cestu k souboru logu.
    """
    out_dir = RESULTS_ROOT / run_id
    out_dir.mkdir(parents=True, exist_ok=True)
    log_path = out_dir / "log.txt"

    cmd = [
        LNS_BIN,
        "-m", map_file,
        "-a", scen_file,
        "-k", str(agent_num),
        "-o", str(out_dir / "out"), # prefix; soubory stejně nepoužijeme
        "--destoryStrategy=SAT",
        "--maxIterations", str(MAX_ITERS),
        "--screen", "0", # ticho, všechno jde do logu
        "--outputPaths", str(out_dir / "paths.txt"),
    ]

    with open(log_path, "w") as fout:
        start = time.time()
        proc = subprocess.run(cmd, stdout=fout, stderr=subprocess.STDOUT)
        elapsed = time.time() - start

    if proc.returncode != 0:
        print(f"[WARN] Běh {run_id} skončil chybou (rc={proc.returncode}).", file=sys.stderr)

    print(f"[INFO] {run_id} hotovo za {elapsed:.1f}s")
    return log_path


def parse_log(log_path: Path):
    """
    Vrátí slovník se všemi metrikami, včetně seznamu SOC po iteracích.
    """
    soc_curve      = []
    sat_improv     = []   # relativní zlepšení, pokud byl SAT použít
    sat_conflicts  = []   # True/False – vznikl konflikt?
    sat_in_action  = False
    # pomocná proměnná pro SAT session: new_cost_tmp
    new_cost_tmp   = None

    final_stats = {}

    with open(log_path, encoding="utf-8") as f:
        for line in f:
            # průběh SOC – ukládáme *pouze* hodnotu po úplném přepočtu
            m_post = RE_SOC_POST.search(line)
            if m_post:
                soc_curve.append(int(m_post.group(1)))

            # --- SAT destroy/repair session ---------------------------------
            # 1) Logger tiskne nejdřív "neighbor.sum_of_costs ..." ( = NEW ),
            #    hned poté "neighbor.old_sum_of_costs ..." ( = OLD ).
            # 2) Zlepšení počítáme teprve, když máme obě čísla
            # 3) Pokud by NEW >= OLD (nemělo by nastat), dvojici ignorujeme
            # -----------------------------------------------------------------
            m_new = RE_NEIGH_NEW.search(line)
            if m_new:
                new_cost_tmp = int(m_new.group(1))
                sat_in_action = True
                conflict_flag = False
                continue

            if sat_in_action and RE_CONFLICT.search(line):
                conflict_flag = True

            m_old = RE_NEIGH_OLD.search(line)
            if m_old and sat_in_action:
                old_cost = int(m_old.group(1))
                if old_cost > 0 and new_cost_tmp is not None and new_cost_tmp <= old_cost:
                    sat_improv.append((old_cost - new_cost_tmp) / old_cost)
                # evidence o tom, zda při seanci vznikl konflikt
                sat_conflicts.append(conflict_flag)
                # reset pro další seanci
                new_cost_tmp = None
                sat_in_action = False
                continue

            # poslední souhrnný řádek
            m_fin = RE_FINAL.search(line)
            if m_fin:
                final_stats = {
                    "runtime"          : float(m_fin.group(1)),
                    "iterations"       : int(m_fin.group(2)),
                    "final_soc"        : int(m_fin.group(3)),
                    "initial_soc"      : int(m_fin.group(4)),
                    "failed_iterations": int(m_fin.group(5)),
                }

    # průměrné relativní zlepšení mezi po sobě jdoucími iteracemi
    if len(soc_curve) >= 2:
        iter_improv = [
            (soc_curve[i - 1] - soc_curve[i]) / soc_curve[i - 1]
            for i in range(1, len(soc_curve))
            if soc_curve[i - 1] > 0
        ]
        iter_improv_mean = sum(iter_improv) / len(iter_improv) if iter_improv else 0.0
    else:
        iter_improv_mean = 0.0

    return {
        **final_stats,
        "soc_curve"       : soc_curve,
        "sat_improv_mean" : sum(sat_improv) / len(sat_improv) if sat_improv else 0.0,
        "sat_conflict_pct": 100 * sum(sat_conflicts) / len(sat_conflicts) if sat_conflicts else 0.0,
        "iter_improv_mean": iter_improv_mean,
    }

def save_curve(curve, out_png: Path, title: str):
    plt.figure()
    plt.plot(curve)
    plt.xlabel("Iterace")
    plt.ylabel("Sum of costs")
    plt.title(title)
    plt.tight_layout()
    plt.savefig(out_png)
    plt.close()

def main():
    RESULTS_ROOT.mkdir(exist_ok=True)

    all_records = []

    for map_file in MAPS:
        map_stem = Path(map_file).stem
        for inst_idx in range(1, INSTANCES_PER_MAP + 1):
            scen = f"{map_stem}-instances/{map_stem}-random-{inst_idx}.scen"

            for k in AGENT_COUNTS:
                run_tag = f"{map_stem}-i{inst_idx}-k{k}-it{MAX_ITERS}-{HEURISTIC_TAG}"
                log_path = run_single(map_file, scen, k, run_tag)
                stats = parse_log(log_path)

                # ulož průběhový graf
                save_curve(
                    stats["soc_curve"],
                    RESULTS_ROOT / run_tag / "soc.png",
                    f"{run_tag}: SOC vs. iterace",
                    )

                record = {
                    "run_id"           : run_tag,
                    "map"              : map_stem,
                    "instance"         : inst_idx,
                    "agents"           : k,
                    **{k: v for k, v in stats.items() if k != "soc_curve"},
                }
                all_records.append(record)

    df = pd.DataFrame(all_records)
    csv_path = RESULTS_ROOT / "results.csv"
    df.to_csv(csv_path, index=False)
    print(f"[OK] Výsledky zapsány do {csv_path}")

if __name__ == "__main__":
    main()