#!/usr/bin/env python3
"""
stats.py – hromadné spouštění MAPF-LNS, parsování logu a tvorba statistik
"""

import os
import re
import subprocess
import sys
import time
from datetime import datetime
from pathlib import Path
from itertools import product

import pandas as pd
import matplotlib.pyplot as plt

# ---------------------------------------------------------------------------
# EXPERIMENT MATRIX – adjust here to generate all desired runs
# ---------------------------------------------------------------------------
PURE_REPLANS   = ["PP", "CBS"]          # pro T1  (100 % PP / 100 % CBS)
PURE_SAT       = True                   # pokud True, vygeneruje i čistých 100 % SAT
MIX_PROBS      = [100, 50, 20, 0]       # pro T2

MAPS              = ["ost003d", "random-32-32-20"]
                     #"lt_gallowstemplar_n", "room-32-32-4", "Paris_1_256"]      # map “stems”
INSTANCES_PER_MAP = 2                                          # scen files 1..N
AGENT_COUNTS      = [100, 500]#, 500]                              # -k values
# wall‑clock time budgets (seconds) – LNS receives identical cut‑off time
CUTOFFS           = [2]# 30]                                  # --cutoffTime (-t)
MAX_ITERS         = [10000]  # high default; customise if you want true iteration caps
SAT_PROBS         = [100, 0]#, 50, 0]#, 20, 0]                             # --satProb
SAT_HEURISTICS    = ["adaptive"]#, "roundRobin", "mostDelayed"]    # --satHeuristic
SUBMAP_SIDES      = [3]#, 5]
LNS_BIN           = "./lns"
RESULTS_ROOT      = Path("results")

# -------- regexes matching [STAT] lines in log -----------------------------
RE_SOC_POST   = re.compile(r"\[STAT\] sum_of_costs after recomputation: (\d+)")
RE_NEIGH_NEW  = re.compile(r"\[STAT\] neighbor\.sum_of_costs before recomputation: (\d+)")
RE_NEIGH_OLD  = re.compile(r"\[STAT\] neighbor\.old_sum_of_costs before recomputation: (\d+)")
RE_SAT_RT     = re.compile(r"\[STAT\] SAT total runtime = ([\d\.eE+-]+) s")
RE_OTH_RT     = re.compile(r"\[STAT\] Other operators runtime = ([\d\.eE+-]+) s")
RE_CONFLICT   = re.compile(r"\[WARNING\] Problem after SAT:")
RE_FINAL      = re.compile(
    r"\[STAT\] .*: runtime = ([\d\.eE+-]+), iterations = (\d+), "
    r"solution cost = (\d+), initial solution cost = (\d+), failed iterations = (\d+)"
)

def run_single(map_file: str, scen_file: str, agent_num: int, run_id: str,
               *, heur_tag: str, submap: int, sat_prob: int, cutoff: int, max_iter: int):
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
        "-t", str(cutoff),
        "--destoryStrategy=SAT",
        f"--satHeuristic={heur_tag}",
        f"--satSubmap={submap}",
        f"--satProb={sat_prob}",
        "--maxIterations", str(max_iter),
        "--satDebug=0",
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


def parse_log(log_path: Path, external_init_soc: int | None = None):
    """
    Vrátí slovník se všemi metrikami …  external_init_soc může přepsat počáteční SOC.
    """
    soc_curve      = []
    sat_improv     = []   # relativní zlepšení, pokud byl SAT použít
    sat_conflicts  = []   # True/False – vznikl konflikt?
    sat_in_action  = False
    # pomocná proměnná pro SAT session: new_cost_tmp
    new_cost_tmp   = None
    conflict_pending = False  # for warnings that appear after session is closed

    final_stats = {}
    # --- první výskyt hodnot pro výpočet počátečního SOC -----------------
    soc_pre_first: int | None   = None   # [DEBUG] sum_of_costs před opětovným přepočtem
    neigh_new_first: int | None = None   # [DEBUG] neighbor.sum_of_costs před opětovným přepočtem
    neigh_old_first: int | None = None   # [DEBUG] neighbor.old_sum_of_costs před opětovným přepočtem

    with open(log_path, encoding="utf-8") as f:
        for line in f:
            # --- pick up SAT and other operator runtimes -----------------
            m_sat = RE_SAT_RT.search(line)
            if m_sat:
                final_stats["sat_runtime"] = float(m_sat.group(1))
                continue  # nothing else of interest on this line

            m_oth = RE_OTH_RT.search(line)
            if m_oth:
                final_stats["other_runtime"] = float(m_oth.group(1))
                continue
            # --------------------------------------------------------------

            # průběh SOC – ukládáme *pouze* hodnotu po úplném přepočtu
            m_post = RE_SOC_POST.search(line)
            if m_post:
                soc_curve.append(int(m_post.group(1)))

            # -------------------------------------------------------------
            # Zachyť první výskyt hodnot před opětovným přepočtem
            # -------------------------------------------------------------
            if soc_pre_first is None:
                m_pre = RE_SOC_POST.search(line)
                if m_pre:
                    soc_pre_first = int(m_pre.group(1))

            if neigh_new_first is None:
                m_nn = RE_NEIGH_NEW.search(line)
                if m_nn:
                    neigh_new_first = int(m_nn.group(1))

            if neigh_old_first is None:
                m_no = RE_NEIGH_OLD.search(line)
                if m_no:
                    neigh_old_first = int(m_no.group(1))

            # --- SAT destroy/repair session ---------------------------------
            # 1) Logger tiskne nejdřív "neighbor.sum_of_costs ..."  (NEW)
            #    => začátek seance (sat_in_action = True)
            # 2) Pak "neighbor.old_sum_of_costs ..."                (OLD)
            #    => konec seance (sat_in_action = False, výsledek uložíme)
            # 3) Upozornění "[WARNING] Problem after SAT:" může přijít
            #    až *po* řádku OLD – musíme tedy případně zpětně přepsat
            #    poslední položku v sat_conflicts.
            # -----------------------------------------------------------------
            m_new = RE_NEIGH_NEW.search(line)
            if m_new:
                new_cost_tmp   = int(m_new.group(1))
                sat_in_action  = True
                conflict_flag  = False       # reset pro aktuální seanci
                continue

            # Zachytíme případný konflikt
            if RE_CONFLICT.search(line):
                if sat_in_action:
                    conflict_flag = True     # ještě jsme ve stejné seanci
                else:
                    # seance už skončila – přepiš poslední zaznamenaný záznam
                    if sat_conflicts:
                        sat_conflicts[-1] = True
                # pokračujeme – warning neobsahuje jiné užitečné regexy
                continue

            m_old = RE_NEIGH_OLD.search(line)
            if m_old and sat_in_action:
                old_cost = int(m_old.group(1))
                if old_cost > 0 and new_cost_tmp is not None and new_cost_tmp <= old_cost:
                    sat_improv.append((old_cost - new_cost_tmp) / old_cost)
                sat_conflicts.append(conflict_flag)
                # reset
                new_cost_tmp  = None
                sat_in_action = False
                continue

            # poslední souhrnný řádek
            m_fin = RE_FINAL.search(line)
            if m_fin:
                final_stats.update({
                    "runtime"          : float(m_fin.group(1)),
                    "iterations"       : int(m_fin.group(2)),
                    "final_soc"        : int(m_fin.group(3)),
                    "initial_soc"      : int(m_fin.group(4)),
                    "failed_iterations": int(m_fin.group(5)),
                })

    # --- initial solution cost ------------------------------------------
    # Přednostně vezmeme hodnotu ze souboru out‑LNS.csv (sloupec
    # „initial solution cost“), pokud je k dispozici.
    if external_init_soc is not None:
        final_stats["initial_soc"] = external_init_soc

    # do křivky přidáme počáteční náklady, aby se počítalo zlepšení i v 1. iteraci
    if final_stats and "initial_soc" in final_stats:
        full_curve = [final_stats["initial_soc"]] + soc_curve
    else:
        full_curve = soc_curve[:]

    if len(full_curve) >= 2:
        iter_improv = [
            (full_curve[i - 1] - full_curve[i]) / full_curve[i - 1]
            for i in range(1, len(full_curve))
            if full_curve[i - 1] > 0
        ]
        iter_improv_mean = sum(iter_improv) / len(iter_improv) if iter_improv else 0.0
    else:
        iter_improv_mean = 0.0

    return {
        **final_stats,
        "soc_curve"       : soc_curve,
        "sat_improv_mean" : sum(sat_improv) / len(sat_improv) if sat_improv else 0.0,
        "iter_improv_mean": iter_improv_mean,
        "sat_conflict_pct": 100 * sum(sat_conflicts) / len(sat_conflicts) if sat_conflicts else 0.0,
        "sat_runtime"   : final_stats.get("sat_runtime", 0.0),
        "other_runtime" : final_stats.get("other_runtime", 0.0),
        "sat_ratio"     : (100.0 * final_stats.get("sat_runtime", 0.0) /
                           max(1e-9, final_stats.get("sat_runtime",0.0)+final_stats.get("other_runtime",0.0))),
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

    for m in MAPS:
        map_path  = f"maps/{m}.map"
        for inst in range(1, INSTANCES_PER_MAP + 1):
            scen_path = f"instances/{m}-instances/{m}-random-{inst}.scen"
            for k, prob, heur, sub, T, iters in product(AGENT_COUNTS, SAT_PROBS, SAT_HEURISTICS, SUBMAP_SIDES, CUTOFFS, MAX_ITERS):
                tag = f"{m}-i{inst}-k{k}-t{T}-it{iters}-{heur}-sub{sub}-p{prob}"
                log_path = run_single(map_path, scen_path, k, tag,
                                      heur_tag=heur, submap=sub, sat_prob=prob,
                                      cutoff=T, max_iter=iters)
                out_dir = log_path.parent
                init_soc = None
                lns_csv = out_dir / "out-LNS.csv"
                if lns_csv.exists():
                    try:
                        df0 = pd.read_csv(lns_csv)
                        if not df0.empty and "initial solution cost" in df0.columns:
                            init_soc = int(df0["initial solution cost"].iloc[0])
                    except Exception as e:
                        print(f"[WARN] Cannot read {lns_csv}: {e}", file=sys.stderr)

                stats = parse_log(log_path, external_init_soc=init_soc)
                # save SOC curve
                save_curve(stats["soc_curve"],
                           out_dir / "soc.png",
                           f"{tag}: SOC vs. iter")
                rec = {"run_id": tag, "map": m, "instance": inst, "agents": k,
                       "sat_prob": prob, "heuristic": heur, "submap": sub,
                       **{k: v for k, v in stats.items() if k != "soc_curve"}}
                all_records.append(rec)

    df = pd.DataFrame(all_records)
    csv_path = RESULTS_ROOT / "results.csv"
    df.to_csv(csv_path, index=False)
    print(f"[OK] Výsledky zapsány do {csv_path}")

    ts_name = RESULTS_ROOT / f"results_{datetime.now():%Y%m%d-%H%M%S}.csv"
    df.to_csv(ts_name, index=False)
    print(f"[OK] Backup written to {ts_name}")

if __name__ == "__main__":
    main()