#!/usr/bin/env python3
"""
tables.py
─────────
Enumerates all benchmark runs for the bachelor-thesis experiments.
Launches `lns` with the chosen parameter matrix.
Parses each run’s log to collect the headline statistics required for
  – T1 : pure replanning strategies (PP / CBS / optional pure SAT).
  – T2 : PP + SAT mixes (varying SAT probability, fallback strategies).
  – T3 : runtime statistics (SAT time, other LNS ops time, ratios).
Produces three CSV files in <repo>/results/:
      ├─ results_all.csv   – every single run
      ├─ results_T1_pure.csv
      ├─ results_T2_mix.csv
      └─ results_T3_runtime.csv
"""

from itertools import product
from pathlib import Path
import re
import subprocess
import time
import sys
import pandas as pd
import matplotlib.pyplot as plt
import shutil

# CONFIGURATION MATRIX
MAPS              = {"ost003d"}
INSTANCES_PER_MAP = 1
AGENT_COUNTS      = [100]
TIMEOUTS          = [5] # 5 (+ margin)
MAX_ITERS         = [5_000] # arbitrary long
SUBMAP_SIDES      = [3]

PURE_REPLANS      = ["PP"]
INCLUDE_PURE_SAT  = True
MIX_PROBS         = [100, 50, 20] 
SAT_HEURISTICS    = ["adaptive"]

FALLBACK_DESTS  = ["Random", "Intersection"]
FALLBACK_ALGOS  = ["PP","CBS"]#, "EECBS"]

SAFE_MARGIN = 2  # seconds added on top of cfg['T'] to forcibly kill hanging runs

'''
MAPS = {"random-32-32-20", "room-32-32-4", "warehouse-10-20-10-2-1",
        "maze-32-32-4", "Paris_1_256"}
INSTANCES_PER_MAP = 5
AGENT_COUNTS      = [100, 300, 500]
TIMEOUTS          = [30]
SUBMAP_SIDES      = [5]
MIX_PROBS         = [50, 20] # 0 a 100 jsou generovány mimo MIX
SAT_HEURISTICS    = ["adaptive"]
FALLBACK_DESTS    = ["Random", "Intersection"]
FALLBACK_ALGOS    = ["PP", "CBS", "EECBS"]
'''

LNS_BIN     = "./lns"          # path to compiled solver
RESULTS_DIR = Path("results").absolute()
RESULTS_DIR.mkdir(exist_ok=True)

# REGEX PATTERNS (must match solver output)

RE_FINAL   = re.compile(
    r"\[STAT\] .*: runtime = ([\d\.eE+-]+), iterations = (\d+), "
    r"solution cost = (\d+), initial solution cost = (\d+), failed iterations = (\d+)"
)
RE_SAT_RT  = re.compile(r"\[STAT\] SAT total runtime = ([\d\.eE+-]+) s")
RE_OTH_RT  = re.compile(r"\[STAT\] Other operators runtime = ([\d\.eE+-]+) s")
RE_SOC_POST = re.compile(r"\[STAT\] sum_of_costs after recomputation: (\d+)")

# CASE GENERATION

cases: list[dict] = []

# T1 ─ pure strategies (PP / CBS / optional pure SAT)
for m, scen_i, k, T, iters, sub in product(
        MAPS, range(1, INSTANCES_PER_MAP + 1),
        AGENT_COUNTS, TIMEOUTS, MAX_ITERS, SUBMAP_SIDES):
    for repl in PURE_REPLANS:
        cases.append(dict(kind="PURE", algo=repl, satProb=0,
                          dest="Intersection",   # any non-SAT strategy
                          map=m, inst=scen_i, k=k, T=T, iters=iters, sub=sub))
    if INCLUDE_PURE_SAT:
        cases.append(dict(kind="PURE", algo="PP", satProb=100,
                          dest="SAT", satHeur="adaptive",
                          map=m, inst=scen_i, k=k, T=T, iters=iters, sub=sub))

# T2 ─ PP + SAT mixes
for m, scen_i, k, T, iters, sub, prob, heur, fb_dest, fb_algo in product(
        MAPS, range(1, INSTANCES_PER_MAP + 1), AGENT_COUNTS,
        TIMEOUTS, MAX_ITERS, SUBMAP_SIDES,
        MIX_PROBS, SAT_HEURISTICS,
        FALLBACK_DESTS, FALLBACK_ALGOS):
    # skip degenerate mixes – 0 % >>> žádný SAT, 100 % >>> čistý SAT už máme v "PURE"
    if prob in (0, 100):
        continue
    cases.append(dict(kind="MIX",
                      algo="PP", # primary non‑SAT replanner (internal)
                      satProb=prob,
                      dest="SAT",
                      satHeur=heur,
                      destFallback=fb_dest,
                      algoFallback=fb_algo,
                      map=m, inst=scen_i, k=k, T=T, iters=iters, sub=sub))

# LAUNCH / PARSE HELPERS

def build_cmd(cfg: dict, out_dir: Path) -> list[str]:
    """Build command-line for a single solver run."""
    cmd = [
        LNS_BIN,
        "-m", f"maps/{cfg['map']}.map",
        "-a", f"instances/{cfg['map']}-instances/{cfg['map']}-random-{cfg['inst']}.scen",
        "-k", str(cfg["k"]),
        "-o", str(out_dir / "out"),
        "-t", str(cfg["T"]),
        "--maxIterations", str(cfg["iters"]),
        "--screen", "0",
        "--outputPaths", str(out_dir / "paths.txt"),
        f"--replanAlgo={cfg['algo']}",
        f"--destoryStrategy={cfg['dest']}",
        f"--satSubmap={cfg['sub']}",
        f"--satProb={cfg['satProb']}",
        "--satDebug=0", # keep logs small
    ]
    # Add fallback strategy parameters when the primary destroy strategy is SAT.
    if cfg['dest'] == 'SAT':
        cmd.append(f"--destoryStrategyFallback={cfg.get('destFallback', 'Intersection')}")
        cmd.append(f"--replanAlgoFallback={cfg.get('algoFallback', 'PP')}")
    if cfg["dest"] == "SAT":
        cmd.append(f"--satHeuristic={cfg['satHeur']}")
    return cmd

def parse_log(log_path: Path) -> tuple[dict, list[int]]:
    """
    Extract headline statistics from the solver log.

    Returns
    -------
    (stats, curve)
        stats … dict with aggregated numbers (no SoC curve)
        curve … list with SOC-after-recomputation values –
                caller may save it to a separate file for plotting
    """
    stats: dict[str, float | int] = {}
    curve: list[int] = []

    with open(log_path, encoding="utf-8") as f:
        for line in f:
            if m := RE_FINAL.search(line):
                stats.update(runtime=float(m[1]),
                             iterations=int(m[2]),
                             final_soc=int(m[3]),
                             initial_soc=int(m[4]),
                             failed_iterations=int(m[5]))
            elif m := RE_SAT_RT.search(line):
                stats["sat_runtime"] = float(m[1])
            elif m := RE_OTH_RT.search(line):
                stats["other_runtime"] = float(m[1])
            elif m := RE_SOC_POST.search(line):
                curve.append(int(m[1]))

    # defaults
    stats.setdefault("sat_runtime", 0.0)
    stats.setdefault("other_runtime", 0.0)
    stats.setdefault("runtime", 0.0)   # total wall‑clock time of the run

    # Derived fields:
    # (1) Share of SAT time *within all operator time* (SAT + “other” LNS ops)
    total_ops = stats["sat_runtime"] + stats["other_runtime"]
    stats["sat_ratio_ops"] = 100.0 * stats["sat_runtime"] / total_ops if total_ops else 0.0

    # (2) Share of SAT time w.r.t. the entire wall‑time reported by the solver
    stats["sat_ratio"] = (
        100.0 * stats["sat_runtime"] / stats["runtime"]
        if stats["runtime"] else 0.0
    )

    # (3) Percentage improvement of final SoC vs. initial SoC (positive = better)
    if "initial_soc" in stats and stats["initial_soc"]:
        stats["soc_improvement_pct"] = 100.0 * (stats["initial_soc"] - stats["final_soc"]) / stats["initial_soc"]
    else:
        stats["soc_improvement_pct"] = 0.0

    return stats, curve

# MAIN LOOP

records = []
for cfg in cases:
    tag = (f"{cfg['map']}-i{cfg['inst']}-k{cfg['k']}-t{cfg['T']}"
           f"-sub{cfg['sub']}-p{cfg['satProb']}"
           f"-{cfg.get('satHeur','')}")
    # add fallback info only for MIX cases (keys present)
    if cfg.get("destFallback"):
        tag += f"-fb{cfg['algoFallback']}-{cfg['destFallback']}"
    out_dir = RESULTS_DIR / tag
    if out_dir.exists():
        shutil.rmtree(out_dir)   # ensure fresh directory for each run
    out_dir.mkdir(exist_ok=True)
    log_file = out_dir / "log.txt"

    start = time.time()
    try:
        proc = subprocess.run(
            build_cmd(cfg, out_dir),
            stdout=log_file.open("w"),
            stderr=subprocess.STDOUT,
            timeout=cfg["T"] + SAFE_MARGIN, # hard wall‑time per run
        )
    except subprocess.TimeoutExpired:
        # kill the still‑running solver and mark as timeout
        proc = subprocess.CompletedProcess(args=[], returncode=-9)
        print(f"[TIMEOUT] {tag} exceeded {cfg['T']+SAFE_MARGIN}s, skipping …",
              file=sys.stderr)
    elapsed = time.time() - start
    if proc.returncode != 0:
        print(f"[WARN] {tag} exited with code {proc.returncode}", file=sys.stderr)

    stats, curve = parse_log(log_file)
    # save SoC curve and create a PNG plot
    if curve:
        (out_dir / "soc.csv").write_text("\n".join(map(str, curve)))
        # quick line chart
        plt.figure()
        plt.plot(curve)
        plt.xlabel("Iteration")
        plt.ylabel("Sum of Costs")
        plt.title(f"{tag}: SOC vs. iteration")
        plt.tight_layout()
        plt.savefig(out_dir / "soc.png")
        plt.close()

    # The algo column in the CSV shows "SAT" whenever dest=="SAT". However, internally it remains --replanAlgo=PP so the binary doesn't run with an invalid value
    csv_algo = "SAT" if cfg["dest"] == "SAT" else cfg["algo"]
    records.append({**cfg, "algo": csv_algo, "run_id": tag, **stats})
    print(f"[INFO] {tag:<70} {elapsed:5.1f}s")

# CSV EXPORT 

df = pd.DataFrame(records)
df.to_csv(RESULTS_DIR / "results_all.csv", index=False)
df[df.kind == "PURE"].to_csv(RESULTS_DIR / "results_T1_pure.csv", index=False)
df[df.kind == "MIX" ].to_csv(RESULTS_DIR / "results_T2_mix.csv",  index=False)

df_runtime = df[[
    "map", "inst", "k",
    "algo", "satProb", "destFallback", "algoFallback",
    "sat_runtime", "other_runtime", "sat_ratio_ops", "sat_ratio"
]]
df_runtime.to_csv(RESULTS_DIR / "results_T3_runtime.csv", index=False)

print("CSV files written to", RESULTS_DIR)