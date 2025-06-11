#!/usr/bin/env python3
"""
tables.py
─────────
• Enumerates all benchmark runs for the bachelor-thesis experiments.
• Launches `lns` with the chosen parameter matrix.
• Parses each run’s log to collect the headline statistics required for
  – T1 : pure re-planning strategies (PP, CBS) and pure SAT,
  – T2 : mixed SAT + PP with different satProb values.
• Produces three CSV files in <repo>/results/:
      ├─ results_all.csv   – every single run
      ├─ results_T1_pure.csv
      └─ results_T2_mix.csv
"""

from itertools import product
from pathlib import Path
import re
import subprocess
import time
import sys
import pandas as pd
import matplotlib.pyplot as plt

# ------------------ CONFIGURATION MATRIX (edit as needed) ------------------ #

MAPS = {
    "warehouse-10-20-10-2-1",
    "random-32-32-20",
    "room-32-32-4",
    "Paris_1_256",
    "ost003d",
}

# --- CONFIGURATION MATRIX (quick sanity run) ---
MAPS              = {"ost003d"}
INSTANCES_PER_MAP = 1
AGENT_COUNTS      = [100]
TIMEOUTS          = [5]          # 5 s wall-time
MAX_ITERS         = [5_000]      # libovolně velké, stejně to utneme časem
SUBMAP_SIDES      = [3]

PURE_REPLANS      = ["PP"]
INCLUDE_PURE_SAT  = True
MIX_PROBS         = [100, 50, 20]     # 100 % SAT   ×  0 % SAT
SAT_HEURISTICS    = ["adaptive"]

SAFE_MARGIN = 3  # seconds added on top of cfg['T'] to forcibly kill hanging runs
'''
INSTANCES_PER_MAP = 5          # scenario-files random-1..random-5
AGENT_COUNTS      = [100, 500]
TIMEOUTS          = [2]       # seconds; thesis advisor suggested 30 s
SAFE_MARGIN = 3  # seconds added on top of cfg['T'] to forcibly kill hanging runs
MAX_ITERS         = [10_000]   # “high enough” – loop is limited by timeout
SUBMAP_SIDES      = [3, 5]

PURE_REPLANS = ["PP", "CBS"]   # T1 without SAT
INCLUDE_PURE_SAT = True
MIX_PROBS    = [100, 50, 20]   # satProb for T2
SAT_HEURISTICS = ["adaptive"]  # you can add ["roundRobin","mostDelayed"]
'''
LNS_BIN     = "./lns"          # path to compiled solver
RESULTS_DIR = Path("results").absolute()
RESULTS_DIR.mkdir(exist_ok=True)

# ------------------ REGEX PATTERNS (must match solver output) -------------- #

RE_FINAL   = re.compile(
    r"\[STAT\] .*: runtime = ([\d\.eE+-]+), iterations = (\d+), "
    r"solution cost = (\d+), initial solution cost = (\d+), failed iterations = (\d+)"
)
RE_SAT_RT  = re.compile(r"\[STAT\] SAT total runtime = ([\d\.eE+-]+) s")
RE_OTH_RT  = re.compile(r"\[STAT\] Other operators runtime = ([\d\.eE+-]+) s")
RE_SOC_POST = re.compile(r"\[STAT\] sum_of_costs after recomputation: (\d+)")

# ────────────────────────── CASE GENERATION ─────────────────────────────── #

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
for m, scen_i, k, T, iters, sub, prob, heur in product(
        MAPS, range(1, INSTANCES_PER_MAP + 1), AGENT_COUNTS,
        TIMEOUTS, MAX_ITERS, SUBMAP_SIDES, MIX_PROBS, SAT_HEURISTICS):
    if prob == 0:
        continue
    cases.append(dict(kind="MIX", algo="PP", satProb=prob,
                      dest="SAT", satHeur=heur,
                      map=m, inst=scen_i, k=k, T=T, iters=iters, sub=sub))

# ──────────────────────── LAUNCH / PARSE HELPERS ────────────────────────── #

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
        "--satDebug=0",                     # keep logs small
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

    return stats, curve

# ───────────────────────────── MAIN LOOP ────────────────────────────────── #

records = []
for cfg in cases:
    tag = (f"{cfg['map']}-i{cfg['inst']}-k{cfg['k']}-t{cfg['T']}"
           f"-sub{cfg['sub']}-p{cfg['satProb']}"
           f"-{cfg.get('satHeur','')}")
    out_dir = RESULTS_DIR / tag
    out_dir.mkdir(exist_ok=True)
    log_file = out_dir / "log.txt"

    start = time.time()
    try:
        proc = subprocess.run(
            build_cmd(cfg, out_dir),
            stdout=log_file.open("w"),
            stderr=subprocess.STDOUT,
            timeout=cfg["T"] + SAFE_MARGIN,          # hard wall‑time per run
        )
    except subprocess.TimeoutExpired:
        # Kill the still‑running solver and mark as timeout
        proc = subprocess.CompletedProcess(args=[], returncode=-9)
        print(f"[TIMEOUT] {tag} exceeded {cfg['T']+SAFE_MARGIN}s, skipping …",
              file=sys.stderr)
    elapsed = time.time() - start
    if proc.returncode != 0:
        print(f"[WARN] {tag} exited with code {proc.returncode}", file=sys.stderr)

    stats, curve = parse_log(log_file)
    # Save SoC curve and create a PNG plot
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
    records.append({**cfg, "run_id": tag, **stats})
    print(f"[INFO] {tag:<70} {elapsed:5.1f}s")

# ───────────────────────────── CSV EXPORT ───────────────────────────────── #

df = pd.DataFrame(records)
df.to_csv(RESULTS_DIR / "results_all.csv", index=False)
df[df.kind == "PURE"].to_csv(RESULTS_DIR / "results_T1_pure.csv", index=False)
df[df.kind == "MIX" ].to_csv(RESULTS_DIR / "results_T2_mix.csv",  index=False)

print("✓ CSV files written to", RESULTS_DIR)