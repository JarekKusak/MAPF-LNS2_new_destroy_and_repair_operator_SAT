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
MAPS = {"Paris_1_256"}#, "random-32-32-20",  "warehouse-20-40-10-2-1", "lt_gallowstemplar_n", "room-64-64-16"}
INSTANCES_PER_MAP = 10
AGENT_COUNTS      = [100]#, 200, 300, 400]#, 200, 300, 400] #[100, 500, 1000] # rozdělit na malé a velké mapy (maximálně 15 s pro initial solution)
TIMEOUTS          = [30] 
SUBMAP_SIDES      = [5]
MIX_PROBS         = [50]#[80, 50, 20] # 0 and 100 is generated outside the MIX
PURE_REPLANS     = ["PP"]  # pure replanners to test
SAT_HEURISTICS    = ["adaptive"]
FALLBACK_DESTS    = ["Adaptive"]
FALLBACK_ALGOS    = ["PP"]
MAX_ITERS         = [1_000_000] # arbitrary long
INCLUDE_PURE_SAT  = True

# ──────────────────────────────────────────────────────────────────────────────
# Runtime limits
#   T  … algorithm budget handed to the solver via -t (soft deadline)
#   Hard kill = T + SAFE_MARGIN:   **only** a last‑resort fuse when the solver
#     neither respects the soft deadline nor terminates within a very generous
#     clean‑up window (I/O, destructors, etc.).
#
# In practice the solver finishes normal runs within ≈10 s of bookkeeping after
# hitting its soft  -t  limit.  We therefore give it a full minute of slack.
# ──────────────────────────────────────────────────────────────────────────────
SAFE_MARGIN = 60        # seconds – generous clean‑up window, still guards hangs
def wall_clock_limit(T):
    """Return wall-clock limit = soft limit T  + SAFE_MARGIN.
    Soft deadline is enforced by solver; hard kill only catches infinite loops."""
    return T + SAFE_MARGIN

 # Path to the compiled solver; adjust per platform
LNS_BIN = "./build-macos/lns"   # use ./build-linux/lns on Linux
RESULTS_DIR = Path("results").absolute()
RESULTS_DIR.mkdir(exist_ok=True)

# REGEX PATTERNS (must match solver output)

RE_SAT_RT  = re.compile(
    r"\[STAT\]\s+SAT total runtime\s*=\s*([\d\.eE+-]+)\s+s")
RE_OTH_RT  = re.compile(
    r"\[STAT\]\s+Other operators runtime\s*=\s*([\d\.eE+-]+)\s+s")
RE_RS = re.compile(
    r"^\[STAT\]\s*wall_runtime\s*=\s*([\d\.eE+-]+).*?"
    r"sat_time\s*=\s*([\d\.eE+-]+).*?"
    r"other_time\s*=\s*([\d\.eE+-]+).*?"
    r"overhead\s*=\s*([\d\.eE+-]+).*?"
    r"sat_calls\s*=\s*(\d+).*?"
    r"sat_iters\s*=\s*(\d+).*?"
    r"sat_ok\s*=\s*(\d+).*?"
    r"sat_fail\s*=\s*(\d+).*?"
    r"final_soc\s*=\s*(\d+).*?"
    r"initial_soc\s*=\s*(\d+).*?"
    r"failed_iterations\s*=\s*(\d+).*?"
    r"outer_iterations\s*=\s*(\d+)"
)
RE_SOC_POST = re.compile(r"\[STAT\] sum_of_costs after recomputation: (\d+)")

RE_SOC_INLINE = re.compile(r"\[SOC\]\s+(\d+)")

# CASE GENERATION

def main():
    cases: list[dict] = []

    # T1 ─ pure strategies (PP / CBS / optional pure SAT)
    for m, scen_i, k, T, iters, sub in product(
            MAPS, range(1, INSTANCES_PER_MAP + 1),
            AGENT_COUNTS, TIMEOUTS, MAX_ITERS, SUBMAP_SIDES):
        for repl in PURE_REPLANS:
            cases.append(dict(kind="PURE", algo=repl, satProb=0,
                              dest="Adaptive",   # any non-SAT strategy
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

    # MAIN LOOP

    # --- parallel execution -------------------------------------------------
    import os
    from concurrent.futures import ProcessPoolExecutor, as_completed

    MAX_PARALLEL = max(1, os.cpu_count() - 2)   # keep some cores free

    print(f"Running {len(cases)} cases with up to {MAX_PARALLEL} workers …",
          file=sys.stderr)

    records: list[dict] = []
    with ProcessPoolExecutor(max_workers=MAX_PARALLEL) as pool:
        future_to_cfg = {pool.submit(run_case, cfg): cfg for cfg in cases}
        for fut in as_completed(future_to_cfg):
            rec = fut.result()
            records.append(rec)
            print("[DONE]", rec["run_id"], file=sys.stderr)

    # quick SOC plots per run
    for rec in records:
        soc_path = RESULTS_DIR / rec["run_id"] / "soc.csv"
        if soc_path.exists():
            curve = pd.read_csv(soc_path, header=None).squeeze("columns")
            plt.figure()
            plt.plot(curve)
            plt.xlabel("Iteration")
            plt.ylabel("Sum of Costs")
            plt.title(f"{rec['run_id']}: SOC vs. iteration")
            plt.tight_layout()
            plt.savefig(soc_path.with_suffix('.png'))
            plt.close()

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
        cmd.append(f"--destoryStrategyFallback={cfg.get('destFallback', 'Adaptive')}") # if destFallback is not set, use default Intersection
        cmd.append(f"--replanAlgoFallback={cfg.get('algoFallback', 'PP')}") # if algoFallback is not set, use default PP
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
            if m := RE_SAT_RT.search(line):
                stats["sat_runtime"] = float(m[1])
            elif m := RE_OTH_RT.search(line):
                stats["other_runtime"] = float(m[1])
            elif m := RE_RS.search(line):
                (stats["runtime"],
                 stats["sat_runtime"],
                 stats["other_runtime"],
                 stats["overhead_runtime"],
                 stats["sat_calls"],
                 stats["sat_iters"],
                 stats["sat_ok"],
                 stats["sat_fail"],
                 stats["final_soc"],
                 stats["initial_soc"],
                 stats["failed_iterations"],
                 stats["outer_iterations"]) = (
                    float(m[1]), float(m[2]), float(m[3]), float(m[4]),
                    int(m[5]), int(m[6]), int(m[7]), int(m[8]),
                    int(m[9]), int(m[10]), int(m[11]), int(m[12])
                )
            elif m := RE_SOC_INLINE.search(line):
                curve.append(int(m[1]))

    # defaults
    stats.setdefault("sat_runtime", 0.0)
    stats.setdefault("other_runtime", 0.0)
    stats.setdefault("runtime", 0.0)   # total wall‑clock time of the run
    stats.setdefault("overhead_runtime", 0.0)
    stats.setdefault("sat_iters", 0)
    stats.setdefault("outer_iterations", 0)
    stats.setdefault("failed_iterations", 0)

    # Derived fields:
    # (1) Share of **successful SAT operator calls** among all SAT calls
    stats["sat_ratio_ops"] = (
        100.0 * stats["sat_runtime"] /
        (stats["sat_runtime"] + stats["other_runtime"])
        if (stats["sat_runtime"] + stats["other_runtime"]) else 0.0
    )

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

# ----------------------------------------------------------------------
# PARALLEL EXECUTION ----------------------------------------------------
def run_case(cfg: dict) -> dict:
    """
    Execute one solver run in its own working directory and
    return the record that will later be appended to the CSV.
    Designed to be picklable for ProcessPoolExecutor.
    """
    tag = (f"{cfg['map']}-i{cfg['inst']}-k{cfg['k']}-t{cfg['T']}"
           f"-sub{cfg['sub']}-p{cfg['satProb']}"
           f"-{cfg.get('satHeur', '')}")
    if cfg.get("destFallback"):
        tag += f"-fb{cfg['algoFallback']}-{cfg['destFallback']}"

    out_dir = RESULTS_DIR / tag
    if out_dir.exists():
        shutil.rmtree(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    log_file = out_dir / "log.txt"

    try:
        subprocess.run(
            build_cmd(cfg, out_dir),
            stdout=log_file.open("w"),
            stderr=subprocess.STDOUT,
            timeout=wall_clock_limit(cfg["T"]),
            check=False,
        )
    except subprocess.TimeoutExpired:
        return {**cfg, "run_id": tag, "state": "timeout"}

    stats, curve = parse_log(log_file)

    # store SoC curve for later plotting
    if curve:
        (out_dir / "soc.csv").write_text("\n".join(map(str, curve)))

    csv_algo = "SAT" if cfg["dest"] == "SAT" else cfg["algo"]
    return {**cfg, "algo": csv_algo, "run_id": tag, "state": "OK", **stats}

# Standard multiprocessing guard for Windows compatibility
if __name__ == "__main__":
    import multiprocessing as mp
    mp.freeze_support()   # for Windows; no‑op elsewhere
    main()