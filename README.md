# MAPF-LNS2
![test_ubuntu](https://github.com/Jiaoyang-Li/MAPF-LNS2/actions/workflows/test_ubuntu.yml/badge.svg)
![test_macos](https://github.com/Jiaoyang-Li/MAPF-LNS2/actions/workflows/test_macos.yml/badge.svg)

MAPF-LNS2: Fast Repairing for Multi-Agent Path Finding via Large Neighborhood Search


MAPF-LNS2 is an efficient algorithm for solving Multi-Agent Path Finding (MAPF). 
More details can be found in our paper at AAAI 2022 [1].

This software is an advanced version and a superset of MAPF-LNS (https://github.com/Jiaoyang-Li/MAPF-LNS); that is, it also contains Anytime Multi-Agent Path Finding via Large Neighborhood Search [2]. 

## Build & Run

Below we give **two** concise build recipes – one for **macOS** (Home‑brew + CMake) and one for **Linux** (APT + Ninja).  
Both assume you already cloned the repository **after** the library folders have been split into  
`libs-macos/` and `libs-linux/` as described in the build notes.

---

### macOS (Apple silicon / Intel)

```bash
# Packages
brew install boost eigen gmp zlib cmake

# Configure + build
cmake -B build-macos -DCMAKE_BUILD_TYPE=Release
cmake --build build-macos -j$(sysctl -n hw.ncpu)

# Example run
build-macos/lns \
  -m random-32-32-20.map \
  -a random-32-32-20-random-1.scen \
  -o test -k 150 --outputPaths=paths.txt \
  --destoryStrategy=SAT --maxIterations 20
```

---

### Linux (Ubuntu 20.04/22.04)

```bash
# Packages
sudo apt update
sudo apt install build-essential ninja-build \
                 libboost-all-dev libeigen3-dev \
                 libgmp-dev zlib1g-dev cmake

# Configure with Ninja (faster)
cmake -S . -B build-linux -G Ninja -DCMAKE_BUILD_TYPE=Release
cmake --build build-linux

# Example run
build-linux/lns \
  -m random-32-32-20.map \
  -a random-32-32-20-random-1.scen \
  -o test -k 150 --outputPaths=paths.txt \
  --destoryStrategy=SAT --maxIterations 20
```

> **Note:**  
> The static libraries from *MAPF‑encodings* are pre‑compiled for each OS:  
> `libs-macos/` is used automatically on macOS, `libs-linux/` on Linux.

---

### Command‑line options

- `-m` : map file (MovingAI format)  
- `-a` : scenario file  
- `-o` : prefix for output files (no extension)  
- `-k` : number of agents to load from the scenario  
- `-t` : runtime limit in seconds  
- `--outputPaths` : file to which the final paths are written  
- `--destoryStrategy` / `--maxIterations` : algorithmic parameters

Run `./lns --help` to see the full list and defaults.

---

## Credits

The software was developed by Jiaoyang Li and Zhe Chen based on [MAPF-LNS](https://github.com/Jiaoyang-Li/MAPF-LNS).

The rule-based MAPF solvers (i.e., PPS, PIBT, and winPIBT) inside the software were borrowed from 
https://github.com/Kei18/pibt/tree/v1.3

MAPF-LNS2 is released under USC – Research License. See license.txt for further details.
 
## References
[1] Jiaoyang Li, Zhe Chen, Daniel Harabor, Peter J. Stuckey and Sven Koenig.
MAPF-LNS2: Fast Repairing for Multi-Agent Path Finding via Large Neighborhood Search
In Proceedings of the AAAI Conference on Artificial Intelligence, (in print), 2022.

[2] Jiaoyang Li, Zhe Chen, Daniel Harabor, Peter J. Stuckey, Sven Koenig. 
Anytime Multi-Agent Path Finding via Large Neighborhood Search. 
In Proceedings of the International Joint Conference on Artificial Intelligence (IJCAI), pages 4127-4135, 2021.
