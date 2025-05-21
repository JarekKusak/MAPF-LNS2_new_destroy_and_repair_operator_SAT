# MAPF-LNS2
![test_ubuntu](https://github.com/Jiaoyang-Li/MAPF-LNS2/actions/workflows/test_ubuntu.yml/badge.svg)
![test_macos](https://github.com/Jiaoyang-Li/MAPF-LNS2/actions/workflows/test_macos.yml/badge.svg)

MAPF-LNS2: Fast Repairing for Multi-Agent Path Finding via Large Neighborhood Search

MAPF-LNS2 is an efficient algorithm for solving Multi-Agent Path Finding (MAPF).  
More details can be found in our AAAI 2022 paper [1].

This software extends [MAPF-LNS](https://github.com/Jiaoyang-Li/MAPF-LNS) and includes 
the Anytime variant presented in IJCAI 2021 [2].

---

## 🔧 Build Instructions

### 📦 Dependencies

Required system packages:
- `libboost-all-dev`
- `libeigen3-dev`
- `libgmp-dev`
- `zlib1g-dev`
  
Install on Ubuntu:
```bash
sudo apt update
sudo apt install libboost-all-dev libeigen3-dev libgmp-dev zlib1g-dev
```

You must also compile and include static libraries from the https://github.com/svancaj/MAPF-encodings repository.

In particular, run:
```bash
make lib
```
Then copy:

- All lib*.a from MAPF-encodings/release/libs/ → into MAPF-LNS2/libs/
- The file MAPF.hpp → into MAPF-LNS2/include/

###  Linux (Recommended with Ninja)

```bash
cmake -S . -B build-linux -G Ninja -DCMAKE_BUILD_TYPE=Release
cmake --build build-linux -j$(nproc)
cd build-linux
./lns -m ../random-32-32-20.map -a ../random-32-32-20-random-1.scen -o test -k 150 --outputPaths=paths.txt --destoryStrategy=SAT --maxIterations 300
```

## Usage
Basic run:

```bash
./lns -m <map> -a <scenario> -o <output_prefix> -k <num_agents> \
      --outputPaths=paths.txt --destoryStrategy=SAT --maxIterations 20
```
Example:

```bash
./lns -m random-32-32-20.map -a random-32-32-20-random-1.scen \
      -o test -k 150 --outputPaths=paths.txt --destoryStrategy=SAT --maxIterations 20
```

Arguments:
- -m: map file from MAPF benchmark
- -a: scenario file
- -o: output prefix
- -k: number of agents
- --outputPaths: write computed paths to this file
- --destoryStrategy: strategy to break the solution (SAT / others)
- --maxIterations: number of LNS iterations

More help:

```bash
./lns --help
```

## Test Data
We provide:

- `random-32-32-20.map`
- `random-32-32-20-random-1.scen`

More MAPF instances available at movingai.com.
The scen format is explained here.

To test k agents, the first k rows of the .scen file are used.

## Credits
Developed by Jiaoyang Li and Zhe Chen, based on MAPF-LNS.
PIBT-based solvers integrated from Kei18/pibt.

Licensed under USC Research License. See license.txt.

## References
[1] Jiaoyang Li, Zhe Chen, Daniel Harabor, Peter J. Stuckey and Sven Koenig.
MAPF-LNS2: Fast Repairing for Multi-Agent Path Finding via Large Neighborhood Search, AAAI 2022.

[2] Jiaoyang Li et al.
Anytime Multi-Agent Path Finding via Large Neighborhood Search, IJCAI 2021.