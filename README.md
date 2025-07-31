# MOCAT-MC

**MOCAT-MC** is a Python conversion of the MIT Orbital Capacity Assessment Toolbox - Monte Carlo (MOCAT-MC), originally written in MATLAB. This toolbox enables Monte Carlo-based analysis of orbital capacity, debris, and satellite population evolution, supporting research and engineering in space situational awareness and orbital risk. This Python conversion maintains the same functionality as the original MATLAB code while providing better integration with modern Python scientific computing ecosystems.

---

## Software Purpose

MOCAT-MC provides a flexible, extensible framework for simulating the evolution of satellite populations and debris in Earth orbit. It is designed for:
- Assessing orbital capacity and congestion
- Modeling satellite launches, collisions, and fragmentation
- Evaluating the impact of megaconstellations
- Supporting policy and engineering studies in orbital risk

---

## Setup and Configuration

### Requirements

Key dependencies include:
- numpy >= 1.21.0
- scipy >= 1.7.0
- astropy >= 5.0.0
- pandas >= 1.3.0
- numba >= 0.56.0 (optional)
- matplotlib >= 3.5.0 (optional)
- tqdm >= 4.62.0 (optional)
- joblib >= 1.1.0 (optional)

### Installation

It is recommended to use a virtual environment:

```bash
python -m venv mocat_env
source mocat_env/bin/activate
pip install -r requirements.txt
```

---
## Usage
### Quick Start

Example scripts are provided in the `Examples/Quick_Start/` directory. To run a basic simulation:

```bash
python Examples/Quick_Start/Quick_Start.py
```

Or for a scenario without launches:

```bash
python Examples/Scenario_No_Launch/Scenario_No_Launch.py
```

### Custom Configuration

```python
from Examples.Quick_Start.setup_MCconfig import setup_MCconfig
from Examples.Quick_Start.main_mc import main_mc

# Setup configuration
cfgMC = setup_MCconfig(seed=1, ICfile='2020.mat')

# Run simulation
nS, nD, nN, nB, mat_sats = main_mc(cfgMC, seed=1)

print(f"Final counts - Satellites: {nS}, Derelicts: {nD}, Debris: {nN}, Rocket Bodies: {nB}")
```

---

## Data

Supporting data files are provided in the `supporting_data/` directory, including:
- TLE historic data (`supporting_data/TLEhistoric/`): CSV files of Two-Line Element sets for satellites, 2000–2023. See the included `readme.txt` for details and data sources.
- Atmospheric density and launch data: `dens_jb2008_032020_022224.mat`, `megaconstellationLaunches.xlsx`

---

## Project Structure

```
MOCAT-MC/
├── Examples/
│   ├── Quick_Start/
│   │   ├── Quick_Start.py          # Main entry point for quick start
│   │   ├── setup_MCconfig.py       # Configuration setup
│   │   ├── initSim.py              # Simulation initialization
│   │   └── main_mc.py              # Main Monte Carlo engine
│   └── Scenario_No_Launch/
│       └── Scenario_No_Launch.py   # No launch scenario example
├── supporting_functions/
│   ├── cfgMC_constants.py          # Physical constants
│   ├── getidx.py                   # Matrix index definitions
│   ├── categorizeObj.py            # Object categorization
│   ├── collision_prob_vec.py       # Collision probability calculation
│   ├── cross_vec.py                # Cross product operations
│   ├── jd2date.py                  # Julian date conversions
│   ├── cube_vec_v3.py              # Cube method collision detection
│   ├── prop_mit_vec.py             # MIT orbital propagator
│   ├── orbcontrol_vec.py           # Orbit control functions
│   ├── getZeroGroups.py            # Zero group analysis
│   ├── fillin_physical_parameters.py # Physical parameter filling
│   └── fillin_atmosphere.py        # Atmospheric model setup
├── supporting_data/                # Data files (.mat, .csv, etc.)
├── requirements.txt                # Python dependencies
└── README.md               # This file
```

---

## Contributors

- Rushil Kukreja, Thomas Jefferson High School for Science and Technology
- Dr. Edward J. Oughton, George Mason University
- Dr. Giovanni Lavezzi, Massachusetts Institute of Technology
- Dr. Enrico Zucchelli, Massachusetts Institute of Technology
- Dr. Daniel Zhang, Lincoln Labs
- Dr. Richard Linares, Massachusetts Institute of Technology