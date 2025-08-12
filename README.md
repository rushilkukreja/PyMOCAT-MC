# PyMOCAT-MC

**PyMOCAT-MC** is a Python implementation of the MIT Orbital Capacity Assessment Toolbox - Monte Carlo (MOCAT-MC), originally written in MATLAB. This toolbox enables Monte Carlo-based analysis of orbital capacity, debris, and satellite population evolution, supporting research and engineering in space situational awareness and orbital risk. The Python implementation maintains full compatibility with the original MATLAB version while leveraging modern Python scientific computing libraries.

---

## Software Purpose
PyMOCAT-MC provides a flexible, extensible framework for simulating the evolution of satellite populations and debris in Earth orbit. It is designed for:
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
- pandas >= 1.3.0
- matplotlib >= 3.5.0
- Python >= 3.8

### Installation

It is recommended to use a virtual environment:

```bash
cd python_implementation
python -m venv mocat_env
source mocat_env/bin/activate  # On Windows: mocat_env\Scripts\activate
pip install numpy scipy pandas matplotlib
```

---
## Usage
### Quick Start

Example scripts are provided in the `examples/` directory:

```bash
cd python_implementation
python examples/Quick_Start/quick_start.py
```

For a scenario without launches:

```bash
python examples/Scenario_No_Launch/scenario_no_launch.py
```

For a realistic 2020s scenario with megaconstellations:

```bash
python examples/Realistic_Scenario/realistic_scenario_2020s.py
```

### Custom Configuration

```python
from mocat_mc import MOCATMC

# Initialize MOCAT-MC
mocat = MOCATMC()

# Setup configuration
cfg_mc = mocat.setup_mc_config(rng_seed=1, ic_file='2020.mat')

# Customize simulation parameters
cfg_mc['tspan'] = 30  # Simulation duration in years
cfg_mc['N_ini_samples'] = 10  # Number of initial orbit samples
cfg_mc['prop']['dt_TL'] = 90  # Time step in days

# Run simulation
results = mocat.run_simulation(cfg_mc)

print(f"Final counts - Satellites: {results['nS']}, Derelicts: {results['nD']}, Debris: {results['nN']}, Rocket Bodies: {results['nB']}")
```

---

## Data

Supporting data files are provided in the `python_implementation/supporting_data/` directory:

- **TLE Historic Data** (`TLEhistoric/`): CSV and MAT files containing Two-Line Element sets for satellites from 2000–2023. See the included `readme.txt` for details and data sources.
- **Atmospheric Density Model**: `dens_jb2008_032020_022224.mat` - JB2008 atmospheric density data
- **Launch Schedules**: `megaconstellationLaunches.xlsx` - Planned megaconstellation launch data

---

## Project Structure

```
PyMOCAT-MC-2/
├── python_implementation/
│   ├── mocat_mc.py                    # Main MOCAT-MC class
│   ├── examples/
│   │   ├── Quick_Start/
│   │   │   ├── quick_start.py         # Basic simulation example
│   │   │   └── quick_start_with_plots.py  # Example with visualizations
│   │   ├── Realistic_Scenario/
│   │   │   └── realistic_scenario_2020s.py  # 2020s megaconstellation scenario
│   │   └── Scenario_No_Launch/
│   │       └── scenario_no_launch.py  # Decay-only scenario
│   ├── supporting_functions/
│   │   ├── cfg_mc_constants.py        # Physical constants
│   │   ├── get_idx.py                 # Matrix index definitions
│   │   ├── categorize_obj.py          # Object categorization
│   │   ├── init_sim.py                # Simulation initialization
│   │   ├── main_mc.py                 # Main Monte Carlo engine
│   │   ├── prop_mit_vec.py            # MIT orbital propagator
│   │   ├── collision_prob_vec.py      # Collision probability
│   │   ├── cube_vec_v3.py             # Cube method for conjunctions
│   │   ├── frag_col_sbm_vec.py        # Collision fragmentation model
│   │   ├── orbcontrol_vec.py          # Orbit control/station-keeping
│   │   ├── fillin_atmosphere.py       # Atmospheric model interface
│   │   ├── fillin_physical_parameters.py  # Object physical properties
│   │   └── [additional utilities...]
│   └── supporting_data/
│       ├── TLEhistoric/               # Historical TLE data (2000-2023)
│       ├── dens_jb2008_032020_022224.mat  # Atmospheric density data
│       └── megaconstellationLaunches.xlsx  # Launch schedules
├── MATLAB_implementation/             # Original MATLAB implementation with additional test files
├── comparison_tests/                  # Python vs MATLAB comparison tests
│   ├── test_all_scenarios.py          # Comprehensive comparison across all scenarios
│   ├── test_quick_scenarios.py        # Quick validation tests for basic scenarios
│   ├── accuracy_error_data.csv        # Accuracy error measurements data
│   └── accuracy_error_data.json       # Accuracy error data in JSON format
├── analysis_figures/                  # Analysis and visualization scripts
│   ├── plot_total_object_counts.py    # Generate total object count comparison plots
│   ├── plot_execution_time.py         # Visualize execution time performance
│   ├── plot_object_type_heatmap.py    # Create heatmap of object type percentages
│   ├── plot_computational_efficiency.py # Analyze computational efficiency metrics
│   ├── measure_accuracy_errors.py     # Calculate accuracy errors between implementations
│   ├── plot_error_analysis.py         # Generate error distribution plots
│   ├── total_object_counts.png        # Output: Object count comparison visualization
│   ├── execution_time_comparison.png  # Output: Performance comparison chart
│   ├── object_type_percentage_heatmap.png # Output: Object type distribution heatmap
│   ├── computational_efficiency_analysis.png # Output: Efficiency analysis chart
│   └── error_box_plots.png            # Output: Error distribution visualization
└── README.md                          # This file
```

---

## Contributors

- Rushil Kukreja, Thomas Jefferson High School for Science and Technology
- Dr. Edward J. Oughton, George Mason University
- Dr. Giovanni Lavezzi, Massachusetts Institute of Technology
- Dr. Enrico M. Zucchelli, Massachusetts Institute of Technology
- Dr. Daniel Jang, Lincoln Labs
- Dr. Richard Linares, Massachusetts Institute of Technology