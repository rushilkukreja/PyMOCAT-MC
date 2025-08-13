# PyMOCAT-MC Test Suite

This directory contains the complete test suite for PyMOCAT-MC.

## Test Organization

### Core Tests (`/python_implementation/tests/`)
Essential tests for the Python implementation:
- `test_basic_functionality_simple.py` - Core functionality tests
- `test_import.py` - Import and data loading tests  
- `minimal_test.py` - Minimal 2-step simulation
- `test_simple_run.py` - Single time-step test

### Comparison Tests (`/comparison_tests/`)
Cross-validation tests between Python and MATLAB:
- `test_all_scenarios.py` - Full scenario testing
- `test_quick_scenarios.py` - Quick comparison tests

### MATLAB Tests (`/MATLAB_implementation/tests/`)
Original MATLAB test suite:
- Various `.m` test files for MATLAB validation

## Running Tests

### Quick Testing (Recommended)
```bash
# Run all essential tests (~11 seconds)
python3 run_tests.py

# Or using bash
./run_tests.sh
```

### Individual Tests
```bash
# Basic functionality (~1s)
python3 python_implementation/tests/test_basic_functionality_simple.py

# Import validation (~2s)  
python3 python_implementation/tests/test_import.py

# Minimal simulation (~5s)
python3 python_implementation/tests/minimal_test.py

# Simple run (~3s)
python3 python_implementation/tests/test_simple_run.py
```

### Comparison Tests
```bash
# Cross-validation with MATLAB (longer running)
python3 comparison_tests/test_quick_scenarios.py
```

## Test Coverage

The test suite validates:
- ✅ Core imports and dependencies
- ✅ Data file access (TLE historic data)
- ✅ Basic simulation pipeline
- ✅ Orbital mechanics calculations
- ✅ Object categorization
- ✅ End-to-end functionality
- ✅ Example script execution

## Expected Results

All tests should pass with output:
```
✅ All tests passed!
```

For troubleshooting, see `../TESTING.md`.