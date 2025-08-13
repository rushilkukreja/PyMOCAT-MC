# Testing PyMOCAT-MC

## Quick Start

### Automated Setup (Recommended)
```bash
# Create a clean environment with all dependencies
python3 setup_environment.py

# Activate the environment (follow instructions from setup script)
source pymocat_env/bin/activate  # On Unix/macOS
# or
pymocat_env\Scripts\activate     # On Windows

# Run tests
python3 run_tests.py
```

### Manual Setup
```bash
# Install dependencies
pip install -r requirements.txt

# Run tests
python3 run_tests.py
```

## Test Runners

### Method 1: Python test runner (Cross-platform)
```bash
python3 run_tests.py
```

### Method 2: Bash test runner (Unix/macOS)
```bash
./run_tests.sh
```

### Method 3: Individual tests
```bash
# Basic functionality (~1s)
python3 python_implementation/tests/test_basic_functionality_simple.py

# Import test (~2s)
python3 python_implementation/tests/test_import.py

# Minimal simulation (~5s)
python3 python_implementation/tests/minimal_test.py

# Simple run test (~3s)
python3 python_implementation/tests/test_simple_run.py
```

## Test Structure

The test suite includes 4 essential tests that provide complete coverage:

1. **Basic Functionality Tests** (`test_basic_functionality_simple.py`)
   - Tests core imports and functions
   - Verifies data structures (get_idx, categorize_obj)
   - Tests MOCATMC class initialization
   - Fastest test (~1 second)

2. **Import Test** (`test_import.py`)
   - Verifies all required modules can be imported
   - Checks TLE data loading (finds 14,207 satellites)
   - Tests MOCATMC initialization
   - Validates file paths (~2 seconds)

3. **Minimal Simulation Test** (`minimal_test.py`)
   - Runs a minimal 2-step simulation
   - Verifies basic simulation pipeline
   - Tests end-to-end functionality (~5 seconds)

4. **Simple Run Test** (`test_simple_run.py`)
   - Runs a single time step simulation
   - Tests core simulation functionality
   - Quick smoke test (~3 seconds)

5. **Example Scripts** (run as integration tests)
   - Quick Start example (full year simulation)
   - Scenario No Launch example
   - Realistic Scenario example

## Requirements

### System Requirements
- Python 3.8 or higher
- NumPy
- SciPy
- Matplotlib (for plotting examples)
- Pandas (optional, for data handling)

### Installing Dependencies

```bash
# Using pip
pip3 install numpy scipy matplotlib pandas

# Or using the requirements file
pip3 install -r requirements.txt
```

## Virtual Environment Setup (Optional)

If you want to use a virtual environment:

```bash
# Create a new virtual environment
python3 -m venv test_env

# Activate it
source test_env/bin/activate  # On Unix/macOS
# or
test_env\Scripts\activate  # On Windows

# Install dependencies
pip install -r requirements.txt

# Run tests
python run_tests.py
```

## Known Issues and Solutions

### Issue 1: Import Errors
If you see `ModuleNotFoundError: No module named 'mocat_mc'`:
- Make sure you're running tests from the repository root directory
- The test scripts automatically add the correct paths

### Issue 2: pytest hangs
If pytest hangs during test collection:
- Use the provided `run_tests.py` script instead
- Or run individual test files directly with Python

### Issue 3: Virtual Environment Issues
If the virtual environment doesn't work:
- Use your system Python 3 installation
- Install dependencies globally with `pip3 install --user`

### Issue 4: Missing TLE Data
If tests fail due to missing data files:
- Ensure the `supporting_data/TLEhistoric/` directory contains `.mat` files
- At minimum, `2020.mat` should be present for tests to run

## Expected Output

When all tests pass, you should see:
```
+ Basic functionality tests: PASSED
+ Import test: PASSED
+ Minimal simulation test: PASSED
+ Simple run test: PASSED
+ Quick Start example: PASSED

All tests passed!
```

## Warnings

You may see some warnings during test execution:
- `RuntimeWarning: invalid value encountered in cast` - This is expected for NaN handling
- `RuntimeWarning: divide by zero` - This occurs in edge cases and is properly handled

These warnings do not indicate test failures.

## Continuous Integration

For CI/CD pipelines, use:
```bash
python3 run_tests.py
```

The script returns:
- Exit code 0: All tests passed
- Exit code 1: One or more tests failed

## Troubleshooting

### Common Issues and Solutions

#### Package: Missing Dependencies
```
Error: ModuleNotFoundError: No module named 'numpy'
```
**Solution:**
```bash
pip install -r requirements.txt
# or use the automated setup
python3 setup_environment.py
```

#### Directory: Wrong Working Directory
```
Error: No such file or directory: 'python_implementation/...'
```
**Solution:**
```bash
# Make sure you're in the repository root
cd PyMOCAT-MC-2
python3 run_tests.py
```

#### Python: Python Version Issues
```
Error: Python 3.8 or higher is required
```
**Solution:**
- Install Python 3.8+ from [python.org](https://python.org)
- On Ubuntu: `sudo apt install python3.8`
- On macOS: `brew install python@3.8`

#### Timeout: Tests Hanging/Timeout
```
Error: Test timed out (>120 seconds)
```
**Solution:**
- Check if virtual environment is corrupted
- Try running: `python3 setup_environment.py`
- Run individual tests to isolate the issue

#### Environment: Virtual Environment Issues
```
Error: Virtual environment activation fails
```
**Solution:**
```bash
# Remove old environments and create fresh one
rm -rf venv test_env verify_env pymocat_env
python3 setup_environment.py
```

#### Data: Missing TLE Data
```
Warning: Could not find 2020.mat file
```
**Solution:**
- Tests will work with synthetic data (reduced functionality)
- This is expected for most users - not a failure
- Full TLE data would need to be provided separately

#### Windows: Windows-Specific Issues
```
Error: 'python3' is not recognized
```
**Solution:**
```cmd
# Use 'python' instead of 'python3' on Windows
python run_tests.py

# Or use the Python launcher
py -3 run_tests.py
```

### Getting Help

If tests still fail:
1. Check Python version: `python3 --version` (should be 3.8+)
2. Try the automated setup: `python3 setup_environment.py`
3. Run with verbose output: `python3 run_tests.py > test_output.log 2>&1`
4. Check the log file for detailed error messages
5. Open an issue on GitHub with the log file

### Platform-Specific Notes

**macOS:**
- Use `python3` command
- May need to install Xcode command line tools: `xcode-select --install`

**Ubuntu/Debian:**
```bash
sudo apt update
sudo apt install python3 python3-pip python3-venv
```

**Windows:**
- Use `python` instead of `python3`
- Install from [python.org](https://python.org) or Microsoft Store
- Use PowerShell or Command Prompt

**CentOS/RHEL:**
```bash
sudo yum install python3 python3-pip
# or on newer versions:
sudo dnf install python3 python3-pip
```