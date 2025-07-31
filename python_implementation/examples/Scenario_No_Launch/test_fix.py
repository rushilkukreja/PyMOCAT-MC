#!/usr/bin/env python3
"""
Test script to verify path fixes and file accessibility
Python equivalent of test_fix.m

This script tests whether required data files can be found
and verifies that the Python implementation has access to
all necessary supporting functions and data.

Based on: Examples/Scenario_No_Launch/test_fix.m
"""

import os
import sys
import numpy as np
from pathlib import Path

# Add parent directories to path
script_dir = Path(__file__).parent
examples_dir = script_dir.parent
python_impl_dir = examples_dir.parent
repo_root = python_impl_dir.parent

sys.path.append(str(python_impl_dir))
sys.path.append(str(repo_root))


def test_file_paths():
    """Test if required data files can be found"""
    print("=== Path Fix Test Script ===")
    print(f"Current working directory: {os.getcwd()}")
    print(f"Script location: {script_dir}")
    print(f"Python implementation root: {python_impl_dir}")
    print(f"Repository root: {repo_root}")
    
    # Test 1: Check if density file can be found
    print("\nTest 1: Checking density file...")
    density_filename = 'dens_jb2008_032020_022224.mat'
    
    possible_paths = [
        density_filename,
        python_impl_dir / 'supporting_data' / density_filename,
        repo_root / 'supporting_data' / density_filename,
        repo_root / 'supporting_data' / 'TLEhistoric' / density_filename,
        script_dir / '..' / '..' / 'supporting_data' / density_filename,
        script_dir / '..' / '..' / '..' / 'supporting_data' / density_filename
    ]
    
    density_found = False
    for path in possible_paths:
        path = Path(path)
        if path.exists():
            print(f"✓ Density file found at: {path.resolve()}")
            density_found = True
            break
        else:
            print(f"✗ Not found: {path}")
    
    if not density_found:
        print("⚠️  Warning: Density file not found in any expected location")
    
    # Test 2: Check TLE data files
    print("\nTest 2: Checking TLE data files...")
    tle_years = [2018, 2019, 2020, 2021, 2022, 2023]
    tle_found_count = 0
    
    for year in tle_years:
        tle_filename = f"{year}.mat"
        tle_paths = [
            python_impl_dir / 'supporting_data' / 'TLEhistoric' / tle_filename,
            repo_root / 'supporting_data' / 'TLEhistoric' / tle_filename,
            script_dir / '..' / '..' / 'supporting_data' / 'TLEhistoric' / tle_filename
        ]
        
        year_found = False
        for path in tle_paths:
            path = Path(path)
            if path.exists():
                print(f"✓ TLE data {year} found at: {path.resolve()}")
                tle_found_count += 1
                year_found = True
                break
        
        if not year_found:
            print(f"✗ TLE data for {year} not found")
    
    print(f"Found TLE data for {tle_found_count}/{len(tle_years)} years")
    
    # Test 3: Check supporting functions
    print("\nTest 3: Checking supporting functions...")
    
    try:
        from mocat_mc import MOCATMC
        print("✓ Main MOCATMC class imported successfully")
    except ImportError as e:
        print(f"✗ Failed to import MOCATMC: {e}")
    
    # Test individual supporting functions
    supporting_functions = [
        'categorizeObj',
        'orbcontrol_vec', 
        'prop_mit_vec',
        'analytic_propagation_vec',
        'collision_prob_vec',
        'frag_col_SBM_vec',
        'cube_vec_v3'
    ]
    
    functions_found = 0
    for func_name in supporting_functions:
        try:
            module = __import__(f'supporting_functions.{func_name}', fromlist=[func_name])
            print(f"✓ {func_name} module imported successfully")
            functions_found += 1
        except ImportError as e:
            print(f"✗ Failed to import {func_name}: {e}")
    
    print(f"Successfully imported {functions_found}/{len(supporting_functions)} supporting functions")
    
    # Test 4: Check Python path configuration
    print("\nTest 4: Checking Python path configuration...")
    print("Current sys.path entries related to MOCAT-MC:")
    
    for i, path in enumerate(sys.path):
        if 'MOCAT-MC' in path or 'mocat' in path.lower():
            print(f"  {i}: {path}")
    
    # Test 5: Check write permissions
    print("\nTest 5: Checking write permissions...")
    
    test_dirs = [
        script_dir,
        python_impl_dir / 'examples',
        python_impl_dir / 'supporting_data'
    ]
    
    for test_dir in test_dirs:
        test_dir = Path(test_dir)
        try:
            test_file = test_dir / 'test_write_permissions.tmp'
            test_file.write_text('test')
            test_file.unlink()  # Delete test file
            print(f"✓ Write permissions OK in: {test_dir}")
        except Exception as e:
            print(f"✗ Write permission issue in {test_dir}: {e}")
    
    # Test 6: Test basic MOCAT-MC functionality
    print("\nTest 6: Testing basic MOCAT-MC functionality...")
    
    try:
        mocat = MOCATMC()
        print("✓ MOCATMC instance created successfully")
        
        # Test index mapping
        idx = mocat.idx
        expected_keys = ['a', 'ecco', 'inclo', 'mass', 'radius', 'objectclass']
        missing_keys = [key for key in expected_keys if key not in idx]
        
        if not missing_keys:
            print("✓ Index mapping contains all expected keys")
        else:
            print(f"✗ Missing index keys: {missing_keys}")
            
    except Exception as e:
        print(f"✗ MOCAT-MC functionality test failed: {e}")
    
    # Test 7: Memory and performance check
    print("\nTest 7: Basic memory and performance check...")
    
    try:
        # Create test arrays similar to what MOCAT-MC uses
        n_objects = 10000
        test_array = np.random.rand(n_objects, 24)  # Typical mat_sats size
        
        # Test basic operations
        altitudes = (test_array[:, 0] - 1) * 6378.137
        mean_alt = np.mean(altitudes)
        
        print(f"✓ Memory test passed: processed {n_objects} objects")
        print(f"  Test calculation result: mean altitude = {mean_alt:.1f} km")
        
    except Exception as e:
        print(f"✗ Memory/performance test failed: {e}")
    
    # Summary
    print("\n" + "="*50)
    print("TEST SUMMARY")
    print("="*50)
    
    if density_found:
        print("✓ Atmospheric density data: Available")
    else:
        print("⚠️  Atmospheric density data: Missing (simulations will use defaults)")
    
    if tle_found_count > 0:
        print(f"✓ TLE data: {tle_found_count} years available")
    else:
        print("✗ TLE data: None found")
    
    if functions_found >= len(supporting_functions) // 2:
        print("✓ Supporting functions: Most available")
    else:
        print("⚠️  Supporting functions: Many missing")
    
    print("\nRecommendations:")
    
    if not density_found:
        print("- Consider running simulations with simplified atmospheric model")
        
    if tle_found_count == 0:
        print("- Ensure TLE .mat files are in supporting_data/TLEhistoric/")
        
    if functions_found < len(supporting_functions):
        print("- Some supporting functions missing - simulations may have reduced functionality")
    
    print("\nTest complete!")


def test_specific_scenario():
    """Test a specific scenario to verify end-to-end functionality"""
    print("\n" + "="*50)
    print("SCENARIO-SPECIFIC TEST")
    print("="*50)
    
    try:
        from mocat_mc import MOCATMC
        
        print("Attempting to run basic MOCAT-MC setup...")
        mocat = MOCATMC()
        
        # Try to setup a basic configuration
        try:
            cfg = mocat.setup_mc_config(1, '2020.mat')
            print("✓ Basic configuration setup successful")
            print(f"  Initial population: {cfg['mat_sats'].shape[0]} objects")
            print(f"  Simulation time steps: {cfg['n_time']}")
            
        except FileNotFoundError as e:
            print(f"⚠️  Configuration setup failed (missing data): {e}")
            print("  This is expected if TLE data files are not available")
            
        except Exception as e:
            print(f"✗ Configuration setup failed: {e}")
            
    except ImportError as e:
        print(f"✗ Cannot run scenario test: {e}")


if __name__ == "__main__":
    # Change to script directory for relative path tests
    os.chdir(script_dir)
    
    # Run all tests
    test_file_paths()
    test_specific_scenario()
    
    print(f"\nTest script completed. Check output above for any issues.")
    print(f"Script location: {__file__}")
    print(f"Working directory: {os.getcwd()}")