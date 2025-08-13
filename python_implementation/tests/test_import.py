#!/usr/bin/env python3
"""Test basic imports and setup"""

print("Testing Python MOCAT-MC imports...")

try:
    import numpy as np
    print("✓ NumPy imported successfully")
    
    import scipy
    print("✓ SciPy imported successfully")
    
    import sys
    import os
    # Add parent directory to path for mocat_mc import
    sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..'))
    
    from mocat_mc import MOCATMC
    print("✓ MOCATMC imported successfully")
    
    # Try to initialize
    mocat = MOCATMC()
    print("✓ MOCATMC initialized")
    
    # Check if we can load the initial conditions
    import scipy.io as sio
    from pathlib import Path
    
    # Find the 2020.mat file using multiple possible paths
    script_dir = Path(__file__).parent
    possible_paths = [
        script_dir / '..' / 'supporting_data' / 'TLEhistoric' / '2020.mat',
        script_dir / '..' / '..' / 'python_implementation' / 'supporting_data' / 'TLEhistoric' / '2020.mat',
        Path('python_implementation') / 'supporting_data' / 'TLEhistoric' / '2020.mat',
        Path('supporting_data') / 'TLEhistoric' / '2020.mat'
    ]
    
    data_found = False
    for ic_path in possible_paths:
        if ic_path.exists():
            try:
                data = sio.loadmat(str(ic_path))
                if 'mat_sats' in data:
                    print(f"✓ Found 2020.mat with {data['mat_sats'].shape[0]} satellites at {ic_path}")
                    data_found = True
                    break
                else:
                    print(f"⚠ Found {ic_path} but no mat_sats data")
            except Exception as e:
                print(f"⚠ Found {ic_path} but couldn't load: {e}")
    
    if not data_found:
        print("⚠ Could not find 2020.mat file - tests will use synthetic data")
    
    print("\nPython setup appears to be working!")
    
except Exception as e:
    print(f"✗ Error: {e}")
    import traceback
    traceback.print_exc()