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
    sys.path.append(os.path.dirname(os.path.abspath(__file__)))
    
    from mocat_mc import MOCATMC
    print("✓ MOCATMC imported successfully")
    
    # Try to initialize
    mocat = MOCATMC()
    print("✓ MOCATMC initialized")
    
    # Check if we can load the initial conditions
    import scipy.io as sio
    ic_path = os.path.join('supporting_data', 'TLEhistoric', '2020.mat')
    if os.path.exists(ic_path):
        data = sio.loadmat(ic_path)
        if 'mat_sats' in data:
            print(f"✓ Found 2020.mat with {data['mat_sats'].shape[0]} satellites")
        else:
            print("✗ mat_sats not found in 2020.mat")
    else:
        print(f"✗ Could not find {ic_path}")
    
    print("\nPython setup appears to be working!")
    
except Exception as e:
    print(f"✗ Error: {e}")
    import traceback
    traceback.print_exc()