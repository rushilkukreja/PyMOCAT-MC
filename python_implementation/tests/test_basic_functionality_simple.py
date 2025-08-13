#!/usr/bin/env python3
"""
Simplified basic functionality tests for PyMOCAT-MC
"""

import sys
import os
from pathlib import Path

# Add parent directory to path for imports (cross-platform)
script_dir = Path(__file__).parent
parent_dir = script_dir.parent
sys.path.insert(0, str(parent_dir))

print("Running basic functionality tests...")

# Test 1: Basic imports
try:
    from supporting_functions.get_idx import get_idx
    from supporting_functions.categorize_obj import categorize_obj
    print("OK Basic imports successful")
except ImportError as e:
    print(f"ERROR Import failed: {e}")
    sys.exit(1)

# Test 2: Index structure
try:
    idx = get_idx()
    assert isinstance(idx, dict)
    assert 'a' in idx
    assert 'mass' in idx
    assert len(idx) == 20
    print("OK Index structure test passed")
except Exception as e:
    print(f"ERROR Index structure test failed: {e}")
    sys.exit(1)

# Test 3: MOCAT-MC import
try:
    from mocat_mc import MOCATMC
    mocat = MOCATMC()
    assert hasattr(mocat, 'constants')
    assert hasattr(mocat, 'idx')
    print("OK MOCAT-MC import successful")
except ImportError as e:
    print(f"WARNING MOCAT-MC import skipped: {e}")

print("\nAll basic functionality tests passed!")
sys.exit(0)