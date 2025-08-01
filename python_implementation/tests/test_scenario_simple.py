#!/usr/bin/env python3
"""
Simple test of scenario_no_launch functionality
"""

import numpy as np
import matplotlib.pyplot as plt
import sys
import os

# Add parent directory to path to import mocat_mc
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from mocat_mc import MOCATMC


def test_simple():
    print("=== Simple Test ===")
    
    # Initialize MOCAT-MC
    mocat = MOCATMC()
    
    # Initial condition file
    ic_file = '2020.mat'
    seed = 1
    
    print('MC configuration starting...')
    cfg_mc = mocat.setup_mc_config(seed, ic_file)
    
    # Just use 1-year simulation for testing
    print(f'Configuration n_time: {cfg_mc["n_time"]}')
    print(f'Initial Population: {cfg_mc["mat_sats"].shape[0]} sats')
    
    # Disable collisions for simplicity
    cfg_mc['skipCollisions'] = 1
    
    print('Starting main_mc...')
    nS, nD, nN, nB, mat_sats = mocat.main_mc(cfg_mc, seed)
    
    print(f'Simulation completed: {nS} satellites, {nD} derelicts, {nN} debris, {nB} rocket bodies')
    
    # Test the deorbit tracking method
    print('Testing deorbit tracking...')
    nS2, nD2, nN2, nB2, mat_sats2, deorbitlist_r = mocat.main_mc_with_deorbit_tracking(cfg_mc, seed)
    
    print(f'Deorbit tracking completed: {nS2} satellites, {nD2} derelicts, {nN2} debris, {nB2} rocket bodies')
    print(f'Deorbit list shape: {deorbitlist_r.shape}')
    print(f'Final cumulative deorbited: {deorbitlist_r[-1]:.0f}')
    
    # Create simple plot
    plt.figure(figsize=(10, 6))
    time_steps = np.arange(len(deorbitlist_r))
    plt.plot(time_steps, deorbitlist_r, 'b-', linewidth=2, label='Cumulative Deorbited')
    plt.xlabel('Time Steps')
    plt.ylabel('Cumulative Deorbited Objects')
    plt.title('Simple Deorbit Test')
    plt.legend()
    plt.grid(True, alpha=0.3)
    plt.tight_layout()
    plt.savefig('python_test_simple_deorbit.png', dpi=150, bbox_inches='tight')
    print('Saved: python_test_simple_deorbit.png')
    plt.show()
    
    print('Simple test completed successfully!')


if __name__ == "__main__":
    test_simple()