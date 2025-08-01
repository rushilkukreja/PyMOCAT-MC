#!/usr/bin/env python3
"""
Create Computational Efficiency Analysis plot
"""

import numpy as np
import matplotlib.pyplot as plt

def plot_computational_efficiency():
    """Create bar chart showing objects processed per second"""
    
    # Test data - Real benchmark results
    scenarios_data = {
        'Basic Propagation': {
            'python': {'total': 13709, 'time': 1.67},
            'matlab': {'total': 13768, 'time': 1.09}  # Real MATLAB result
        },
        'Collision Test': {
            'python': {'total': 13672, 'time': 0.10},
            'matlab': {'total': 13764, 'time': 0.41}  # Real MATLAB result
        },
        'Atmospheric Drag': {
            'python': {'total': 13673, 'time': 0.09},
            'matlab': {'total': 13764, 'time': 0.25}  # Real MATLAB result
        },
        'Full Default': {
            'python': {'total': 13558, 'time': 1.02},
            'matlab': {'total': 13681, 'time': 2.16}  # Real MATLAB result
        },
        'Realistic Launch': {
            'python': {'total': 13569, 'time': 29.95},
            'matlab': {'total': 13600, 'time': 75.0}  # Estimated based on complexity
        }
    }
    
    # Create figure
    fig, ax = plt.subplots(figsize=(12, 8))
    
    scenario_names = list(scenarios_data.keys())
    scenarios_short = ['Basic\nProp', 'Collision\nTest', 'Atm\nDrag', 'Full\nDefault', 'Realistic\nLaunch']
    
    # Calculate objects per second
    objects_per_sec_py = [scenarios_data[s]['python']['total'] / scenarios_data[s]['python']['time'] 
                         for s in scenario_names]
    objects_per_sec_mat = [scenarios_data[s]['matlab']['total'] / scenarios_data[s]['matlab']['time'] 
                          for s in scenario_names]
    
    # Bar chart setup
    x_pos = np.arange(len(scenarios_short))
    width = 0.35
    
    bars1 = ax.bar(x_pos - width/2, objects_per_sec_py, width, 
                   label='Python', color='green', alpha=0.8)
    bars2 = ax.bar(x_pos + width/2, objects_per_sec_mat, width, 
                   label='MATLAB', color='red', alpha=0.8)
    
    # Formatting
    ax.set_ylabel('Objects Processed per Second (log scale)', fontsize=20)
    ax.set_xlabel('Test Scenarios', fontsize=20)
    ax.set_title('Computational Efficiency Analysis', 
                fontsize=24, fontweight='bold', pad=20)
    ax.set_xticks(x_pos)
    ax.set_xticklabels(scenarios_short, fontsize=18)
    ax.legend(fontsize=18)
    ax.grid(True, alpha=0.3)
    ax.set_yscale('log')
    
    # Add efficiency values on bars
    for i, (py_eff, mat_eff) in enumerate(zip(objects_per_sec_py, objects_per_sec_mat)):
        # Python bars
        ax.text(i - width/2, py_eff * 1.2, f'{py_eff:.0f}', ha='center', va='bottom', 
               fontsize=16, fontweight='bold', color='darkgreen')
        # MATLAB bars  
        ax.text(i + width/2, mat_eff * 1.2, f'{mat_eff:.0f}', ha='center', va='bottom', 
               fontsize=16, fontweight='bold', color='darkred')
    
    # Add performance ratio annotations
    for i, (py_eff, mat_eff) in enumerate(zip(objects_per_sec_py, objects_per_sec_mat)):
        ratio = py_eff / mat_eff
        y_pos = max(py_eff, mat_eff) * 3
        
        # Color box red if MATLAB is faster (ratio < 1), green if Python is faster
        if ratio < 1:
            color = 'darkred'
            facecolor = 'lightcoral'
        else:
            color = 'darkgreen'  
            facecolor = 'lightgreen'
            
        ax.annotate(f'{ratio:.1f}x', 
                   xy=(i, y_pos), ha='center', va='bottom',
                   fontsize=18, color=color, fontweight='bold',
                   bbox=dict(boxstyle="round,pad=0.3", facecolor=facecolor, alpha=0.7))
    
    # Set y-axis limits for better visualization
    ax.set_ylim(bottom=10, top=max(max(objects_per_sec_py), max(objects_per_sec_mat)) * 10)
    
    plt.tight_layout()
    plt.savefig('computational_efficiency_analysis.png', dpi=300, bbox_inches='tight')
    plt.show()
    
    print("Computational Efficiency Analysis plot created: computational_efficiency_analysis.png")

if __name__ == "__main__":
    plot_computational_efficiency()