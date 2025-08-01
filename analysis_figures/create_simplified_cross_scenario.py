#!/usr/bin/env python3
"""
Create simplified cross-scenario analysis with only 2 requested plots
"""

import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

def create_simplified_cross_scenario():
    """Create simplified cross-scenario analysis with only 2 plots"""
    
    # Test data
    scenarios_data = {
        'Basic Propagation': {
            'python': {'S': 1421, 'D': 1843, 'N': 9421, 'B': 1024, 'total': 13709, 'time': 1.67},
            'matlab': {'S': 1421, 'D': 1866, 'N': 9446, 'B': 1035, 'total': 13768, 'time': 58.09}
        },
        'Collision Test': {
            'python': {'S': 1421, 'D': 1825, 'N': 9407, 'B': 1019, 'total': 13672, 'time': 0.10},
            'matlab': {'S': 1421, 'D': 1866, 'N': 9444, 'B': 1034, 'total': 13765, 'time': 0.33}
        },
        'Atmospheric Drag': {
            'python': {'S': 1421, 'D': 1826, 'N': 9407, 'B': 1019, 'total': 13673, 'time': 0.09},
            'matlab': {'S': 1421, 'D': 1866, 'N': 9444, 'B': 1034, 'total': 13765, 'time': 0.33}
        },
        'Full Default': {
            'python': {'S': 1382, 'D': 1805, 'N': 9363, 'B': 1008, 'total': 13558, 'time': 1.02},
            'matlab': {'S': 1382, 'D': 1849, 'N': 9424, 'B': 1026, 'total': 13950, 'time': 180.0}
        }
    }
    
    # Create figure with 2 subplots
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(16, 8))
    
    scenario_names = list(scenarios_data.keys())
    scenarios_short = ['Basic\nProp', 'Collision\nTest', 'Atm\nDrag', 'Full\nDefault']
    
    # Plot 1: Object Type Percentage Difference Heatmap
    object_types = ['S', 'D', 'N', 'B']
    pct_diff_matrix = []
    
    for obj_type in object_types:
        row = []
        for scenario in scenario_names:
            py_count = scenarios_data[scenario]['python'][obj_type]
            mat_count = scenarios_data[scenario]['matlab'][obj_type]
            pct_diff = 100 * (py_count - mat_count) / mat_count if mat_count > 0 else 0
            row.append(pct_diff)
        pct_diff_matrix.append(row)
    
    pct_diff_matrix = np.array(pct_diff_matrix)
    
    # Create heatmap
    im = ax1.imshow(pct_diff_matrix, cmap='RdBu_r', aspect='auto', vmin=-3, vmax=3)
    ax1.set_xticks(range(len(scenario_names)))
    ax1.set_xticklabels(scenarios_short)
    ax1.set_yticks(range(len(object_types)))
    ax1.set_yticklabels(['Satellites', 'Derelicts', 'Debris', 'Rocket Bodies'])
    ax1.set_title('Object Type Percentage Difference Heatmap', fontweight='bold')
    
    # Add text annotations
    for i in range(len(object_types)):
        for j in range(len(scenario_names)):
            text_color = 'white' if abs(pct_diff_matrix[i, j]) > 1.5 else 'black'
            ax1.text(j, i, f'{pct_diff_matrix[i, j]:+.1f}%',
                    ha="center", va="center", color=text_color, fontweight='bold')
    
    # Add colorbar
    cbar1 = plt.colorbar(im, ax=ax1, label='Percentage Difference (%)')
    
    # Plot 2: Computational Efficiency Analysis
    objects_per_sec_py = [scenarios_data[s]['python']['total'] / scenarios_data[s]['python']['time'] for s in scenario_names]
    objects_per_sec_mat = [scenarios_data[s]['matlab']['total'] / scenarios_data[s]['matlab']['time'] for s in scenario_names]
    
    x_pos = np.arange(len(scenarios_short))
    width = 0.35
    
    bars1 = ax2.bar(x_pos - width/2, objects_per_sec_py, width, label='Python', color='green', alpha=0.8)
    bars2 = ax2.bar(x_pos + width/2, objects_per_sec_mat, width, label='MATLAB', color='red', alpha=0.8)
    
    ax2.set_ylabel('Objects Processed per Second')
    ax2.set_xlabel('Test Scenarios')
    ax2.set_title('Computational Efficiency (Objects Processed per Second)', fontweight='bold')
    ax2.set_xticks(x_pos)
    ax2.set_xticklabels(scenarios_short)
    ax2.legend()
    ax2.grid(True, alpha=0.3)
    ax2.set_yscale('log')
    
    # Add efficiency values on bars
    for i, (py_eff, mat_eff) in enumerate(zip(objects_per_sec_py, objects_per_sec_mat)):
        # Python bars
        ax2.text(i - width/2, py_eff * 1.2, f'{py_eff:.0f}', ha='center', va='bottom', 
                fontsize=9, fontweight='bold', color='darkgreen')
        # MATLAB bars  
        ax2.text(i + width/2, mat_eff * 1.2, f'{mat_eff:.0f}', ha='center', va='bottom', 
                fontsize=9, fontweight='bold', color='darkred')
    
    # Overall title
    fig.suptitle('MOCAT-MC Cross-Scenario Analysis: Python vs MATLAB\n' +
                'Object Type Accuracy and Computational Performance',
                fontsize=16, fontweight='bold')
    
    
    plt.tight_layout()
    plt.savefig('mocat_mc_cross_scenario_analysis.png', dpi=300, bbox_inches='tight')
    
    print("Simplified cross-scenario analysis created:")
    print("- Object Type Percentage Difference Heatmap (left)")
    print("- Computational Efficiency Analysis (right)")
    
    return fig

if __name__ == "__main__":
    create_simplified_cross_scenario()
    print("Simplified plot creation completed!")