#!/usr/bin/env python3
"""
Modify the summary analysis plot per user requirements
"""

import numpy as np
import matplotlib.pyplot as plt

def create_modified_summary_plot():
    """Create modified summary plot with requested changes"""
    
    fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, figsize=(16, 12))
    
    # Data
    scenarios = ['Basic\nPropagation', 'Collision\nTest', 'Atmospheric\nDrag', 'Full\nDefault']
    python_totals = [13709, 13672, 13673, 13558]
    matlab_totals = [13768, 13765, 13765, 13950]
    python_times = [1.67, 0.10, 0.09, 1.02]
    matlab_times = [58.09, 0.33, 0.33, 180.0]
    
    # Plot 1: Total objects comparison (NO PERCENTAGES)
    x = np.arange(len(scenarios))
    width = 0.35
    
    bars1 = ax1.bar(x - width/2, python_totals, width, label='Python', color='steelblue', alpha=0.8)
    bars2 = ax1.bar(x + width/2, matlab_totals, width, label='MATLAB', color='orange', alpha=0.8, hatch='///')
    
    ax1.set_title('Total Object Counts', fontweight='bold')
    ax1.set_ylabel('Total Objects')
    ax1.set_xticks(x)
    ax1.set_xticklabels(scenarios)
    ax1.legend()
    ax1.grid(True, alpha=0.3)
    
    # NO percentage annotations - removed as requested
    
    # Plot 2: Performance comparison (INCREASED Y-AXIS MAX)
    bars1 = ax2.bar(x - width/2, python_times, width, label='Python', color='green', alpha=0.8)
    bars2 = ax2.bar(x + width/2, matlab_times, width, label='MATLAB', color='red', alpha=0.8, hatch='///')
    
    ax2.set_title('Execution Time Comparison', fontweight='bold')
    ax2.set_ylabel('Time (seconds, log scale)')
    ax2.set_yscale('log')
    ax2.set_xticks(x)
    ax2.set_xticklabels(scenarios)
    ax2.legend()
    ax2.grid(True, alpha=0.3)
    
    # Increase y-axis max to 150% of current max
    current_max = max(matlab_times)
    new_max = current_max * 1.5
    ax2.set_ylim(bottom=0.01, top=new_max)
    
    # Add speedup factors
    for i, (py, mat) in enumerate(zip(python_times, matlab_times)):
        speedup = mat / py
        y_pos = max(py, mat) * 2
        ax2.annotate(f'{speedup:.1f}x', xy=(i, y_pos), ha='center', va='bottom',
                    color='green', fontweight='bold')
    
    # Plot 3: Satellite count accuracy (unchanged)
    python_sats = [1421, 1421, 1421, 1382]
    matlab_sats = [1421, 1421, 1421, 1382]
    
    ax3.scatter(python_sats, matlab_sats, s=150, alpha=0.7, color='purple')
    
    # Perfect agreement line
    min_sat = min(min(python_sats), min(matlab_sats))
    max_sat = max(max(python_sats), max(matlab_sats))
    ax3.plot([min_sat, max_sat], [min_sat, max_sat], 'k--', alpha=0.5, label='Perfect Agreement')
    
    ax3.set_xlabel('Python Satellite Count')
    ax3.set_ylabel('MATLAB Satellite Count')
    ax3.set_title('Satellite Count Accuracy\n(Perfect Agreement)', fontweight='bold')
    ax3.legend()
    ax3.grid(True, alpha=0.3)
    
    # Add scenario labels
    for i, scenario in enumerate(['T1', 'T2', 'T3', 'T4']):
        ax3.annotate(scenario, xy=(python_sats[i], matlab_sats[i]), 
                    xytext=(5, 5), textcoords='offset points', fontweight='bold')
    
    # Plot 4: Difference magnitude by object type (unchanged)
    object_types = ['Satellites', 'Derelicts', 'Debris', 'Rocket Bodies']
    
    # Average absolute differences across all tests
    sat_diffs = [0, 0, 0, 0]  # Satellites match exactly
    der_diffs = [23, 41, 40, 44]  # Derelict differences
    deb_diffs = [25, 37, 37, 61]  # Debris differences  
    rb_diffs = [11, 15, 15, 18]   # Rocket body differences
    
    avg_diffs = [np.mean(sat_diffs), np.mean(der_diffs), np.mean(deb_diffs), np.mean(rb_diffs)]
    colors_obj = ['#2E8B57', '#CD853F', '#DC143C', '#4169E1']
    
    bars = ax4.bar(object_types, avg_diffs, color=colors_obj, alpha=0.8)
    ax4.set_title('Average Absolute Differences\nby Object Type', fontweight='bold')
    ax4.set_ylabel('Average |Python - MATLAB|')
    ax4.set_xticks(range(len(object_types)))
    ax4.set_xticklabels(object_types, rotation=45, ha='right')
    ax4.grid(True, alpha=0.3)
    
    # Add value labels
    for bar, diff in zip(bars, avg_diffs):
        height = bar.get_height()
        ax4.annotate(f'{diff:.0f}', xy=(bar.get_x() + bar.get_width()/2, height),
                    xytext=(0, 3), textcoords="offset points", ha='center', va='bottom',
                    fontweight='bold')
    
    plt.suptitle('MOCAT-MC Python vs MATLAB: Summary Analysis\n' +
                'Key Performance and Accuracy Metrics',
                fontsize=16, fontweight='bold')
    
    plt.tight_layout()
    plt.savefig('mocat_mc_summary_analysis.png', dpi=300, bbox_inches='tight')
    
    print("Modified summary analysis plot created:")
    print("- Removed percentages from top-left Total Object Counts plot")
    print("- Increased y-axis max to 150% for top-right Execution Time plot")
    print("- All other elements unchanged")
    
    return fig

if __name__ == "__main__":
    create_modified_summary_plot()
    print("Plot modification completed!")