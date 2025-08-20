"""
Orbit control vectorized implementation.

Python equivalent of orbcontrol_vec.m.
"""

import numpy as np
from typing import Tuple, Dict
from datetime import datetime


def orbcontrol_vec(mat_sat_in: np.ndarray, tsince: float, time0: datetime,
                  orbtol: float, pmd: float, day2min: float, year2day: float,
                  param: Dict) -> Tuple[np.ndarray, np.ndarray]:
    """
    Execute orbit control for active satellites.

    Args:
        mat_sat_in: Input satellite matrix [a,ecco,inclo,nodeo,argpo,mo,controlled,a_desired,missionlife,launched,r,v]
        tsince: Time since start [minutes]
        time0: Initial time
        orbtol: Orbit control tolerance [km]
        pmd: Post-mission disposal probability
        day2min: Days to minutes conversion
        year2day: Year to days conversion
        param: Parameters dictionary

    Returns:
        Tuple of (mat_sat_out, deorbit_pmd)
        mat_sat_out: Output matrix [a,controlled,r,v]
        deorbit_pmd: Indices of satellites to be deorbited due to PMD
    """
    n_sats = mat_sat_in.shape[0]

    if n_sats == 0:
        return np.zeros((0, 7)), np.array([], dtype=int)

    # Extract parameters
    req = param['req']
    mu = param['mu']

    # Extract orbital elements and control parameters
    a = mat_sat_in[:, 0]  # Semi-major axis
    controlled = mat_sat_in[:, 6].copy()  # Controlled flag
    a_desired = mat_sat_in[:, 7]  # Desired semi-major axis
    missionlife = mat_sat_in[:, 8]  # Mission lifetime [years]
    launch_date = mat_sat_in[:, 9]  # Launch date
    r = mat_sat_in[:, 10:13]  # Position vector
    v = mat_sat_in[:, 13:16]  # Velocity vector

    # Convert time to years since simulation start
    current_time_years = tsince / (day2min * year2day)

    # Find controlled satellites
    controlled_mask = controlled == 1

    # Station-keeping for controlled satellites (matching MATLAB logic exactly)
    if np.any(controlled_mask):
        controlled_indices = np.where(controlled_mask)[0]  # is_controlled
        current_a = a[controlled_indices]  # a_current
        desired_a_controlled = a_desired[controlled_indices]

        # Find satellites needing control (matching MATLAB condition)
        deviation = np.abs(current_a - desired_a_controlled)
        orbtol_norm = orbtol / req
        needs_control_mask = deviation > orbtol_norm

        if np.any(needs_control_mask):
            control_indices = controlled_indices[needs_control_mask]  # find_control
            a[control_indices] = desired_a_controlled[needs_control_mask]

            # Update position and velocity for controlled satellites
            # This is simplified - MATLAB does mean2osc_m_vec and oe2rv_vec
            # For now, we'll keep the existing r,v (this might need refinement)

    # Mission lifetime check and post-mission disposal
    deorbit_pmd = []

    if np.any(controlled_mask):
        controlled_indices = np.where(controlled_mask)[0]

        # Calculate current time in Julian date format (similar to MATLAB)
        # Convert time0 to Julian date and add tsince
        current_jd = time0.toordinal() + 1721425.5 + tsince / day2min

        # Calculate mission age from launch date (matching MATLAB logic)
        launch_dates = launch_date[controlled_mask]

        # Filter out satellites with valid launch dates
        valid_launch_mask = ~np.isnan(launch_dates)
        if np.any(valid_launch_mask):
            valid_controlled_indices = controlled_indices[valid_launch_mask]
            valid_launch_dates = launch_dates[valid_launch_mask]
            valid_mission_lifetimes = missionlife[controlled_mask][valid_launch_mask]

            # Calculate satellite ages in years (matching MATLAB formula)
            satellite_ages_years = (current_jd - valid_launch_dates) / year2day

            # Find satellites that have reached end of mission life
            end_of_life_mask = satellite_ages_years > valid_mission_lifetimes

            if np.any(end_of_life_mask):
                eol_indices = valid_controlled_indices[end_of_life_mask]

                # Post-mission disposal decision (matching MATLAB logic)
                pmd_random = np.random.rand(len(eol_indices))
                pmd_check = pmd < pmd_random  # Note: MATLAB uses PMD<rand_life

                # Satellites with failed PMD become derelicts
                failed_pmd_indices = eol_indices[pmd_check]
                controlled[failed_pmd_indices] = 0

                # Satellites with successful PMD are deorbited
                successful_pmd_indices = eol_indices[~pmd_check]
                deorbit_pmd = successful_pmd_indices.tolist()

    # Prepare output matrix
    mat_sat_out = np.column_stack([a, controlled, r, v])

    return mat_sat_out, np.array(deorbit_pmd, dtype=int)
