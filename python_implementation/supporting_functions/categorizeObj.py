"""
Separate satellites into satellite, derelict, nonusable/debris, rocket body

Categorizes space objects based on their object class and control status.

Args:
    objint_st: array of object indices (object class)
    cont_st: array of control indices (controlled=1, uncontrolled=0)

Returns:
    nS: number of satellites (object class 1, controlled)
    nD: number of derelicts (object class 1, uncontrolled)  
    nN: number of debris objects (object classes 3,4,6,7,8,9+)
    nB: number of rocket bodies (object class 5)

Object classes:
  P(payload), PMRO(debris), Pfrag(debris), Pdeb(debris), RB(rocket body),
  RBMRO(debris), RBfrag(debris), RBdeb(debris), Deb(debris), OtherDeb(debris),
  Unkwn(debris), untracked(debris)

Categories:
  S: satellite, object index 1 and control index 1
  D: derelict, object index 1 and control index 0
  N: nonusable/debris, object index 3,4,6,7,8,9..
  B: rocket body, object index 5
"""

import numpy as np


def categorizeObj(objint_st, cont_st):
    """
    Separate satellites into satellite, derelict, nonusable/debris, rocket body
    
    Args:
        objint_st: array of object indices (object class)
        cont_st: array of control indices (controlled=1, uncontrolled=0)
        
    Returns:
        nS: number of satellites (object class 1, controlled)
        nD: number of derelicts (object class 1, uncontrolled)
        nN: number of debris objects (object classes 3,4,6,7,8,9+)
        nB: number of rocket bodies (object class 5)
    """
    
    # Convert to numpy arrays for vectorized operations
    objint_st = np.asarray(objint_st)
    cont_st = np.asarray(cont_st)
    
    # Count satellites: object class 1 and controlled
    nS = np.sum((objint_st == 1) & (cont_st == 1))
    
    # Count derelicts: object class 1 and uncontrolled
    nD = np.sum((objint_st == 1) & (cont_st == 0))
    
    # Count debris: object classes 3,4,6,7,8,9+ (excluding 1,2,5)
    nN = np.sum((objint_st == 3) | (objint_st == 4) | (objint_st >= 6))
    
    # Count rocket bodies: object class 5
    nB = np.sum(objint_st == 5)
    
    return nS, nD, nN, nB