# -*- coding: utf-8 -*-
"""
Created on Wed Oct 27 15:02:30 2021

@author: fata1r18
"""

import os,numpy as np
from matplotlib import pyplot as plt
from ctREFPROP.ctREFPROP import REFPROPFunctionLibrary

os.environ['RPPREFIX'] = r'H:\IP 202122\REFPROP 10'

RP = REFPROPFunctionLibrary(os.environ['RPPREFIX']) #Initialise REFPROP
RP.SETPATHdll(os.environ['RPPREFIX'])
MASS_BASE_SI = RP.GETENUMdll(0,"MASS BASE SI").iEnum

def Rankine(fluid,T_cond,T_evap,eff_pump,eff_turb,comp,prop):
    
    # State calculations
    r1 = RP.REFPROPdll(fluid,"TQ","P;H;S",MASS_BASE_SI, 0,0,T_cond,0, comp)
    P_1 = r1.Output[0]
    H_1 = r1.Output[1]
    S_1 = r1.Output[2]
    
    r3 = RP.REFPROPdll(fluid,"TQ","P;H;S",MASS_BASE_SI, 0,0,T_evap,1, comp)
    P_3 = r3.Output[0]
    H_3 = r3.Output[1]
    S_3 = r3.Output[2]
    
    r2s = RP.REFPROPdll(fluid,"SP","T;H",MASS_BASE_SI, 0,0,S_1,P_3, comp)
    H_2s = r2s.Output[1]
    
    r4s = RP.REFPROPdll(fluid,"SP","T;H",MASS_BASE_SI, 0,0,S_3,P_1, comp)
    H_4s = r4s.Output[1]
    
    H_2 = (H_2s - H_1)/eff_pump + H_1
    r2 = RP.REFPROPdll(fluid,"HP","T;S",MASS_BASE_SI, 0,0,H_2,P_3, comp)
    T_2 = r2.Output[0]
    S_2 = r2.Output[1]
    
    H_4 = H_3 - eff_turb * (H_3 - H_4s)
    r4 = RP.REFPROPdll(fluid,"HP","T;S",MASS_BASE_SI, 0,0,H_4,P_1, comp)
    T_4 = r2.Output[0]
    S_4 = r2.Output[1]
    
    if prop == "H":
        return H_1, H_2, H_3, H_4
    elif prop == "P":
        return P_1, P_3, P_3, P_1
    elif prop == "S":
        return S_1, S_2, S_3, S_4
    elif prop == "T":
        return T_cond, T_2, T_evap, T_4

def Rankine_eff(H_1,H_2,H_3,H_4):
    eff_TH = 1 - (H_4 - H_1)/(H_3 - H_2)
    return eff_TH

def mixture(fluid1,fluid2,comp1,comp2,T_cond,T_evap,eff_pump,eff_turb):
    fluid_mix = fluid1 + " * " + fluid2
    comp_mix = [comp1, comp2]
    H_1, H_2, H_3, H_4 = Rankine(fluid_mix,T_cond,T_evap,eff_pump,eff_turb,comp_mix,"H")
    efficiency = Rankine_eff(H_1, H_2, H_3, H_4)
    print(efficiency)
    
mixture("ISOBUTAN","IPENTANE",0.4,0.6,273+30,273+150,1,1)