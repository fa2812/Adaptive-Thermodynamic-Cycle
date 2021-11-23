# -*- coding: utf-8 -*-
"""
Created on Mon Oct 18 16:38:35 2021

@author: fata1r18
"""

import os,numpy as np
from matplotlib import pyplot as plt
from ctREFPROP.ctREFPROP import REFPROPFunctionLibrary

os.environ['RPPREFIX'] = r'H:\IP 202122\REFPROP 10'

RP = REFPROPFunctionLibrary(os.environ['RPPREFIX']) #Initialise REFPROP
RP.SETPATHdll(os.environ['RPPREFIX'])
MASS_BASE_SI = RP.GETENUMdll(0,"MASS BASE SI").iEnum

def efficiency(fluid,T_cond,T_evap,eff_pump,eff_turb):
    
    # State calculations
    r1 = RP.REFPROPdll(fluid,"TQ","P;H;S",MASS_BASE_SI, 0,0,T_cond,0, [1.0])
    P_1 = r1.Output[0]
    H_1 = r1.Output[1]
    S_1 = r1.Output[2]
    
    r3 = RP.REFPROPdll(fluid,"TQ","P;H;S",MASS_BASE_SI, 0,0,T_evap,1, [1.0])
    P_3 = r3.Output[0]
    H_3 = r3.Output[1]
    S_3 = r3.Output[2]
    
    r2s = RP.REFPROPdll(fluid,"SP","T;H",MASS_BASE_SI, 0,0,S_1,P_3, [1.0])
    H_2s = r2s.Output[1]
    
    r4s = RP.REFPROPdll(fluid,"SP","T;H",MASS_BASE_SI, 0,0,S_3,P_1, [1.0])
    H_4s = r4s.Output[1]
    
    H_2 = (H_2s - H_1)/eff_pump + H_1
    
    H_4 = H_3 - eff_turb * (H_3 - H_4s)
    
    eff_TH = 1 - (H_4 - H_1)/(H_3 - H_2)
    return eff_TH

T_H_range = [320, 340, 360, 380, 400, 423]
T_H_2 = [x+5 for x in T_H_range]
T_H_2[2] = 362
T_H_2[5] = 420
eff_Toluene = []
eff_R123 = []

for T in T_H_range:    
    result_Tol = efficiency('TOLUENE',303,T,1,1)
    eff_Toluene.append(result_Tol)
    
for T in T_H_2:
    result_R123 = efficiency('R123',303,T,1,1)
    eff_R123.append(result_R123)

plt.title('Variation in Cycle Efficiency as a function of T_H')
plt.xlabel('T_H (K)')
plt.ylabel('Cycle Efficiency')
plt.plot(T_H_range,eff_Toluene,'vk', label='Toluene')
plt.plot(T_H_2,eff_R123,'^k', label='R123')
plt.legend()
plt.xlim([280, 440])
plt.yticks(np.arange(0.00,0.24,0.04))