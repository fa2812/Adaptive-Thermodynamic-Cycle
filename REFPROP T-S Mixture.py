# -*- coding: utf-8 -*-
"""
Created on Sat Oct 30 14:20:49 2021

@author: fata1r18
"""
# pip install ctREFPROP

import os,numpy as np
from matplotlib import pyplot as plt
from ctREFPROP.ctREFPROP import REFPROPFunctionLibrary

os.environ['RPPREFIX'] = r'H:\IP 202122\REFPROP 10'

RP = REFPROPFunctionLibrary(os.environ['RPPREFIX']) #Initialise REFPROP
RP.SETPATHdll(os.environ['RPPREFIX'])
MASS_BASE_SI = RP.GETENUMdll(0,"MASS BASE SI").iEnum

#Fluid properties

fluid = 'ISOBUTAN * IPENTANE'
comp = [0.5,0.5]
T_cond = 273 + 30
T_evap = 273 + 150
eff_pump = 1
eff_turb = 1

r_crit = RP.REFPROPdll(fluid,"","TC",MASS_BASE_SI, 0,0,0,0, comp)
T_crit = r_crit.Output[0]

# State calculations
r1 = RP.REFPROPdll(fluid,"TQ","P;H;S",MASS_BASE_SI, 0,0,T_cond,0, comp)
P_1 = r1.Output[0]
H_1 = r1.Output[1]
S_1 = r1.Output[2]

r3 = RP.REFPROPdll(fluid,"TQ","P;H;S",MASS_BASE_SI, 0,0,T_evap,1, comp)
P_3 = r3.Output[0]
H_3 = r3.Output[1]
S_3 = r3.Output[2]

r22 = RP.REFPROPdll(fluid,"PQ","T;S",MASS_BASE_SI, 0,0,P_3,0, comp)
T_22 = r22.Output[0]
S_22 = r22.Output[1]

r2s = RP.REFPROPdll(fluid,"SP","T;H",MASS_BASE_SI, 0,0,S_1,P_3, comp)
T_2s = r2s.Output[0]
H_2s = r2s.Output[1]

r44 = RP.REFPROPdll(fluid,"PQ","T;S",MASS_BASE_SI, 0,0,P_1,1, comp)
T_44 = r44.Output[0]
S_44 = r44.Output[1]

r4s = RP.REFPROPdll(fluid,"SP","T;H",MASS_BASE_SI, 0,0,S_3,P_1, comp)
T_4s = r4s.Output[0]
H_4s = r4s.Output[1]

H_2 = (H_2s - H_1)/eff_pump + H_1
S_2 = S_1/eff_pump
r2 = RP.REFPROPdll(fluid,"HS","T",MASS_BASE_SI, 0,0,H_2,S_2, comp)
T_2 = r2.Output[0]

H_4 = H_3 - eff_turb * (H_3 - H_4s)
S_4 = S_3/eff_turb
r4 = RP.REFPROPdll(fluid,"HS","T",MASS_BASE_SI, 0,0,H_4,S_4, comp)
T_4 = r4.Output[0]

#print(P_1,H_1,S_1)
#print(T_2s, H_2s)
#print(P_3, H_3, S_3)
#print(T_4s, H_4s)

eff_TH = 1 - (H_4 - H_1)/(H_3 - H_2)
print('Cycle efficiency:', eff_TH)

# Plotting Graph

# State Points

S_states = [S_1, S_2, S_22, S_3, S_4, S_44, S_1]
T_states = [T_cond, T_2, T_22, T_evap, T_4, T_44, T_cond]
#print (S_states, T_states)

# Saturation Line
r0 = RP.REFPROPdll(fluid,"SQ","P;H;T",MASS_BASE_SI, 0,0,0,0, comp)
T_s0 = r0.Output[2]
S_sat = []
Temp_range_1 = np.linspace(T_s0,T_crit,500)
Temp_range_2 = np.linspace(T_crit,T_s0,500)
Temp_graph = list(Temp_range_1)

for T_1 in Temp_range_1:
    R = RP.REFPROPdll(fluid,"TQ","P;H;S",MASS_BASE_SI, 0,0,T_1,0, comp)
    S_sat.append(R.Output[2])

for T_2 in Temp_range_2:
    R = RP.REFPROPdll(fluid,"TQ","P;H;S",MASS_BASE_SI, 0,0,T_2,1, comp)
    S_sat.append(R.Output[2])   
    Temp_graph.append(T_2)

plt.title('T-S Diagram of ORC (' + fluid + ')')
for i in range(6):
    if i < 2:
        plt.annotate(str(i+1),(S_states[i], T_states[i]))
    elif i == 2:
        plt.annotate('22',(S_states[i], T_states[i]))
    elif 2 < i < 5:
        plt.annotate(str(i),(S_states[i], T_states[i]))
    elif i == 5:
        plt.annotate('44',(S_states[i], T_states[i]))
plt.xlabel('s (J / kg K)')
plt.ylabel('T (K)')
plt.plot(S_sat,Temp_graph)
plt.plot(S_states, T_states, 'o-')
#print(S_ent)