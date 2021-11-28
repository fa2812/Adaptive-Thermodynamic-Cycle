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

def Rankine(fluid,T_cond,T_evap,T_coolant_in,T_heater_in,eff_pump,eff_turb,comp,output):
    
    # Fluid Specifications
    T_sh = 10                   # Superheating temperature (K)
    T_amb = 273+27
    m_dot_ORC = 5               # Mass flow rate of working fluid(s)
    m_dot_coolant = 120         # Mass flow rate of coolant
    m_dot_heater = 25           # Mass flow rate of evaporator fluid
    cp_coolant = 1000           # Cp of Air (J/kg K)
    cp_heater = 4186 
    T_3_sh = T_evap + T_sh
    
    # State calculations
    r_crit = RP.REFPROPdll(fluid,"","TC",MASS_BASE_SI, 0,0,0,0, comp)
    T_crit = r_crit.Output[0]
    
    # Inputs check
    if T_evap > T_crit:
        print("T_crit = " + str(T_crit))
        print("Error: Evaporator temperature is higher than critical temperature.")
        return
    
    r1 = RP.REFPROPdll(fluid,"TQ","P;H;S",MASS_BASE_SI, 0,0,T_cond,0, comp)
    P_1 = r1.Output[0]
    H_1 = r1.Output[1]
    S_1 = r1.Output[2]
    
    r3 = RP.REFPROPdll(fluid,"TQ","P;H;S",MASS_BASE_SI, 0,0,T_evap,1, comp)
    P_3 = r3.Output[0]
    H_3 = r3.Output[1]
    S_3 = r3.Output[2]
    
    r3_sh = RP.REFPROPdll(fluid,"TP","H;S",MASS_BASE_SI, 0,0,T_3_sh,P_3, comp)
    H_3_sh = r3_sh.Output[0]
    S_3_sh = r3_sh.Output[1]
    
    r22 = RP.REFPROPdll(fluid,"PQ","T;S",MASS_BASE_SI, 0,0,P_3,0, comp)
    T_22 = r22.Output[0]
    S_22 = r22.Output[1]
    
    r2s = RP.REFPROPdll(fluid,"SP","T;H",MASS_BASE_SI, 0,0,S_1,P_3, comp)
    T_2s = r2s.Output[0]
    H_2s = r2s.Output[1]
    
    r44 = RP.REFPROPdll(fluid,"PQ","T;S",MASS_BASE_SI, 0,0,P_1,1, comp)
    T_44 = r44.Output[0]
    S_44 = r44.Output[1]
    
    r4s = RP.REFPROPdll(fluid,"SP","T;H",MASS_BASE_SI, 0,0,S_3_sh,P_1, comp)
    T_4s = r4s.Output[0]
    H_4s = r4s.Output[1]
    
    H_2 = (H_2s - H_1)/eff_pump + H_1
    S_2 = S_1/eff_pump
    r2 = RP.REFPROPdll(fluid,"HS","T",MASS_BASE_SI, 0,0,H_2,S_2, comp)
    T_2 = r2.Output[0]
    
    H_4 = H_3_sh - eff_turb * (H_3_sh - H_4s)
    S_4 = S_3_sh/eff_turb
    r4 = RP.REFPROPdll(fluid,"HS","T",MASS_BASE_SI, 0,0,H_4,S_4, comp)
    T_4 = r4.Output[0]
    
    T_coolant_range = []
    S_cc = np.linspace(S_1, S_4, 11)
    H_prev_cc = H_1
    T_prev_cc = T_coolant_in
    for i in S_cc:
        r_cc = RP.REFPROPdll(fluid,"SP","H",MASS_BASE_SI, 0,0,i,P_1, comp)
        T_out_cc = T_prev_cc + (m_dot_ORC * (r_cc.Output[0] - H_prev_cc))/(m_dot_coolant*cp_coolant)
        T_coolant_range.append(T_out_cc)
        H_prev_cc = r_cc.Output[0]
        T_prev_cc = T_out_cc
    
    # Heater Curve
    T_heater_range = []
    S_hc = np.linspace(S_3_sh, S_2, 11)
    H_prev_hc = H_3_sh
    T_prev_hc = T_heater_in
    for i in S_hc:
        r_hc = RP.REFPROPdll(fluid,"SP","H",MASS_BASE_SI, 0,0,i,P_3, comp)
        T_out_hc = T_prev_hc - (m_dot_ORC * (H_prev_hc - r_hc.Output[0]))/(m_dot_heater*cp_heater)
        T_heater_range.append(T_out_hc)
        H_prev_hc = r_hc.Output[0]
        T_prev_hc = T_out_hc
    
    if output == "H":
        return H_1, H_2, H_3_sh, H_4, "Enthalpies (J / kg): "
    elif output == "P":
        return P_1, P_3, P_3, P_1, "Pressures (Pa): "
    elif output == "S":
        return S_1, S_2, S_3_sh, S_4, "Entropies (J / kg K): "
    elif output == "T":
        return T_cond, T_2, T_evap, T_4, "Temperatures (K): "
    elif output == "cycle eff":
        eff_TH = 1 - (H_4 - H_1)/(H_3_sh - H_2)
        return eff_TH, "Cycle efficiency = "
    elif output == "net work":
        return (H_3_sh - H_4 - (H_2 - H_1)), "Net work (J): "
    elif output == "overall eff":
        eff_OA = (H_3_sh - H_4 - (H_2 - H_1))/(m_dot_heater * cp_heater * (T_heater_in-T_amb))
        return eff_OA, "Overall efficiency: "
    elif output == "graph":
        S_states = [S_1, S_2, S_22, S_3, S_3_sh, S_4, S_44, S_1]
        T_states = [T_cond, T_2, T_22, T_evap, T_3_sh, T_4, T_44, T_cond]
        
        return S_states, T_states, fluid, comp, T_crit, S_cc, S_hc, T_coolant_range, T_heater_range

def graph(S_states, T_states, fluid, comp, T_crit, S_cc, S_hc, T_coolant_range, T_heater_range, options):    
    graph_index = [1,2,22,33,3,4,44]
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
    for i in range(7):
        plt.annotate(str(graph_index[i]),(S_states[i], T_states[i]))
    plt.xlabel('s (J / kg K)')
    plt.ylabel('T (K)')
    plt.plot(S_sat,Temp_graph)
    plt.plot(S_states, T_states, 'o-')
    plt.plot(S_cc,T_coolant_range)
    plt.plot(S_hc,T_heater_range)
    if options == "limit":
        plt.xlim([S_states[0]-100,S_states[5]+100])
        plt.ylim([T_states[3]-100,T_crit+20])

def reader(args, out, graph_options):
    for i in out:
        if len(i) == 1: 
            var1, var2, var3, var4, title = Rankine(*args, i)
            print(title + str(round(var1,2)) + ", " + str(round(var2,2)) + ", " 
                  + str(round(var3,2)) + ", " + str(round(var4,2)))
        elif i != "graph":
            var = Rankine(*args, i)
            print(str(var[1]) + str(var[0]))
        elif i == "graph":
            graph_inputs = Rankine(*args, i)
            graph(*graph_inputs, graph_options)
        else:
            return "Invalid Output"
            
def values(func, args, out):
    var = func(*args, out)
    return var

inputs_1 = ["R134A * R245FA",273+20,273+75,273+15,273+120,1,1,[0.5,0.5]]
inputs_2 = ["R245FA",273+20,273+75,273+15,273+120,1,1,[1]]
inputs_3 = ["R134A",273+20,273+75,273+15,273+120,1,1,[1]]

#reader(inputs_1,["cycle eff"],0)
#reader(inputs_1,["graph"],0)
reader(inputs_2,["graph"],0)
reader(inputs_3,["graph"],"limit")


# =============================================================================
# eff_results = []
# comp_range = np.linspace(0,1,31)
# 
# for i in comp_range:
#     inputs_1[-1] = [i,1-i]
#     var = values(Rankine, inputs_1, "overall eff")
#     eff_results.append(var[0])
# 
# inputs_1[-1] = [0.5,0.5]
# fig1, ax1 = plt.subplots()    
# ax1.plot(comp_range,eff_results,'o')
# plt.title('Overall Efficiency against composition of R245FA')
# plt.ylabel('Overall Efficiency')
# plt.xlabel('x_R245FA')
# =============================================================================

