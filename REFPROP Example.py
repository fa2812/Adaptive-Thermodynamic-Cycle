# -*- coding: utf-8 -*-
"""
Created on Wed Nov 11 11:26:21 2020
 
@author: esr1u09
"""
# ------------------------------------------------------------
# Initialisation
# ------------------------------------------------------------
# pip install ctREFPROP

import os,numpy as np
from ctREFPROP.ctREFPROP import REFPROPFunctionLibrary

os.environ['RPPREFIX'] = r'H:\IP 202122\REFPROP 10'

#Fluid properties:
RP = REFPROPFunctionLibrary(os.environ['RPPREFIX']) #Initialise REFPROP
RP.SETPATHdll(os.environ['RPPREFIX'])
MASS_BASE_SI = RP.GETENUMdll(0,"MASS BASE SI").iEnum
fluid = 'AIR'
 
r = RP.REFPROPdll(fluid,"TP","D",MASS_BASE_SI, 0,0,300.,101325., [1.0])
assert(r.ierr == 0)
Density=r.Output[0]
print(Density)
 