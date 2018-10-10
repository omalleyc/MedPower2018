# -*- coding: utf-8 -*-
"""
Created on Sun May 13 13:43:27 2018

@author: omalleyc
"""

import LoadData
import DC_OPF_Function
import DC_OPF_RO_Function
import OGF_SS_Function
import Coupled_Function


# Choose A Case
[Elecdata,Gasdata]=LoadData.Load_Case1() 
[Elecdata,Gasdata]=LoadData.Load_Case2()
[Elecdata,Gasdata]=LoadData.Load_Case5()

#[Elecdata,Gasdata]=LoadData.Load_Case2_v() # Variation on Case 2
#[Elecdata,Gasdata]=LoadData.Load_Case5_v() # Variation on Case 5


# Add the shedding Weight
Elecdata.Gen=Elecdata.Gen.assign(ShedWeight=Elecdata.Gen.CostCoeff_1.rank(ascending=True))
Elecdata.Gen=Elecdata.Gen.assign(AddWeight=Elecdata.Gen.CostCoeff_1.rank(ascending=False))

# Solve Initial DC OPF
DC_OPF_Function.DC_OPF(Elecdata)

# Solve the robust Optimization
DC_OPF_RO_Function.DC_OPF_RO(Elecdata)

# Solve the OGF Problem without the Gas Generators
OGF_SS_Function.OGF_SS(Gasdata)

# Use the Gas Generator data and solve the OGF problem
GasGenData=Elecdata.Gen[Elecdata.Gen['FuelType']=='Gas']
GasGenData=OGF_SS_Function.OGF_SS(Gasdata,GasGenData)
Elecdata.Gen=Elecdata.Gen.combine_first(GasGenData)

# Solve OPG Again to find the redispatch non gas generation
DC_OPF_Function.DC_OPF(Elecdata)


# Solve Coupled Problem for Reference
Coupled_Function.CoupledOPFOGF(Elecdata,Gasdata)

# Print Some results to the screen
print('=====================================================================')
print('\n\t\t\tDispatch of Generators\n')
print('=====================================================================')
print(Elecdata.Gen[['Coupled_RES','OGF_P_output','DC_OPF_RES','SS_OGF_Res_shed']].round(decimals=1))
print('=====================================================================')
print('\n\t\tResults of Robust Optimization\n')
print('=====================================================================')
print(Elecdata.Gen[['DC_RO_OPF_P_LB','DC_RO_OPF_P_UB','DC_OPF_RES']].round(decimals=1))
print('=====================================================================')
print('\n\t\tOGF Load Shedding/Adding Actions\n')
print('=====================================================================')
print(Elecdata.Gen[['P_Add','P_Shed','DC_OPF_RES']].round(decimals=1))


Gen=Elecdata.Gen
Obj_OPF=Elecdata.Gen.DC_OPF_RES.multiply(Elecdata.Gen.CostCoeff_1).sum().round(decimals=1)
Obj_Coupled = Elecdata.Gen.Coupled_RES.multiply(Elecdata.Gen.CostCoeff_1).sum().round(decimals=1)
Obj_OPF_OGF= Elecdata.Gen.OGF_P_output.multiply(Elecdata.Gen.CostCoeff_1).sum().round(decimals=1)
print('=====================================================================')
print('\n\t\tObjective Functions\n')
print('=====================================================================')
print('Original OPF Objective ='+str(Obj_OPF))
print('OPF after OGF Objective ='+str(Obj_OPF_OGF))
print('Coupled Objective ='+str(Obj_Coupled))
