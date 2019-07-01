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

import numpy as np
import pandas as pd


# Choose A Case
[Elecdata,Gasdata]=LoadData.Load_Case1() 
[Elecdata,Gasdata]=LoadData.Load_Case2()
[Elecdata,Gasdata]=LoadData.Load_Case5()

[Elecdata,Gasdata]=LoadData.Load_CaseGB()

#[Elecdata,Gasdata]=LoadData.Load_Case2_v() # Variation on Case 2
#[Elecdata,Gasdata]=LoadData.Load_Case5_v() # Variation on Case 5


# Solve Coupled Problem for Reference
initial_from='Flat_Start' # Flat_Start or G_Optimal
Coupled_Function.CoupledOPFOGF(Elecdata,Gasdata,initial_from)

# Add the shedding Weight
Elecdata.Gen=Elecdata.Gen.assign(ShedWeight=Elecdata.Gen.CostCoeff_1.rank(ascending=True))
Elecdata.Gen=Elecdata.Gen.assign(AddWeight=Elecdata.Gen.CostCoeff_1.rank(ascending=False))

# Solve Initial DC OPF
DC_OPF_Function.DC_OPF(Elecdata)

# Solve the robust Optimization
print('\n\nFinding the feasibility region\n')

# can read the flex set from a file if it was previously created
DC_OPF_RO_Function.DC_OPF_RO(Elecdata)
Elecdata.Gen.to_csv('Generators_w_Flex.csv')
#Elecdata.Gen=pd.read_csv('Generators_w_Flex.csv',index_col=0)

#DC_OPF_RO_Function.DC_OPF_RO_Col_Constraint_Gen(Elecdata)


# Solve the OGF Problem without the Gas Generators
print('\n\nSolving Gas network without Gas Generators\n')
OGF_SS_Function.OGF_SS(Gasdata)

# Use the Gas Generator data and solve the OGF problem
print('\n\nSolving Gas network with Gas Generators\n')
GasGenData=Elecdata.Gen[Elecdata.Gen['FuelType']=='Gas']

# Change the initial guess of the NLP if wanted
#Gasdata.Params.Initialize_from='Coupled'

Gas_output=OGF_SS_Function.OGF_SS(Gasdata,GasGenData)
GasGenData=Gas_output.GasGenData
Elecdata.Gen=Elecdata.Gen.combine_first(GasGenData)

# Solve OPG Again to find the redispatch non gas generation
DC_OPF_Function.DC_OPF(Elecdata)


# Print Some results to the screen
Gas_index=Elecdata.Gen.FuelType=='Gas'
print('=====================================================================')
print('\n\t\t\tDispatch of Generators\n')
print('=====================================================================')
print(Elecdata.Gen[['Coupled_RES','SS_OGF_P_output','DC_OPF_RES']].round(decimals=1))
print('=====================================================================')
print('\n\t\t\tDispatch of Gas Generators\n')
print('=====================================================================')
print(Elecdata.Gen.loc[Gas_index,['Coupled_RES','SS_OGF_P_output','DC_OPF_RES']].round(decimals=1))
print('=====================================================================')
print('\n\t\t\t Change in outpout of Gas Generators\n')
print('=====================================================================')
print(Elecdata.Gen.loc[Gas_index,['SS_OGF_Res_shed','SS_OGF_Res_add']].round(decimals=1))
print('=====================================================================')
print('\n\t\tResults of Robust Optimization\n')
print('=====================================================================')
print(Elecdata.Gen.loc[Gas_index,['DC_RO_OPF_P_LB','DC_OPF_RES','DC_RO_OPF_P_UB']].round(decimals=1))


Obj_OPF     = Elecdata.Gen.DC_OPF_RES.multiply(Elecdata.Gen.CostCoeff_1).sum().round(decimals=1)
Obj_Coupled = Elecdata.Gen.Coupled_RES.multiply(Elecdata.Gen.CostCoeff_1).sum().round(decimals=1)
Obj_OPF_OGF = Elecdata.Gen.SS_OGF_P_output.multiply(Elecdata.Gen.CostCoeff_1).sum().round(decimals=1)
print('=====================================================================')
print('\n\t\tObjective Functions\n')
print('=====================================================================')
print('Original OPF Objective ='+str(Obj_OPF))
print('OPF after OGF Objective ='+str(Obj_OPF_OGF))
print('Coupled Objective ='+str(Obj_Coupled))


# Check what solution is if start from g_optimal result
#initial_from='G_Optimal' # Flat_Start or G_Optimal
#Coupled_Function.CoupledOPFOGF(Elecdata,Gasdata,initial_from)

