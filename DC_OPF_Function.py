# -*- coding: utf-8 -*-
"""
Created on Sun May 13 13:46:47 2018

@author: omalleyc
"""

# -*- coding: utf-8 -*-
"""
Created on Fri May 11 15:51:02 2018

@author: omalleyc
"""


import pandas as pd
import pyomo.environ as pm
from pyomo.opt import SolverStatus, TerminationCondition
import numpy as np


def DC_OPF(Elecdata):
    Gen    = Elecdata.Gen
    Branch = Elecdata.Branch
    Bus    = Elecdata.Bus
    
    # If the OGF has been run and there are fixed values for the GFPPs change 
    # the limits
    if 'SS_OGF_P_output' in Elecdata.Gen.columns:
        print('Fixing output of Gas Generators')
        fn_max = lambda row:  np.where(row.FuelType!='Gas',row.Pmax_MW,row.SS_OGF_P_output).item(0)
        fn_min = lambda row:  np.where(row.FuelType!='Gas',row.Pmin_MW,row.SS_OGF_P_output).item(0)
        Gen=Gen.assign(Pmax_MW=Gen.apply(fn_max,axis=1).values)
        Gen=Gen.assign(Pmin_MW=Gen.apply(fn_min,axis=1).values)
    
    
    
    
    m = pm.ConcreteModel()
    
    m.gen_set    = pm.Set(initialize=Gen.index.tolist())
    m.branch_set = pm.Set(initialize=Branch.index.tolist())
    m.bus_set    = pm.Set(initialize=Bus.index.tolist())
    
    m.P   = pm.Var(m.gen_set,bounds = lambda m,i : (Gen.Pmin_MW[i],Gen.Pmax_MW[i]), initialize=lambda m,i : (Gen.Pmin_MW[i]))
    m.Pij = pm.Var(m.branch_set,bounds = lambda m,i : (-Branch.RateA_MVA[i], Branch.RateA_MVA[i]), initialize=0)
    m.th  = pm.Var(m.bus_set,bounds = (-np.pi,np.pi),initialize=0)
    
    
    def PowerBal_constr(model,i):
        PowerBal = -Bus.PD_MW[i] \
        + sum(m.P[k]   for k in Gen[Gen.Gen_Bus==i].index.tolist()) \
        - sum(m.Pij[k] for k in Branch[Branch.From_Bus==i].index.tolist()) \
        + sum(m.Pij[k] for k in Branch[Branch.To_Bus==i].index.tolist()) \
        ==0
       
        return PowerBal          
    m.PowerBal_constr=pm.Constraint(m.bus_set,rule=PowerBal_constr)
    
    def Branch_Flow(model,i):
        From_ix=Branch.From_Bus[i]
        To_ix=Branch.To_Bus[i]
        B=Elecdata.Params.BaseMVA/Branch.BR_X_PU[i]
        return m.Pij[i]==(B)*(m.th[From_ix]-m.th[To_ix])
    
    m.branchflow=pm.Constraint(m.branch_set,rule=Branch_Flow)
    
    def Slackbus(model,i):
        if Bus.Bus_Type[i]==3:
            return m.th[i]==0
        else:
            return pm.Constraint.Skip
    
    m.slack=pm.Constraint(m.bus_set,rule=Slackbus)
    
    
    def DC_OPF_Obj(model):
        return sum(m.P[i]*Gen.CostCoeff_1[i]for i in m.gen_set) # Constant term doesnt matter for optimization
    
    m.objective = pm.Objective(rule=DC_OPF_Obj, sense=pm.minimize, doc='Define objective function') 
    
    opt = pm.SolverFactory('gurobi')
    opt.options['print_level']=0     
    
    results = opt.solve(m, tee=True)
    
    status=0
    if (results.solver.status == SolverStatus.ok) and (results.solver.termination_condition == TerminationCondition.optimal):
        print('Model Solved to Optimality')
        status=1
    # Do something when the solution in optimal and feasible
    elif(results.solver.termination_condition == TerminationCondition.infeasible):
        print('Model is infeasible')        
    # Do something when model in infeasible
    else:
        print('Solver Status: ',  results.solver.status)

    
    
    
    Gen_Res     = dict([[i,np.round(m.P[i].value,decimals=3)]   for i in m.P])
    Branch_Res  = dict([[i,np.round(m.Pij[i].value,decimals=3)]   for i in m.Pij])
    Bus_Res     = dict([[i,np.round(m.th[i].value,decimals=3)]   for i in m.th])
    
    if 'SS_OGF_P_output' in Elecdata.Gen.columns:
        Elecdata.Gen = Elecdata.Gen.assign(SS_OGF_P_output=pd.Series(Gen_Res))
        
    else:
        Elecdata.Gen    =  Elecdata.Gen.assign    (DC_OPF_RES=pd.Series(Gen_Res))
        Elecdata.Bus    =  Elecdata.Bus.assign    (DC_OPF_RES=pd.Series(Bus_Res))
        Elecdata.Branch =  Elecdata.Branch.assign (DC_OPF_RES=pd.Series(Branch_Res))
    

    Elecdata.Params.DC_OPF_Obj = m.objective()
    
    Elecdata.Params.status=status
    
    return m.objective()

# Code for getting Duals

## Gurobi 
#m.rc      = pmc.Suffix(direction=pmc.Suffix.IMPORT)
#m.dual    = pmc.Suffix(direction=pmc.Suffix.IMPORT_EXPORT)
#opt.options['extract_reduced_costs']=True
#for i in data:
#    print('Fixed Dual of ',i+1,'is ',m.dual.get(m.mynewconstraint[i]))
#print('Objective is',m.objective.expr())
##print(m.dual.get(m.mynewconstraint1))
##print(m.dual.get(m.mynewconstraint2))
##print(m.dual.get(m.mynewconstraint3))
#m.rc.pprint()
#m.dual.pprint()
