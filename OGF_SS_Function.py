# -*- coding: utf-8 -*-
"""
Created on Sun May 13 16:13:45 2018

@author: omalleyc
"""

import pandas as pd
import pyomo.environ as pm
import numpy as np
from pyomo.opt import SolverStatus, TerminationCondition


def OGF_SS(Gasdata,Gen=pd.DataFrame()):
            
    Pipes=Gasdata.Pipes
    Compressors=Gasdata.Compressors
    Nodes=Gasdata.Nodes
    Params=Gasdata.Params
    
    m = pm.ConcreteModel()
        
    # Sets
    m.node_set = pm.Set(initialize=Nodes.index.tolist())
    m.comp_set = pm.Set(initialize=Compressors.index.tolist())
    m.pipe_set = pm.Set(initialize=Pipes.index.tolist())
    
    
    m.pi       = pm.Var(m.node_set ,bounds = lambda m,i : (Nodes.Pmin[i], Nodes.Pmax[i]), initialize=1)
    m.mij      = pm.Var(m.pipe_set ,bounds = (None,None),initialize=0)
    m.c        = pm.Var(m.comp_set ,bounds =  lambda m,i : (Compressors.Cmin[i],Compressors.Cmax[i]),initialize =1)
    m.mc_in    = pm.Var(m.comp_set ,bounds = (None,None),initialize=0) # Gas flow into compressor, i.e. out of node
    m.mc_out   = pm.Var(m.comp_set ,bounds = (None,None),initialize=0) # Gas flow out of compressor i.e., into node
    m.mc       = pm.Var(m.comp_set ,bounds = (None,None),initialize=0)
    m.GasLoad  = pm.Var(m.node_set ,bounds = (0,None)   ,initialize=0)
    
    if not(Gen.empty):
        
        m.Gen_set         = pm.Set(initialize = Gen.index.tolist())

        m.Gen_shed = pm.Var(m.Gen_set , bounds = lambda m,i : (0,Gen.P_Shed[i]))
        m.Gen_add  = pm.Var(m.Gen_set ,  bounds = lambda m,i : (0,Gen.P_Add[i]))
        
    
    def Weymouth(model, k):
        From=Pipes.loc[k].From
        To=Pipes.loc[k].To
        return m.pi[From]**2 - m.pi[To]**2 == Pipes.SS_Constant[k]*np.abs(m.mij[k])*m.mij[k]
    m.weymouth = pm.Constraint(m.pipe_set, rule=Weymouth)
    

    def Compressor_Ratio(model, k):
        From = Compressors.loc[k].From
        To   = Compressors.loc[k].To
        return m.pi[From]*m.c[k]== m.pi[To]
    m.comp_ratio = pm.Constraint(m.comp_set, rule=Compressor_Ratio)
    

    def Compressor_Balance(model, k):
        return m.mc_in[k]==m.mc[k]+m.mc_out[k]
    m.comp_bal = pm.Constraint(m.comp_set, rule=Compressor_Balance)
    

    def Compressor_Consumption(model, k):
        return m.mc[k]==m.mc_out[k]*Compressors.K[k]*(m.c[k]**Compressors.exponent[k]-1)
    m.comp_cons = pm.Constraint(m.comp_set, rule=Compressor_Consumption)

    
    def NodalLoad(model,i):
        LHS=Nodes.Load[i]-m.GasLoad[i]
        if not(Gen.empty):
            LHS=LHS + sum( Gen.Power_to_Gas_Norm[k]*Gen.DC_OPF_RES[k] for k in m.Gen_set if Gen.Gas_Node[k]==i) +\
                      sum(-Gen.Power_to_Gas_Norm[k]*m.Gen_shed[k]   for k in m.Gen_set if Gen.Gas_Node[k]==i) +\
                      sum( Gen.Power_to_Gas_Norm[k]*m.Gen_add[k]    for k in m.Gen_set if Gen.Gas_Node[k]==i)
        return LHS==0     
    
    m.nodal_load=pm.Constraint(m.node_set,rule=NodalLoad)

    def NodalBalance(model,i):
        if i==Params.Slacknode: # Slack Node
            return m.pi[i]==Params.Slackpressure
        
        LHS = m.GasLoad[i]\
              + sum(m.mij[k]     for k in Pipes[Pipes.From==i].index.tolist()) \
              - sum(m.mij[k]     for k in Pipes[Pipes.To==i].index.tolist())
              
        if not(Compressors.empty):
            LHS = LHS \
                  + sum(m.mc_in[k]   for k in Compressors[Compressors.From==i].index.tolist()) \
                  - sum(m.mc_out[k]  for k in Compressors[Compressors.To==i].index.tolist()) 
                  
        return LHS==0
    m.nodalbalance=pm.Constraint(m.node_set,rule=NodalBalance)
    


    def objective_rule(m):
        Temp = 0.0
        #Temp= sum(m.mc[k] for k in m.comp_set) 
        
        if not(Gen.empty):
            Temp=Temp+ \
            sum(1/Gen.CostCoeff_1[k]*m.Gen_shed[k] for k in m.Gen_set)+ \
            -sum(1/Gen.CostCoeff_1[k]*m.Gen_add[k] for k in m.Gen_set) 
            
        return Temp
    
    m.objective = pm.Objective(rule=objective_rule, sense=pm.minimize, doc='Define objective function') 
    
    opt = pm.SolverFactory('ipopt')
    opt.options['print_level']=0 
    # Optimize
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
        
        
    Compressor_c_Res       = dict([[i,np.round(m.c[i].value,decimals=8)]   for i in m.c])
    Compressor_mc_Res      = dict([[i,np.round(m.mc[i].value,decimals=8)]   for i in m.mc])
    Compressor_min_Res     = dict([[i,np.round(m.mc_in[i].value,decimals=8)]   for i in m.mc_in])
    Compressor_mout_Res    = dict([[i,np.round(m.mc_out[i].value,decimals=8)]   for i in m.mc_out])
    
    Nodes_pi_Res       = dict([[i,np.round(m.pi[i].value,decimals=8)]   for i in m.pi])
    Pipes_mij_Res       = dict([[i,np.round(m.mij[i].value,decimals=8)]   for i in m.mij])
    
    c_str    = 'SS_OGF_Res_c'
    mc_str   = 'SS_OGF_Res_mc'
    min_str  = 'SS_OGF_Res_min'
    mout_str = 'SS_OGF_Res_mout'
    
    mij_str  = 'SS_OGF_Res_mij'
    pi_str   = 'SS_OGF_Res_pi'
    
    
    if not(Gen.empty):
        FinalOutput=Gen.DC_OPF_RES
        
        Gen_shed = dict([[i,np.round(m.Gen_shed[i].value,decimals=8)]   for i in m.Gen_shed])
        Gen_add = dict([[i,np.round(m.Gen_add[i].value,decimals=8)]   for i in m.Gen_add])
        
        Gen=Gen.assign(SS_OGF_Res_shed=pd.Series(Gen_shed))
        Gen=Gen.assign(SS_OGF_Res_add=pd.Series(Gen_add))
        
        FinalOutput=FinalOutput.add(-pd.Series(Gen_shed),fill_value=0)
        FinalOutput=FinalOutput.add(pd.Series(Gen_add),fill_value=0)      
                 
        Gen=Gen.assign(OGF_P_output=FinalOutput)
        
        c_str    = c_str    +'_wgen'
        mc_str   = mc_str   +'_wgen'
        min_str  = min_str  +'_wgen'
        mout_str = mout_str +'_wgen'
        
        mij_str  = mij_str  +'_wgen'
        pi_str   = pi_str   +'_wgen'
        
    Temp= {c_str    : pd.Series(Compressor_c_Res),
           mc_str   : pd.Series(Compressor_mc_Res),
           min_str  : pd.Series(Compressor_min_Res),
           mout_str : pd.Series(Compressor_mout_Res)}
    Gasdata.Compressors    =  Gasdata.Compressors.assign(**Temp)
    
    Gasdata.Pipes    =  Gasdata.Pipes.assign(**{mij_str : pd.Series(Pipes_mij_Res)})
    Gasdata.Nodes    =  Gasdata.Nodes.assign(**{pi_str  : pd.Series(Nodes_pi_Res)})
    
    Gasdata.Params.status=status
    
    if not(Gen.empty):
        return Gen

    
    
    

    













