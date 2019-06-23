# -*- coding: utf-8 -*-
"""
Created on Sun Jul 22 11:54:39 2018

@author: omalleyc
"""


import pandas as pd
import numpy as np
import pyomo.environ as pm
from pyomo.opt import SolverStatus, TerminationCondition


class expando(object):
    pass

def CoupledOPFOGF(Elecdata,Gasdata,initial_from):
            
    Gen    = Elecdata.Gen
    Branch = Elecdata.Branch
    Bus    = Elecdata.Bus
    
    
    Pipes=Gasdata.Pipes
    Compressors=Gasdata.Compressors
    Nodes=Gasdata.Nodes
    Params=Gasdata.Params
    
    
    Initial=Initial_Values(Elecdata,Gasdata,initial_from)
    
    m = pm.ConcreteModel()
    
        
    m.gen_set    = pm.Set(initialize=Gen.index.tolist())
    m.branch_set = pm.Set(initialize=Branch.index.tolist())
    m.bus_set    = pm.Set(initialize=Bus.index.tolist())
    
    m.P   = pm.Var(m.gen_set,bounds = lambda m,i : (Gen.Pmin_MW[i],Gen.Pmax_MW[i]), initialize=lambda m,i : Initial.Gen.loc[i,'P_init'])
    m.Pij = pm.Var(m.branch_set,bounds = lambda m,i : (-Branch.RateA_MVA[i], Branch.RateA_MVA[i]), initialize=lambda m,i : Initial.Branch.loc[i,'Pij_init'])
    m.th  = pm.Var(m.bus_set,bounds = (-np.pi,np.pi),initialize=lambda m,i : Initial.Bus.loc[i,'th_init'])
        
    # Sets
    m.node_set = pm.Set(initialize=Nodes.index.tolist())
    m.comp_set = pm.Set(initialize=Compressors.index.tolist())
    m.pipe_set = pm.Set(initialize=Pipes.index.tolist())
    
    m.pi       = pm.Var(m.node_set ,bounds = lambda m,i : (Nodes.Pmin[i], Nodes.Pmax[i]), initialize= lambda m,i : Initial.Nodes.loc[i,'pi_init'])
    m.mij      = pm.Var(m.pipe_set ,bounds = (None,None),initialize= lambda m,i : Initial.Pipes.loc[i,'mij_init'])
    m.c        = pm.Var(m.comp_set ,bounds =  lambda m,i : (Compressors.Cmin[i],Compressors.Cmax[i]),initialize= lambda m,i : Initial.Comps.loc[i,'c_init'])
    m.mc_in    = pm.Var(m.comp_set ,bounds = (None,None),initialize= lambda m,i : Initial.Comps.loc[i,'mc_in_init']) # Gas flow into compressor, i.e. out of node
    m.mc_out   = pm.Var(m.comp_set ,bounds = (None,None),initialize= lambda m,i : Initial.Comps.loc[i,'mc_out_init']) # Gas flow out of compressor i.e., into node
    m.mc       = pm.Var(m.comp_set ,bounds = (0,None),initialize= lambda m,i : Initial.Comps.loc[i,'mc_init'])
    m.GasLoad  = pm.Var(m.node_set ,bounds = (None,None)   ,initialize= lambda m,i : Initial.Nodes.loc[i,'GasLoad_init'])
    
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
        LHS=Nodes.Load[i]-m.GasLoad[i] + \
        sum(Gen.Power_to_Gas_Norm[k]*m.P[k] for k in m.gen_set  if Gen.Gas_Node[k]==i)
        
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
        Temp= sum(m.P[i]*Gen.CostCoeff_1[i]for i in m.gen_set)
        #sum(m.mc[k] for k in m.comp_set) +\

        return Temp
    
    m.objective = pm.Objective(rule=objective_rule, sense=pm.minimize, doc='Define objective function') 
    
    opt = pm.SolverFactory('ipopt')
    
    # Optimize
    results = opt.solve(m, tee=True)
    opt.options['tol']=1e-14
    opt.options['nlp_scaling_method']='none'
    
    if (results.solver.status == SolverStatus.ok) and (results.solver.termination_condition == TerminationCondition.optimal):
        print('Model Solved to Optimality')
    # Do something when the solution in optimal and feasible
    elif(results.solver.termination_condition == TerminationCondition.infeasible):
        print('Model is infeasible')        
        return
    # Do something when model in infeasible
    else:
        print('Solver Status: ',  results.solver.status)
        return
    
    
    Gen_Res     = dict([[i,np.round(m.P[i].value,decimals=8)]   for i in m.P])
    Branch_Res  = dict([[i,np.round(m.Pij[i].value,decimals=8)]   for i in m.Pij])
    Bus_Res     = dict([[i,np.round(m.th[i].value,decimals=8)]   for i in m.th])
    
    Elecdata.Gen    =  Elecdata.Gen.assign    (Coupled_RES=pd.Series(Gen_Res))
    Elecdata.Bus    =  Elecdata.Bus.assign    (Coupled_RES=pd.Series(Bus_Res))
    Elecdata.Branch =  Elecdata.Branch.assign (Coupled_RES=pd.Series(Branch_Res))
    
    Compressor_c_Res       = dict([[i,np.round(m.c[i].value,decimals=8)]   for i in m.c])
    Compressor_mc_Res      = dict([[i,np.round(m.mc[i].value,decimals=8)]   for i in m.mc])
    Compressor_min_Res     = dict([[i,np.round(m.mc_in[i].value,decimals=8)]   for i in m.mc_in])
    Compressor_mout_Res    = dict([[i,np.round(m.mc_out[i].value,decimals=8)]   for i in m.mc_out])
    
    Nodes_load_Res       = dict([[i,np.round(m.GasLoad[i].value,decimals=8)]   for i in m.GasLoad])
    
    Nodes_pi_Res       = dict([[i,np.round(m.pi[i].value,decimals=3)]   for i in m.pi])
    Pipes_mij_Res       = dict([[i,np.round(m.mij[i].value,decimals=3)]   for i in m.mij])
    
    c_str    = 'Coupled_Res_c'
    mc_str   = 'Coupled_Res_mc'
    min_str  = 'Coupled_Res_min'
    mout_str = 'Coupled_Res_mout'
    
    mij_str  = 'Coupled_Res_mij'
    pi_str   = 'Coupled_Res_pi'
    load_str = 'Coupled_Res_load'
    
        
    Temp= {c_str    : pd.Series(Compressor_c_Res),
           mc_str   : pd.Series(Compressor_mc_Res),
           min_str  : pd.Series(Compressor_min_Res),
           mout_str : pd.Series(Compressor_mout_Res)}
    Gasdata.Compressors    =  Gasdata.Compressors.assign(**Temp)
    
    Gasdata.Pipes    =  Gasdata.Pipes.assign(**{mij_str : pd.Series(Pipes_mij_Res)})
    
    Temp={pi_str  : pd.Series(Nodes_pi_Res),
          load_str: pd.Series(Nodes_load_Res)}
    Gasdata.Nodes    =  Gasdata.Nodes.assign(**Temp)
    
    return m.objective()


def Initial_Values(Elecdata,Gasdata,initial_from):
    
    Bus=pd.DataFrame(index=Elecdata.Bus.index.tolist())
    Branch=pd.DataFrame(index=Elecdata.Branch.index.tolist())
    Gen=pd.DataFrame(index=Elecdata.Gen.index.tolist())
    
    Nodes=pd.DataFrame(index=Gasdata.Nodes.index.tolist())
    Pipes=pd.DataFrame(index=Gasdata.Pipes.index.tolist())
    Comps=pd.DataFrame(index=Gasdata.Compressors.index.tolist())
    # Default is nonoe
    Bus['th_init']=None
    Branch['Pij_init']=None
    Gen['P_init']=None
    Nodes['pi_init']=None
    Nodes['GasLoad_init']=None
    Comps['mc_in_init']=None
    Comps['mc_out_init']=None
    Comps['mc_init']=None
    Comps['c_init']=None
    Pipes['mij_init']=None
    
    if initial_from=='Flat_Start':
        print('Initializing from Flat Start')
        Bus['th_init']=0
        Branch['Pij_init']=0
        Gen['P_init']=Elecdata.Gen['Pmin_MW']
        Nodes['pi_init']=1
        Nodes['GasLoad_init']=0
        Pipes['mij_init']=0
        Comps['c_init']=1
        Comps['mc_init']=0
        Comps['mc_in_init']=0
        Comps['mc_out_init']=0
    elif initial_from=='G_Optimal':
        print('Initializing from SS_OGF')
        # ize from SS_OGF
        Bus['th_init']=Elecdata.Bus['SS_OGF_th_output']
        Branch['Pij_init']=Elecdata.Branch['SS_OGF_Pij_output']
        Gen['P_init']=Elecdata.Gen['SS_OGF_P_output']
        Nodes['pi_init']=Gasdata.Nodes['SS_OGF_Res_pi_wgen']
        Nodes['GasLoad_init']=Gasdata.Nodes['SS_OGF_Res_load_wgen']
        Pipes['mij_init']=Gasdata.Pipes['SS_OGF_Res_mij_wgen']
        Comps['c_init']=Gasdata.Compressors['SS_OGF_Res_c_wgen']
        Comps['mc_init']=Gasdata.Compressors['SS_OGF_Res_mc_wgen']
        Comps['mc_in_init']=Gasdata.Compressors['SS_OGF_Res_min_wgen']
        Comps['mc_out_init']=Gasdata.Compressors['SS_OGF_Res_mout_wgen']
    else:
        print('Initialized with none')
        
        
    Initial=expando()
    Initial.Nodes=Nodes
    Initial.Pipes=Pipes
    Initial.Comps=Comps
    Initial.Bus=Bus
    Initial.Gen=Gen
    Initial.Branch=Branch
    return Initial
