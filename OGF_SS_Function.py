# -*- coding: utf-8 -*-
"""
Created on Sun May 13 16:13:45 2018

@author: omalleyc
"""

import pandas as pd
import pyomo.environ as pm
import numpy as np
from pyomo.opt import SolverStatus, TerminationCondition

class expando(object):
    pass

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
    
    initial_from=Gasdata.Params.Initialize_from
    Initial=Initial_Values(Gasdata,initial_from)

    m.pi       = pm.Var(m.node_set ,bounds = lambda m,i : (Nodes.Pmin[i], Nodes.Pmax[i]), initialize= lambda m,i : Initial.Nodes.loc[i,'pi'])
    m.mij      = pm.Var(m.pipe_set ,bounds = (None,None),  initialize=lambda m,i :Initial.Pipes.loc[i,'mij'])
    m.c        = pm.Var(m.comp_set ,bounds =  lambda m,i : (Compressors.Cmin[i],Compressors.Cmax[i]),initialize =lambda m,i :Initial.Comps.loc[i,'c'])
    m.mc_in    = pm.Var(m.comp_set ,bounds = (None,None),  initialize=lambda m,i :Initial.Comps.loc[i,'mc_in']) # Gas flow into compressor, i.e. out of node
    m.mc_out   = pm.Var(m.comp_set ,bounds = (None,None),  initialize=lambda m,i :Initial.Comps.loc[i,'mc_out']) # Gas flow out of compressor i.e., into node
    m.mc       = pm.Var(m.comp_set ,bounds = (0,None),     initialize=lambda m,i :Initial.Comps.loc[i,'mc'])
    m.GasLoad  = pm.Var(m.node_set ,bounds = (None,None),  initialize=lambda m,i :Initial.Nodes.loc[i,'GasLoad'])
    
    if not(Gen.empty):
        
        m.Gen_set         = pm.Set(initialize = Gen.index.tolist())

        m.Gen_shed = pm.Var(m.Gen_set , bounds = lambda m,i : (0,Gen.DC_OPF_RES[i]-Gen.DC_RO_OPF_P_LB[i]),initialize=0)
        m.Gen_add  = pm.Var(m.Gen_set ,  bounds = lambda m,i : (0,Gen.DC_RO_OPF_P_UB[i]-Gen.DC_OPF_RES[i]),initialize=0)




    
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
            LHS=LHS + sum( Gen.Power_to_Gas_Norm[k]*(Gen.DC_OPF_RES[k]-m.Gen_shed[k]+m.Gen_add[k])
                for k in m.Gen_set if Gen.Gas_Node[k]==i)
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
    
    if not(Gen.empty):
        def Net_Shedding(model):
            LHS=0
            LHS = LHS + sum(m.Gen_shed[k]-m.Gen_add[k]   for k in m.Gen_set)
            return LHS>=0
        m.net_shedding=pm.Constraint(rule=Net_Shedding)


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
    opt.options['print_level']=5
    opt.options['tol']=1e-14
    opt.options['constr_viol_tol']=1e-10
    opt.options['nlp_scaling_method']='none'
    # Optimize
    results = opt.solve(m, tee=True)
 
    
    status=0
    if (results.solver.status == SolverStatus.ok) and (results.solver.termination_condition == TerminationCondition.optimal):
        print('Model Solved to Optimality')
        status=1
    # Do something when the solution in optimal and feasible
    elif(results.solver.termination_condition == TerminationCondition.infeasible):
        print('Model is infeasible')        
        return
    # Do something when model in infeasible
    else:
        print('Solver Status: ',  results.solver.status)
        return
        
        
    Compressor_c_Res       = dict([[i,np.round(m.c[i].value,decimals=8)]   for i in m.c])
    Compressor_mc_Res      = dict([[i,np.round(m.mc[i].value,decimals=8)]   for i in m.mc])
    Compressor_min_Res     = dict([[i,np.round(m.mc_in[i].value,decimals=8)]   for i in m.mc_in])
    Compressor_mout_Res    = dict([[i,np.round(m.mc_out[i].value,decimals=8)]   for i in m.mc_out])
    
    Nodes_pi_Res       = dict([[i,np.round(m.pi[i].value,decimals=8)]   for i in m.pi])
    Nodes_load_Res       = dict([[i,np.round(m.GasLoad[i].value,decimals=8)]   for i in m.GasLoad])
    
    Pipes_mij_Res       = dict([[i,np.round(m.mij[i].value,decimals=8)]   for i in m.mij])
    
    c_str    = 'SS_OGF_Res_c'
    mc_str   = 'SS_OGF_Res_mc'
    min_str  = 'SS_OGF_Res_min'
    mout_str = 'SS_OGF_Res_mout'
    
    mij_str  = 'SS_OGF_Res_mij'
    pi_str   = 'SS_OGF_Res_pi'
    load_str = 'SS_OGF_Res_load'
    
    
    if not(Gen.empty):
        FinalOutput=Gen.DC_OPF_RES
        
        Gen_shed = dict([[i,np.round(m.Gen_shed[i].value,decimals=8)]   for i in m.Gen_shed])
        Gen_add = dict([[i,np.round(m.Gen_add[i].value,decimals=8)]   for i in m.Gen_add])
        
        Gen=Gen.assign(SS_OGF_Res_shed=pd.Series(Gen_shed))
        Gen=Gen.assign(SS_OGF_Res_add=pd.Series(Gen_add))
        
        
        FinalOutput=FinalOutput.add(-pd.Series(Gen_shed),fill_value=0)
        FinalOutput=FinalOutput.add(pd.Series(Gen_add),fill_value=0)      
                 
        Gen=Gen.assign(SS_OGF_P_output=FinalOutput)
        
        c_str    = c_str    +'_wgen'
        mc_str   = mc_str   +'_wgen'
        min_str  = min_str  +'_wgen'
        mout_str = mout_str +'_wgen'
        
        mij_str  = mij_str  +'_wgen'
        pi_str   = pi_str   +'_wgen'
        load_str = load_str + '_wgen'
        
    Temp= {c_str    : pd.Series(Compressor_c_Res),
           mc_str   : pd.Series(Compressor_mc_Res),
           min_str  : pd.Series(Compressor_min_Res),
           mout_str : pd.Series(Compressor_mout_Res),}
    Gasdata.Compressors    =  Gasdata.Compressors.assign(**Temp)
    
    Gasdata.Pipes    =  Gasdata.Pipes.assign(**{mij_str : pd.Series(Pipes_mij_Res)})
    
    Temp={load_str : pd.Series(Nodes_load_Res),
          pi_str   : pd.Series(Nodes_pi_Res)}
    Gasdata.Nodes    =  Gasdata.Nodes.assign(**Temp)
    
    Gasdata.Params.status=status
    
    Output=expando()
    
    Output.obj=m.objective()
    
    
    if not(Gen.empty):
        Output.GasGenData=Gen
    
    return Output

def Initial_Values(Gasdata,initial_from):
    
    
    Nodes=pd.DataFrame(index=Gasdata.Nodes.index.tolist())
    Pipes=pd.DataFrame(index=Gasdata.Pipes.index.tolist())
    Comps=pd.DataFrame(index=Gasdata.Compressors.index.tolist())
    # Default is nonoe
    Nodes['pi']=None
    Nodes['GasLoad']=None
    Comps['mc_in']=None
    Comps['mc_out']=None
    Comps['mc']=None
    Comps['c']=None
    Pipes['mij']=None
    
    if initial_from=='SS_OGF_wgen':
        print('Initializing from SS_OGF_wgen')
        Nodes['pi']=Gasdata.Nodes['SS_OGF_Res_pi_wgen']
        Nodes['GasLoad']=Gasdata.Nodes['SS_OGF_Res_load_wgen']
        Pipes['mij']=Gasdata.Pipes['SS_OGF_Res_mij_wgen']
        Comps['c']=Gasdata.Compressors['SS_OGF_Res_c_wgen']
        Comps['mc']=Gasdata.Compressors['SS_OGF_Res_mc_wgen']
        Comps['mc_in']=Gasdata.Compressors['SS_OGF_Res_min_wgen']
        Comps['mc_out']=Gasdata.Compressors['SS_OGF_Res_mout_wgen']
    elif initial_from=='Flat_Start':
        print('Initializing from Flat Start')
        Nodes['pi']=1
        Nodes['GasLoad']=0
        Pipes['mij']=0
        Comps['c']=1
        Comps['mc']=0
        Comps['mc_in']=0
        Comps['mc_out']=0
    elif initial_from=='Coupled':
        print('Initializing from Coupled')
        Nodes['pi']=Gasdata.Nodes['Coupled_Res_pi']
        Nodes['GasLoad']=Gasdata.Nodes['Coupled_Res_load']
        Pipes['mij']=Gasdata.Pipes['Coupled_Res_mij']
        Comps['c']=Gasdata.Compressors['Coupled_Res_c']
        Comps['mc']=Gasdata.Compressors['Coupled_Res_mc']
        Comps['mc_in']=Gasdata.Compressors['Coupled_Res_min']
        Comps['mc_out']=Gasdata.Compressors['Coupled_Res_mout']
    elif initial_from=='SS_OGF':
        print('Initializing from SS_OGF')
        # ize from SS_OGF
        Nodes['pi']=Gasdata.Nodes['SS_OGF_Res_pi']
        Nodes['GasLoad']=Gasdata.Nodes['SS_OGF_Res_load']
        Pipes['mij']=Gasdata.Pipes['SS_OGF_Res_mij']
        Comps['c']=Gasdata.Compressors['SS_OGF_Res_c']
        Comps['mc']=Gasdata.Compressors['SS_OGF_Res_mc']
        Comps['mc_in']=Gasdata.Compressors['SS_OGF_Res_min']
        Comps['mc_out']=Gasdata.Compressors['SS_OGF_Res_mout']
    else:
        print('Initialized with none')
        
        
    Initial=expando()
    Initial.Nodes=Nodes
    Initial.Pipes=Pipes
    Initial.Comps=Comps
    return Initial












