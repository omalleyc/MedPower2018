# -*- coding: utf-8 -*-
"""
Created on Fri May 11 11:40:04 2018

@author: omalleyc
"""

import pandas as pd
import pyomo.environ as pm
import numpy as np
import itertools
from pyomo.opt import SolverStatus, TerminationCondition

def DC_OPF_RO(Elecdata):
    
    Gen=Elecdata.Gen
    Branch=Elecdata.Branch
    Bus=Elecdata.Bus

    GasGenerators=Gen[Gen.FuelType=='Gas'].index.tolist()
    
    TruthTable=[i for i in itertools.product([0,1], repeat=len(GasGenerators))]
    TruthTable=pd.DataFrame(TruthTable,columns=GasGenerators)
    print('Created Truth Table')
    
    nscen=len(TruthTable)
    
    m = pm.ConcreteModel()
    
    m.gen_set    = pm.Set(initialize=Gen.index.tolist())
    m.branch_set = pm.Set(initialize=Branch.index.tolist())
    m.bus_set    = pm.Set(initialize=Bus.index.tolist())
    
    
    m.scen_set   = pm.RangeSet(0,nscen-1)
    m.gasgen_set = pm.Set(initialize=GasGenerators,within=m.gen_set,ordered=True)
    
    m.P_UB = pm.Var(m.gasgen_set,bounds= lambda m,i : (Gen.DC_OPF_RES[i],Gen.Pmax_MW[i]),initialize=lambda m,i : Gen.DC_OPF_RES[i] )
    m.P_LB = pm.Var(m.gasgen_set,bounds= lambda m,i : (Gen.Pmin_MW[i],Gen.DC_OPF_RES[i]),initialize=lambda m,i : Gen.DC_OPF_RES[i] )
    m.P    = pm.Var(m.gen_set    ,m.scen_set,bounds = lambda m,i,s : (Gen.Pmin_MW[i], Gen.Pmax_MW[i]), initialize=1)
    m.Pij  = pm.Var(m.branch_set ,m.scen_set,bounds = lambda m,i,s : (-Branch.RateA_MVA[i], Branch.RateA_MVA[i]), initialize=1)
    m.th   = pm.Var(m.bus_set   ,m.scen_set,bounds = (-np.pi,np.pi),initialize=0)
    m.P_Shed = pm.Var(m.gasgen_set,bounds=(0,None))
    m.P_Add  = pm.Var(m.gasgen_set,bounds=(0,None))
    

    def PowerBal_constr(model,i,s):
        PowerBal = -Bus.PD_MW[i] \
        + sum(m.P[k,s]   for k in Gen[Gen.Gen_Bus==i].index.tolist()) \
        - sum(m.Pij[k,s] for k in Branch[Branch.From_Bus==i].index.tolist()) \
        + sum(m.Pij[k,s] for k in Branch[Branch.To_Bus==i].index.tolist()) \
        ==0
        return PowerBal          
    m.PowerBal_constr=pm.Constraint(m.bus_set,m.scen_set,rule=PowerBal_constr)
    
    
    
    def Branch_Flow(model,i,s):
        From_ix=Branch.From_Bus[i]
        To_ix=Branch.To_Bus[i]
        B=Elecdata.Params.BaseMVA/Branch.BR_X_PU[i]
        return m.Pij[i,s]==(B)*(m.th[From_ix,s]-m.th[To_ix,s])    
    m.branchflow=pm.Constraint(m.branch_set,m.scen_set,rule=Branch_Flow)
    
    def SlackNode(model,i,s):
        if Bus.Bus_Type[i]==3:
            return m.th[i,s]==0
        else:
            return pm.Constraint.Skip
    m.slack=pm.Constraint(m.bus_set,m.scen_set,rule=SlackNode)
    
    def Power_UB_LB(model,i,s):
            Binary_LB=TruthTable[i][s] # Truth table index
            Binary_UB=1-Binary_LB
            return m.P[i,s]== Binary_UB*m.P_UB[i] + Binary_LB*m.P_LB[i]
    m.Power_UB_LB=pm.Constraint(m.gasgen_set,m.scen_set,rule=Power_UB_LB)
    
    def Costs(model,s):
        return sum(Gen.CostCoeff_1[k]*m.P[k,s] for k in m.gen_set) <= Elecdata.Params.C_max_RO
    m.Costs=pm.Constraint(m.scen_set,rule=Costs)
    
    def Shed(model,i):
        return m.P_Shed[i]==Gen.DC_OPF_RES[i]-m.P_LB[i]
    m.shed=pm.Constraint(m.gasgen_set,rule=Shed)
    
    def Add(model,i):
        return m.P_Add[i]==m.P_UB[i]-Gen.DC_OPF_RES[i]
    m.add=pm.Constraint(m.gasgen_set,rule=Add)


    def DC_OPF_Obj(model):
        return sum(Gen.AddWeight[i]*m.P_Add[i]+ Gen.ShedWeight[i]*m.P_Shed[i] for i in m.gasgen_set)
        
    m.objective = pm.Objective(rule=DC_OPF_Obj, sense=pm.maximize, doc='Define objective function') 

    
    opt = pm.SolverFactory('ipopt')
    opt.options['print_level']=5
    
    # Optimize
    print('Starting Optimization')
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


    
    
    UB_Res     = dict([[i,np.round(m.P_UB[i].value,decimals=8)]   for i in m.P_UB])
    LB_Res     = dict([[i,np.round(m.P_LB[i].value,decimals=8)]   for i in m.P_LB])
    
    P_Shed     = dict([[i,np.round(m.P_Shed[i].value,decimals=8)]   for i in m.P_Shed])
    P_Add      = dict([[i,np.round( m.P_Add[i].value,decimals=8)]   for i in m.P_Add ])
    
    Elecdata.Gen    =  Elecdata.Gen.assign(DC_RO_OPF_P_UB=pd.Series(UB_Res))
    Elecdata.Gen    =  Elecdata.Gen.assign(DC_RO_OPF_P_LB=pd.Series(LB_Res))
    
    #Elecdata.Gen    =  Elecdata.Gen.assign(P_Shed = pd.Series(P_Shed))
    #Elecdata.Gen    =  Elecdata.Gen.assign(P_Add  = pd.Series(P_Add ))
    #Elecdata.Gen    =  Elecdata.Gen.combine(pd.DataFrame.from_dict(LB_Res,orient='index',columns=['DC_RO_OPF_P_LB']))
    
    # Unload results
  
    
    Pc=pd.DataFrame(columns=[Gen.index.tolist()+['Cost','Max Cost']],index=TruthTable.index.tolist())
    for s in m.scen_set:
        for g in m.gen_set:
            Pc.loc[s][g]=m.P[g,s].value
        Pc.loc[s]['Cost']=m.Costs[s].body()
        Pc.loc[s]['Max Cost']=m.Costs[s].upper()
    
    Pc.index=['Case'+str(x) for x in range(0,len(Pc))]
    
    Pc=Pc.astype(float).round(decimals=1)
    
    Elecdata.Params.Cases=Pc
    Elecdata.Params.RO_OPF_Obj=m.objective()
    
    Elecdata.Params.status=status
    
def DC_OPF_RO_Col_Constraint_Gen(Elecdata):
    Gen=Elecdata.Gen
    GasGenerators=Gen[Gen.FuelType=='Gas'].index.tolist()
    TruthTable=[i for i in itertools.product([0,1], repeat=len(GasGenerators))]
    FullTruthTable=pd.DataFrame(TruthTable,columns=GasGenerators)
    
    
    TruthTable_iter=[len(GasGenerators)*[1],len(GasGenerators)*[0]]
    TruthTable_iter=pd.DataFrame(TruthTable_iter,columns=GasGenerators)
    DC_OPF_RO_TruthTable(Elecdata,TruthTable_iter)
#    Flag=True
#    while Flag:
#    # Solve with initial truth table
#        DC_OPF_RO_TruthTable(Elecdata,TruthTable_iter)
#        
#        
#        Cur_TruthTable=FullTruthTable.iloc[0]
#        
#        m=DC_OPF_RO_TruthTable_ConstraintCheck(Elecdata,Cur_TruthTable)
#        
#        opt = pm.SolverFactory('ipopt')
#        opt.options['print_level']=0
#        
#        Infeasible_Rows=[]
#        for r_i,row in FullTruthTable.iterrows():
#            
#            for i in row.index:
#                m.TruthTable[i]=row[i]
#            
#            results = opt.solve(m, tee=True)
#    
#            if (results.solver.status == SolverStatus.ok) and (results.solver.termination_condition == TerminationCondition.optimal):
#                print('Model {0} Solved to Optimality'.format(r_i))
#            else:
#                Infeasible_Rows.append(r_i)
#        
#        if len(Infeasible_Rows)==0:
#            Flag=False
#        else:
#            for new_row in Infeasible_Rows:
#                TruthTable_iter.iloc[len(TruthTable_iter)+1]=FullTruthTable.iloc[new_row]


def DC_OPF_RO_TruthTable_ConstraintCheck(Elecdata,Cur_TruthTable):
    
    Gen=Elecdata.Gen
    Branch=Elecdata.Branch
    Bus=Elecdata.Bus
    GasGenerators=Gen[Gen.FuelType=='Gas'].index.tolist()
    
    m = pm.ConcreteModel()
    
    m.gen_set    = pm.Set(initialize=Gen.index.tolist())
    m.branch_set = pm.Set(initialize=Branch.index.tolist())
    m.bus_set    = pm.Set(initialize=Bus.index.tolist())
    
    
    m.gasgen_set = pm.Set(initialize=GasGenerators,within=m.gen_set,ordered=True)
    
    m.TruthTable = pm.Param(m.gasgen_set,mutable=True, initialize=lambda m,i : (Cur_TruthTable[i]))
    m.P_UB = pm.Param(m.gasgen_set,mutable=False, initialize= lambda m,i : (Gen.DC_RO_OPF_P_UB[i]))
    m.P_LB = pm.Param(m.gasgen_set,mutable=False, initialize= lambda m,i : (Gen.DC_RO_OPF_P_LB[i]))
    #m.P_UB = pm.Var(m.gasgen_set,bounds= lambda m,i : (Gen.DC_RO_OPF_P_UB[i],Gen.DC_RO_OPF_P_UB[i]),initialize=lambda m,i : Gen.DC_OPF_RES[i] )
    #m.P_LB = pm.Var(m.gasgen_set,bounds= lambda m,i : (Gen.Pmin_MW[i],Gen.DC_OPF_RES[i]),initialize=lambda m,i : Gen.DC_OPF_RES[i] )
    m.P    = pm.Var(m.gen_set    ,bounds = lambda m,i : (Gen.Pmin_MW[i], Gen.Pmax_MW[i]), initialize=1)
    m.Pij  = pm.Var(m.branch_set ,bounds = lambda m,i : (-Branch.RateA_MVA[i], Branch.RateA_MVA[i]), initialize=1)
    m.th   = pm.Var(m.bus_set   ,bounds = (-np.pi,np.pi),initialize=0)
    m.P_Shed = pm.Var(m.gasgen_set,bounds=(0,None))
    m.P_Add  = pm.Var(m.gasgen_set,bounds=(0,None))
    

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
    
    def SlackNode(model,i):
        if Bus.Bus_Type[i]==3:
            return m.th[i]==0
        else:
            return pm.Constraint.Skip
    m.slack=pm.Constraint(m.bus_set,rule=SlackNode)
    
    def Power_UB_LB(model,i):
            Binary_LB=m.TruthTable[i] # Truth table index
            Binary_UB=1-Binary_LB
            return m.P[i]== Binary_UB*m.P_UB[i] + Binary_LB*m.P_LB[i]
    m.Power_UB_LB=pm.Constraint(m.gasgen_set,rule=Power_UB_LB)
    
    def Costs(model):
        return sum(Gen.CostCoeff_1[k]*m.P[k] for k in m.gen_set) <= Elecdata.Params.C_max_RO
    m.Costs=pm.Constraint(rule=Costs)
    
    def Shed(model,i):
        return m.P_Shed[i]==Gen.DC_OPF_RES[i]-m.P_LB[i]
    m.shed=pm.Constraint(m.gasgen_set,rule=Shed)
    
    def Add(model,i):
        return m.P_Add[i]==m.P_UB[i]-Gen.DC_OPF_RES[i]
    m.add=pm.Constraint(m.gasgen_set,rule=Add)


    def DC_OPF_Obj(model):
        return sum(Gen.AddWeight[i]*m.P_Add[i]+ Gen.ShedWeight[i]*m.P_Shed[i] for i in m.gasgen_set)
        
    m.objective = pm.Objective(rule=DC_OPF_Obj, sense=pm.maximize, doc='Define objective function') 

    return m


def DC_OPF_RO_TruthTable(Elecdata,TruthTable):
    
    Gen=Elecdata.Gen
    Branch=Elecdata.Branch
    Bus=Elecdata.Bus

    GasGenerators=Gen[Gen.FuelType=='Gas'].index.tolist()
    

    
    nscen=len(TruthTable)
    
    m = pm.ConcreteModel()
    
    m.gen_set    = pm.Set(initialize=Gen.index.tolist())
    m.branch_set = pm.Set(initialize=Branch.index.tolist())
    m.bus_set    = pm.Set(initialize=Bus.index.tolist())
    
    
    m.scen_set   = pm.RangeSet(0,nscen-1)
    m.gasgen_set = pm.Set(initialize=GasGenerators,within=m.gen_set,ordered=True)
    
    m.P_UB = pm.Var(m.gasgen_set,bounds= lambda m,i : (Gen.DC_OPF_RES[i],Gen.Pmax_MW[i]),initialize=lambda m,i : Gen.DC_OPF_RES[i] )
    m.P_LB = pm.Var(m.gasgen_set,bounds= lambda m,i : (Gen.Pmin_MW[i],Gen.DC_OPF_RES[i]),initialize=lambda m,i : Gen.DC_OPF_RES[i] )
    m.P    = pm.Var(m.gen_set    ,m.scen_set,bounds = lambda m,i,s : (Gen.Pmin_MW[i], Gen.Pmax_MW[i]), initialize=1)
    m.Pij  = pm.Var(m.branch_set ,m.scen_set,bounds = lambda m,i,s : (-Branch.RateA_MVA[i], Branch.RateA_MVA[i]), initialize=1)
    m.th   = pm.Var(m.bus_set   ,m.scen_set,bounds = (-np.pi,np.pi),initialize=0)
    m.P_Shed = pm.Var(m.gasgen_set,bounds=(0,None))
    m.P_Add  = pm.Var(m.gasgen_set,bounds=(0,None))
    

    def PowerBal_constr(model,i,s):
        PowerBal = -Bus.PD_MW[i] \
        + sum(m.P[k,s]   for k in Gen[Gen.Gen_Bus==i].index.tolist()) \
        - sum(m.Pij[k,s] for k in Branch[Branch.From_Bus==i].index.tolist()) \
        + sum(m.Pij[k,s] for k in Branch[Branch.To_Bus==i].index.tolist()) \
        ==0
        return PowerBal          
    m.PowerBal_constr=pm.Constraint(m.bus_set,m.scen_set,rule=PowerBal_constr)
    
    
    
    def Branch_Flow(model,i,s):
        From_ix=Branch.From_Bus[i]
        To_ix=Branch.To_Bus[i]
        B=Elecdata.Params.BaseMVA/Branch.BR_X_PU[i]
        return m.Pij[i,s]==(B)*(m.th[From_ix,s]-m.th[To_ix,s])    
    m.branchflow=pm.Constraint(m.branch_set,m.scen_set,rule=Branch_Flow)
    
    def SlackNode(model,i,s):
        if Bus.Bus_Type[i]==3:
            return m.th[i,s]==0
        else:
            return pm.Constraint.Skip
    m.slack=pm.Constraint(m.bus_set,m.scen_set,rule=SlackNode)
    
    def Power_UB_LB(model,i,s):
            Binary_LB=TruthTable[i][s] # Truth table index
            Binary_UB=1-Binary_LB
            return m.P[i,s]== Binary_UB*m.P_UB[i] + Binary_LB*m.P_LB[i]
    m.Power_UB_LB=pm.Constraint(m.gasgen_set,m.scen_set,rule=Power_UB_LB)
    
    def Costs(model,s):
        return sum(Gen.CostCoeff_1[k]*m.P[k,s] for k in m.gen_set) <= Elecdata.Params.C_max_RO
    m.Costs=pm.Constraint(m.scen_set,rule=Costs)
    
    def Shed(model,i):
        return m.P_Shed[i]==Gen.DC_OPF_RES[i]-m.P_LB[i]
    m.shed=pm.Constraint(m.gasgen_set,rule=Shed)
    
    def Add(model,i):
        return m.P_Add[i]==m.P_UB[i]-Gen.DC_OPF_RES[i]
    m.add=pm.Constraint(m.gasgen_set,rule=Add)


    def DC_OPF_Obj(model):
        return sum(Gen.AddWeight[i]*m.P_Add[i]+ Gen.ShedWeight[i]*m.P_Shed[i] for i in m.gasgen_set)
        
    m.objective = pm.Objective(rule=DC_OPF_Obj, sense=pm.maximize, doc='Define objective function') 

    
    opt = pm.SolverFactory('ipopt')
    opt.options['print_level']=5
    
    # Optimize
    print('Starting Optimization')
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


    
    
    UB_Res     = dict([[i,np.round(m.P_UB[i].value,decimals=8)]   for i in m.P_UB])
    LB_Res     = dict([[i,np.round(m.P_LB[i].value,decimals=8)]   for i in m.P_LB])
    
    P_Shed     = dict([[i,np.round(m.P_Shed[i].value,decimals=8)]   for i in m.P_Shed])
    P_Add      = dict([[i,np.round( m.P_Add[i].value,decimals=8)]   for i in m.P_Add ])
    
    Elecdata.Gen    =  Elecdata.Gen.assign(DC_RO_OPF_P_UB=pd.Series(UB_Res))
    Elecdata.Gen    =  Elecdata.Gen.assign(DC_RO_OPF_P_LB=pd.Series(LB_Res))
    
    Elecdata.Gen    =  Elecdata.Gen.assign(P_Shed = pd.Series(P_Shed))
    Elecdata.Gen    =  Elecdata.Gen.assign(P_Add  = pd.Series(P_Add ))
    #Elecdata.Gen    =  Elecdata.Gen.combine(pd.DataFrame.from_dict(LB_Res,orient='index',columns=['DC_RO_OPF_P_LB']))
    
    # Unload results
  
    
    Pc=pd.DataFrame(columns=[Gen.index.tolist()+['Cost','Max Cost']],index=TruthTable.index.tolist())
    for s in m.scen_set:
        for g in m.gen_set:
            Pc.loc[s][g]=m.P[g,s].value
        Pc.loc[s]['Cost']=m.Costs[s].body()
        Pc.loc[s]['Max Cost']=m.Costs[s].upper()
    
    Pc.index=['Case'+str(x) for x in range(0,len(Pc))]
    
    Pc=Pc.astype(float).round(decimals=1)
    
    Elecdata.Params.Cases=Pc
    Elecdata.Params.RO_OPF_Obj=m.objective()
    
    Elecdata.Params.status=status