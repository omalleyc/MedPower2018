import pandas as pd
import numpy as np


class expando(object):
    pass

def Load_CaseGB():
    # 2 Elec Bus and 2 Pipe with compressor
    filename_e='GB_Elec'   
    filename_g='GB_Gas'
    filename_c='Coupling/Coupling_GB_Reduced'
    
    
    Elecdata=ElecData(filename_e)
    Gasdata=GasData(filename_g)
    Gasdata.Params.Initialize_from='Flat_Start'
    
    # Increase the gas load until the system is stressed
    Gasdata.Nodes.Load=2.15*Gasdata.Nodes.Load
    CouplingData(filename_c,Elecdata,Gasdata)
    
    return Elecdata,Gasdata

def Load_Case1():
    # 2 Elec Bus and 2 Pipe with compressor
    filename_e='Bus2'   
    filename_g='Gas2NodeCompressor'
    filename_c='Coupling/Coupling_Case1'
    
    
    Elecdata=ElecData(filename_e)
    Gasdata=GasData(filename_g)
    Gasdata.Params.Initialize_from='Flat_Start'
    CouplingData(filename_c,Elecdata,Gasdata)
    
    return Elecdata,Gasdata

def Load_Case2():
    filename_e='Bus3'
    filename_g='Gas3NodeCompressor'
    filename_c='Coupling/Coupling_Case2'
    
    Elecdata=ElecData(filename_e)
    Gasdata=GasData(filename_g)
    Gasdata.Params.Initialize_from='Flat_Start'
    CouplingData(filename_c,Elecdata,Gasdata)
    
    return Elecdata,Gasdata


def Load_Case2_v():
    # Variation on Case2 with altered loads
    filename_e='Bus3'
    filename_g='Gas3NodeCompressor'
    filename_c='Coupling/Coupling_Case2'
    
    Elecdata=ElecData(filename_e)
    Gasdata=GasData(filename_g)
    Gasdata.Params.Initialize_from='Flat_Start'
    CouplingData(filename_c,Elecdata,Gasdata)
    
    #Elecdata.Branch.RateA_MVA=1e3
    Gasdata.Nodes.loc['Node4','Load']=4.3/Gasdata.Params.NormalizedM
    Gasdata.Nodes.loc['Node5','Load']=7/Gasdata.Params.NormalizedM
    return Elecdata,Gasdata


def Load_Case5():
    filename_e='case5' 
    filename_g='GasRadial'
    filename_c='Coupling/Coupling_Radial'
    
    Elecdata=ElecData(filename_e)
    Gasdata=GasData(filename_g)
    Gasdata.Params.Initialize_from='Flat_Start'
    CouplingData(filename_c,Elecdata,Gasdata)
    
    Elecdata.Bus['PD_MW']=0.75* Elecdata.Bus['PD_MW']


    return Elecdata,Gasdata

def Load_Case5_v():
    # Variation on Case 5 with alterered parameters
    filename_e='case5' 
    filename_g='GasRadial'
    filename_c='Coupling/Coupling_Radial'
    
    Elecdata=ElecData(filename_e)
    Gasdata=GasData(filename_g)
    Gasdata.Params.Initialize_from='Flat_Start'
    CouplingData(filename_c,Elecdata,Gasdata)
    
    # Adjustment to Original
    Elecdata.Gen.loc['Gen5','CostCoeff_1']=50
    #
    Elecdata.Gen.Pmax_MW=[1000,100,100,50,1000]
    #
    Elecdata.Branch.RateA_MVA=1e6
    
    Elecdata.Bus.PD_MW=0
    
    Elecdata.Bus.loc['Bus5','PD_MW']=1140

    return Elecdata,Gasdata



def ElecData(filename):
    # Read from Files that are created from Matpower cases
    
    Gen=pd.read_csv(filename+'/'+filename+'_gen.csv')
    Bus=pd.read_csv(filename+'/'+filename+'_bus.csv')
    Branch=pd.read_csv(filename+'/'+filename+'_branch.csv')
    
    Params=expando()
    
    with open(filename+'/'+filename+'_other.csv', 'r') as f:
        for line in f:
            data=line.strip().split()
            Name=data[0]
            Value=data[1]
            try:
                setattr(Params,Name , float(Value))
            except ValueError:
                setattr(Params,Name , (Value))
    
    # Convert Numbered components to names to avoid confusion later
    Gen.index=['Gen'+str(x) for x in range(1,len(Gen)+1)]
    Gen.Gen_Bus=['Bus'+str(x) for x in Gen.Gen_Bus.tolist()]
    
    Bus.index=['Bus'+str(x) for x in Bus.Bus_Number.tolist()]
    Bus.drop(['Bus_Number'],axis=1,inplace=True)
    
    Branch.index=['Branch'+str(x) for x in range(1,len(Branch)+1)]
    Branch.From_Bus=['Bus'+str(x) for x in Branch.From_Bus.tolist()]
    Branch.To_Bus=['Bus'+str(x) for x in Branch.To_Bus.tolist()]
    
    Branch.RateA_MVA.replace(0,1e6,inplace=True)
    Branch.RateB_MVA.replace(0,1e6,inplace=True)
    Branch.RateC_MVA.replace(0,1e6,inplace=True)
    
    Elecdata=expando()
    
    Elecdata.Bus=Bus
    Elecdata.Gen=Gen
    Elecdata.Branch=Branch
    Elecdata.Params=Params
    
    # Drop the Unused Parameters
    To_be_dropped=['QG_MVar','QMax_MVar','VG_PU','PG_MW',
                       'Qmin_MVar','MBase_MVA',
                       'StatusONOFF','PC1_MW',
                       'PC2_MW','QC1MIN_MVAR',
                       'QC1MAX_MVAR','QC2MIN_MVAR',
                       'QC2MAX_MVAR','RAMP_AGC_MW_MIN',
                       'RAMP_10_MW','RAMP_30_MW',
                       'RAMP_Q_MVAR_MIN',
                       'AreaParticipationFactor',
                       'CostModel','StartupCost_USD',
                       'Shutdown_USD','No_CostParameters'
                       ]
    for x in To_be_dropped:
        if x in Elecdata.Gen.columns :
            Elecdata.Gen.drop(x,axis=1,inplace=True)
    
    Elecdata.Branch.drop(['BR_R_PU','BR_B_PU',
                          'RateB_MVA','RateC_MVA',
                          'Tap','Shift',
                          'BR_Status','ANGMIN',
                          'ANGMAX'],axis=1,inplace=True)
    
    Elecdata.Bus.drop(['QD_MVar','GS_MW',
                       'BS_MVar','BusArea',
                       'Vm','Va','BaseKV',
                       'LossZone'],axis=1,inplace=True)
    return Elecdata


def GasData(filename):
    
   # Folder='Gas3NodeCompressor'
    Pipes=pd.read_csv(filename+'/Pipes.csv')
    Nodes=pd.read_csv(filename+'/Nodes.csv')
    Compressors=pd.read_csv(filename+'/Compressors.csv')

    Params=expando()
    
    with open(filename+'/Parameters.txt', 'r') as f:
        for line in f:
            data=line.strip().split()
            Name=data[0]
            Value=data[1]
            try:
                setattr(Params,Name , float(Value))
            except ValueError:
                setattr(Params,Name , (Value))
    

    # Convert Bus Numbers to strings to avoid confusion later
    Pipes.index=['Pipe'+str(x) for x in range(1,len(Pipes)+1)]
    Pipes.From=['Node'+str(x) for x in Pipes.From.tolist()]
    Pipes.To  =['Node'+str(x) for x in Pipes.To.tolist()]
    
    Nodes.index=['Node'+str(x) for x in Nodes.Node.tolist()]
    Nodes.drop(['Node'],axis=1,inplace=True)
    
    Compressors.index=['Compressor'+str(x) for x in range(1,len(Compressors)+1)]
    Compressors.From =['Node'+str(x) for x in Compressors.From.tolist()]
    Compressors.To   =['Node'+str(x) for x in Compressors.To.tolist()]
        
    # Convert to Per Unit System for less numerical difficulties
    """
    P = P0 * pi 
    Pressure [MPa] = Pressure Base [MPA] * per unit pressure
    
    M = M0 * mij 
    Mass [kg/s] = Mass Base [kg/s] * per unit mass
    
    Pi^2- Pj^2 = (L*f*c^2 /2*D*A*A )*Mij*abs(Mij)
    
    Pi^2- Pj^2 = (K)*Mij*abs(Mij)
    
    pi^2- pj^2 = (M0^2/P0^2)(K)*mij*abs(mij)
    
    
    """
    Pipes['Area']=np.pi*(Pipes.Diameter/2.0)**2
    Num=Pipes.Length*Pipes.Friction*Params.Speedofsound**2
    Den=(2*Pipes.Diameter*Pipes.Area*Pipes.Area)
    PerUnit=(Params.NormalizedM/Params.NormalizedP)**2
    Pipes['SS_Constant']=PerUnit* Num/Den
    
    Nodes.Load=Nodes.Load/Params.NormalizedM
    Params.Slackpressure=Params.Slackpressure/Params.NormalizedP
    Params.Slacknode= 'Node'+str(int(Params.Slacknode))
    
    Gasdata=expando()
    
    Gasdata.Pipes=Pipes
    Gasdata.Nodes=Nodes
    Gasdata.Compressors=Compressors
    Gasdata.Params=Params
    
    # Drop the unused Parameters
    Gasdata.Pipes.drop(['Length','Friction','Diameter','Area'],axis=1,inplace=True)
    
    return Gasdata
    
def CouplingData(filename,Elecdata,Gasdata):
    
    Coupling=pd.read_csv(filename+'.csv')
    

    
    
    Coupling.index=['Gen'+str(x) for x in Coupling.Gen_No.tolist()]
    Coupling.drop('Gen_No',axis=1,inplace=True)
    Coupling.Gas_Node=['Node'+str(x) for x in Coupling.Gas_Node.tolist()]
    
    Coupling=Coupling.assign(FuelType='Gas')
    
    Coupling.loc[:,'Power_to_Gas_Norm']=1/(Coupling.Eff*Gasdata.Params.SpecificEnergy*Gasdata.Params.NormalizedM)
    Coupling.Power_to_Gas_Norm=Coupling.Power_to_Gas_Norm.replace(np.inf,0)
    
    Elecdata.Gen=Elecdata.Gen.join(Coupling)
    
    
    return 
    
    
    
    
    
    
    
    
    
    
    
