import os
datapath = os.path.join(os.path.dirname(__file__),'data_tables')
os.environ["DATAPATH"] = datapath
from ._RNAstructure import RNA,Dynalign_object,Multilign_object #,HybridRNA,Dynalign_object,
