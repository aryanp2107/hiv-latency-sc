"""
config.py 

"""
import os
from pathlib import Path

# paths
BASE_DIR    = Path("/scratch/maid/DP2_PBMC")       
RAW_DIR     = BASE_DIR / "PBMC_rawData"            
META_PATH   = BASE_DIR / "metadata.csv"      


SAMPLE_DIR_COL = "UniqID"

#Dev mode
DEV = os.environ.get("SC_DEV", "0") == "1"
DEV_N_SAMPLES = 3
DEV_MAX_CELLS = 3000    

CHECKPOINTS = BASE_DIR / "checkpoints_dev" if DEV else "checkpoints"      
CHECKPOINTS.mkdir(parents=True, exist_ok=True)

SEED = 0       

# QC: candidate per-cell thresholds 
MIN_GENES = 200    
MAX_MITO  = 15    

# QC: whole-sample drop list 
DROP_SAMPLES = [
    "KPN_Wk0", "KPN_Wk16", "KPN_Wk72", "KPN_Wk192", "KPN_Wk240",
    "KI6_Wk0", "KI6_Wk16", "KI6_Wk72", "KI6_Wk304",
]

# mitochondrial genes 
MITO_GENES = [
    "ND1", "ND2", "ND3", "ND4", "ND4L", "ND5", "ND6",
    "COX1", "COX2", "COX3", "ATP6", "ATP8", "CYTB",
]

# downstream params
BATCH_KEY = "batch"   
N_HVG     = 2000      
N_PCS     = 30    

SAMPLE_COL = "sample"
MIN_CELLS_PER_GENE = 3
DROP_MITO_GENES = TRUE