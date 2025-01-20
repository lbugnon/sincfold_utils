"""Classical secondary structure predictors. The predictors need to be installed separately."""
import os
from datetime import datetime
import shutil 
from sincfold_utils.io import dot2bp

def linearpartition(fname, mode="V", lppath="../LinearPartition/linearpartition"):
    
    if mode=="C":
        # CONTRAFOLD
        os.system(f'cat {fname}.fasta | {lppath} -M > {fname}.dot 2>{fname}.log')
    elif mode=="V":
        # VIENNA
        os.system(f'cat {fname}.fasta | {lppath} -V -M > {fname}.dot 2>{fname}.log')
            
    # Reading prediction
    lines = open(f"{fname}.dot").readlines()
    prediction = lines[2].strip()
    score = float(open(f"{fname}.log").read().strip().split(":")[-1].replace("kcal/mol", ""))

    return prediction, score

def rnafold(fname, temp=37):

    os.system(f"RNAfold -T {temp} --noPS {fname}.fasta > {fname}.dot")
    lines = open(f"{fname}.dot").readlines()
    L = len(lines[1].strip())
    prediction = lines[2].strip()
    score = float(prediction[L:].strip(" ()"))
    prediction = prediction[:L]
    
    return prediction, score

def rnastructure(fname, temp=298.15, source_path="../RNAstructure/", install_path="/usr/local/RNAstructure/"):
    # Compute structure
    os.system(f"export DATAPATH={source_path}data_tables; {install_path}Fold \
        --bracket --MFE -t {temp} {fname}.fasta {fname}.dot")
    
    lines = open(f"{fname}.dot").readlines()
    score = float(lines[0].split(" ")[2]) # energy
    structure = lines[-1].strip()
    
    return structure, score

def fold(seq, method, args={}, convert_to_bp=True):
    """
    seq : str
        RNA sequence
    lppath : str
        Path to the LinearPartition executable
    mode : str
        C or V.
        
    Output:
        lpscore: Free Energy of Ensemble if "V", Log Partition Coefficient if "C" 
    """
    # Writing sequence to file
    fname = str(datetime.now()).replace(" ","_")
    method = method.lower()
    
    with open(f"{fname}.fasta", "w") as ofile: 
        ofile.write(f">{fname}\n{seq}")

    #==========================================
    if method == "rnastructure":
        prediction, score = rnastructure(fname, **args)
    elif method == "linearpartition":
        prediction, score = linearpartition(fname, **args)
    elif method == "rnafold":
        prediction, score = rnafold(fname, **args)
        
    # Removing files
    os.remove(f'{fname}.fasta')
    os.remove(f'{fname}.dot')
    if os.path.isfile(f'{fname}.log'):
        os.remove(f'{fname}.log')
    
    if convert_to_bp:
        prediction = dot2bp(prediction)
    
    return prediction, score