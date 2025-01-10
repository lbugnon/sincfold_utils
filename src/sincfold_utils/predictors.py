"""Classical secondary structure predictors. The predictors need to be installed separately."""
import os
from datetime import datetime
import shutil 
from sincfold_utils.io import dot2bp

def linearpartition(seq, lppath, mode="C", convert_to_bp=True):
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
    
    fasta_name = f"{fname}.fasta"
    with open(fasta_name, "w") as ofile: 
        ofile.write(f">id\n{seq}\n")

    #==========================================
    if mode=="C":
        # CONTRAFOLD
        os.system(f'cat {fname}.fasta | {lppath} -M > {fname}.dot 2>{fname}.log')
    elif mode=="V":
        # VIENNA
        os.system(f'cat {fname}.fasta | {lppath} -V -M > {fname}.dot 2>{fname}.log')
    #==========================================

    # Reading prediction
    with open(f'{fname}.dot', 'r') as fp:
        lines = fp.read().strip().split('\n')
        prediction = lines[-1]

    lpscore = float(open(f"{fname}.log").read().strip().split(":")[-1].replace("kcal/mol", ""))
    
    # Removing files
    os.remove(f'{fname}.fasta')
    os.remove(f'{fname}.dot')
    os.remove(f'{fname}.log')
    
    if convert_to_bp:
        prediction = dot2bp(prediction)
    
    return prediction, lpscore