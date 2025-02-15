import os 
import numpy as np 
from datetime import datetime

from sincfold_utils import __path__ as sincfold_utils_path
import subprocess as sp 
from sincfold_utils.constants import MATCHING_BRACKETS

CT2DOT_CALL = f"export DATAPATH={sincfold_utils_path[0]}/tools/RNAstructure/data_tables; {sincfold_utils_path[0]}/tools/RNAstructure/ct2dot"

# TODO requires validation
def read_ct(ctfile):
    """Read ct file, return sequence and base_pairs"""
    seq, bp = [], []
    
    k = 1
    for p, line in enumerate(open(ctfile)):
        if p == 0:
            try:
                seq_len = int(line.split()[0])
            except ValueError:
                # >seq length: N extra info
                if line.split(":")[0] == ">seq length":
                    seq_len = int(line.split(":")[1].split()[0])
            
            continue 

        if line[0] == "#" or len(line.strip()) == 0:
            # comment
            continue

        line = line.split()
        if len(line) != 6 or not line[0].isnumeric() or not line[4].isnumeric:
            # header
            continue

        n1, n2 = int(line[0]), int(line[4])
        if k != n1: # add missing nucleotides as N
            seq += ["N"] * (n1-k)
        seq.append(line[1])
        k = len(seq) + 1
        if n2 > 0 and (n1 < n2):
            bp.append([n1, n2])

    assert len(seq) == seq_len, f"ct file format error\n{seq_len}\n{seq}\n{len(seq)}"
    return "".join(seq), bp


def write_ct(fname, seqid, seq, base_pairs):
    """Write ct file from sequence and base pairs. Base_pairs should be 1-based and unique per nt"""
    base_pairs_dict = {}
    for bp in base_pairs:
        base_pairs_dict[bp[0]] = bp[1]
        base_pairs_dict[bp[1]] = bp[0]

    with open(fname, "w") as fout:
        fout.write(f"{len(seq)} {seqid}\n")
        for k, n in enumerate(seq):
            fout.write(f"{k+1} {n} {k} {k+2} {base_pairs_dict.get(k+1, 0)} {k+1}\n")

def mat2dot(mat):
    bp = []
    for i in range(mat.shape[0]):
        for j in range(i+1, mat.shape[1]):
            if mat[i, j] == 1:
                bp.append([i+1, j+1])
    return bp2dot(bp, mat.shape[0])

def bp2mat(base_pairs, L):
    """base_pairs: 1-indexed base pairs"""
    mat = np.zeros((L, L))
    for bp in base_pairs:
        mat[bp[0]-1, bp[1]-1] = 1
        mat[bp[1]-1, bp[0]-1] = 1
        
    return mat


def bp2dot(bp, L):
    
    fname = str(datetime.now()).replace(" ","_")
    write_ct(f"{fname}.ct", "seq", "A"*L, bp)
    dotbracket = ct2dot(f"{fname}.ct")
    os.remove(f"{fname}.ct")
    return dotbracket

def fold2bp(struc, xop="(", xcl=")"):
    """Get base pairs from one page folding (using only one type of brackets).
    BP are 1-indexed"""
    openxs = []
    bps = []
    if struc.count(xop) != struc.count(xcl):
        return False
    for i, x in enumerate(struc):
        if x == xop:
            openxs.append(i)
        elif x == xcl:
            if len(openxs) > 0:
                bps.append([openxs.pop() + 1, i + 1])
            else:
                return False
    return bps

def dot2bp(struct):
    bp = []
    if not set(struct).issubset(
        set(["."] + [c for par in MATCHING_BRACKETS for c in par])
    ):
        return False

    for brackets in MATCHING_BRACKETS:
        if brackets[0] in struct:
            bpk = fold2bp(struct, brackets[0], brackets[1])
            if bpk:
                bp = bp + bpk
            else:
                return False
    return list(sorted(bp))

def ct2dot(ct_file):
    if not os.path.isfile(ct_file) or os.path.splitext(ct_file)[1] != ".ct":
        raise ValueError("ct2dot requires a .ct file")
    dotbracket = ""
    if CT2DOT_CALL:
        res = sp.run(f"{CT2DOT_CALL} {ct_file} 1 tmp.dot", shell=True, capture_output=True)
        try: 
            dotbracket = open("tmp.dot").readlines()[2].strip()
            os.remove("tmp.dot")
        except FileNotFoundError: 
            print(res.stderr)
    else:
        print("Dotbracket conversion only available on linux")
    return dotbracket

