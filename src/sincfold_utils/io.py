import os 
from sincfold_utils import __path__ as sincfold_utils_path
import subprocess as sp 

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


# TODO it would be nice to have an standalone function
def bp2dot(seq, bp):
    write_ct("tmp.ct", "seq", seq, bp)
    dotbracket = ct2dot("tmp.ct")
    os.remove("tmp.ct")
    return dotbracket

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

