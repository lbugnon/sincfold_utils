import os
import tempfile
import numpy as np

from sincfold_utils import __path__ as sincfold_utils_path
import subprocess as sp
from sincfold_utils.constants import STANDARD_BRACKETS
from sincfold_utils.structure import validate_pairs, pairs2db, db2pairs

CT2DOT_CALL = f"export DATAPATH={sincfold_utils_path[0]}/tools/RNAstructure/data_tables; {sincfold_utils_path[0]}/tools/RNAstructure/ct2dot"


def _ct_header_length(line):
    """Sequence length from a ct header: 'N title', 'N ENERGY = x title', 'N dG = x title' or '>seq length: N ...'"""
    tokens = line.split()
    if tokens and tokens[0].isdigit():
        return int(tokens[0])
    if line.split(":")[0].strip() == ">seq length":
        return int(line.split(":")[1].split()[0])
    raise ValueError(f"ct file format error: invalid header '{line.strip()}'")


def _ct_row(line):
    """(index, nucleotide, partner) if line is a ct data row, else None"""
    cols = line.split()
    if len(cols) < 6 or not all(c.lstrip("-").isdigit() for c in [cols[0]] + cols[2:6]):
        return None
    return int(cols[0]), cols[1], int(cols[4])


def read_ct(ctfile, multi=False):
    """Read a ct file, return sequence and base pairs (1-based, i<j).

    Columns are whitespace separated (tabs or spaces). Only columns 1 (index), 2 (nucleotide) and 5 (partner)
    are used, the natural numbering in column 6 is ignored. Gaps in the numbering are filled with N.
    Pairing must be symmetric, otherwise ValueError is raised.
    With multi=True, all the structures in the file are returned as a list of (title, sequence, base_pairs)."""
    lines = [l for l in open(ctfile).read().splitlines() if l.strip() and not l.lstrip().startswith("#")]
    blocks, k = [], 0
    while k < len(lines):
        header = lines[k]
        seq_len = _ct_header_length(header)
        k += 1
        seq, partner = [], {}
        while k < len(lines) and len(seq) < seq_len:
            row = _ct_row(lines[k])
            if row is None:
                break
            n1, nt, n2 = row
            if n1 <= len(seq) or n1 > seq_len:
                raise ValueError(f"{ctfile}: nucleotide index {n1} out of order or range")
            seq += ["N"] * (n1 - len(seq) - 1)
            seq.append(nt)
            if n2 != 0:
                partner[n1] = n2
            k += 1

        if len(seq) != seq_len:
            raise ValueError(f"{ctfile}: header length {seq_len} but {len(seq)} nucleotides")
        for i, j in partner.items():
            if j < 1 or j > seq_len or j == i or partner.get(j) != i:
                raise ValueError(f"{ctfile}: inconsistent pairing {i}-{j}")
        bp = [[i, j] for i, j in sorted(partner.items()) if i < j]
        blocks.append((" ".join(header.split()[1:]), "".join(seq), bp))

    if len(blocks) == 0:
        raise ValueError(f"{ctfile}: empty ct file")
    if multi:
        return blocks
    return blocks[0][1], blocks[0][2]


def write_ct(fname, seqid, seq, base_pairs):
    """Write ct file from sequence and base pairs. Base_pairs should be 1-based and unique per nt"""
    base_pairs = validate_pairs(base_pairs, len(seq))
    base_pairs_dict = {}
    for i, j in base_pairs:
        base_pairs_dict[i] = j
        base_pairs_dict[j] = i

    L = len(seq)
    with open(fname, "w") as fout:
        fout.write(f"{L} {seqid}\n")
        for k, n in enumerate(seq, start=1):
            nxt = k + 1 if k < L else 0
            fout.write(f"{k} {n} {k-1} {nxt} {base_pairs_dict.get(k, 0)} {k}\n")


def read_seq(seqfile):
    """Read a .seq file (RNAstructure / ArchiveII format): ';' comment lines, a title line and the
    sequence terminated by '1'. Returns title and sequence (case is preserved)."""
    lines = [l.strip() for l in open(seqfile) if not l.startswith(";") and l.strip()]
    title, seq = lines[0], "".join(lines[1:]).replace(" ", "")
    if not seq.endswith("1"):
        raise ValueError(f"{seqfile}: sequence must end with '1'")
    return title, seq[:-1]


def mat2dot(mat, brackets=STANDARD_BRACKETS):
    """Dot-bracket from a binary (upper triangular or symmetric) connection matrix"""
    ind = np.argwhere(np.triu(mat, k=1) == 1) + 1
    return pairs2db(ind.tolist(), mat.shape[0], brackets)


def bp2mat(base_pairs, L):
    """base_pairs: 1-indexed base pairs"""
    mat = np.zeros((L, L))
    for bp in base_pairs:
        mat[bp[0]-1, bp[1]-1] = 1
        mat[bp[1]-1, bp[0]-1] = 1

    return mat


def bp2dot(seq, bp, brackets=STANDARD_BRACKETS):
    """Convert base pairs to dot-bracket notation (pure python, supports pseudoknots).
    seq can be the sequence or its length."""
    L = seq if isinstance(seq, int) else len(seq)
    return pairs2db(bp, L, brackets)


def dot2bp(struct):
    """Base pairs (1-based, sorted) from a dot-bracket string with any bracket type ()[]{}<>AaBb...
    Raises ValueError for invalid structures."""
    return db2pairs(struct)


def ct2dot(ct_file, structure=1, brackets=STANDARD_BRACKETS):
    """Dot-bracket from a ct file (pure python, supports pseudoknots). structure is the 1-based index of the
    structure in multi-structure ct files. Use brackets=RNASTRUCTURE_BRACKETS to get the RNAstructure ct2dot
    bracket order ( < { [ A B ..."""
    blocks = read_ct(ct_file, multi=True)
    if not 1 <= structure <= len(blocks):
        raise ValueError(f"{ct_file} has {len(blocks)} structures, requested {structure}")
    _, seq, bp = blocks[structure - 1]
    return pairs2db(bp, len(seq), brackets)


def ct2dot_rnastructure(ct_file, structure=1):
    """Dot-bracket using the bundled RNAstructure ct2dot binary (bracket order ( < { [ A B ...).
    Kept for comparison with RNAstructure, ct2dot does not need it."""
    if not os.path.isfile(ct_file) or os.path.splitext(ct_file)[1] != ".ct":
        raise ValueError("ct2dot requires a .ct file")
    with tempfile.TemporaryDirectory() as tmp:
        out = os.path.join(tmp, "out.dot")
        res = sp.run(f"{CT2DOT_CALL} {ct_file} {structure} {out}", shell=True, capture_output=True)
        lines = open(out).read().splitlines() if os.path.isfile(out) else []
    if len(lines) < 3 or len(lines[2].strip()) == 0:
        raise RuntimeError(f"ct2dot failed on {ct_file}: {res.stdout.decode()} {res.stderr.decode()}")
    return lines[2].strip()
