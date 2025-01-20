import subprocess
import numpy as np
from sklearn.metrics import f1_score, precision_score, recall_score
from grakel import Graph
from grakel.kernels import WeisfeilerLehman, VertexHistogram, WeisfeilerLehmanOptimalAssignment, ShortestPath

from sincfold_utils.constants import MATCHING_BRACKETS


def get_graph_kernel(kernel, n_iter=5, normalize=True):
    if kernel == 'WeisfeilerLehman':
        return WeisfeilerLehman(n_iter=n_iter,
                                normalize=normalize,
                                base_graph_kernel=VertexHistogram)
    elif kernel == 'WeisfeilerLehmanOptimalAssignment':
        return WeisfeilerLehmanOptimalAssignment(n_iter=n_iter,
                                                 normalize=normalize)
    elif kernel == 'ShortestPath':
        return ShortestPath(normalize=normalize)
    
def mat2graph(matrix, node_labels=None):
    if node_labels is None:
        node_labels = {s: str(s) for s in range(matrix.shape[0])}
    graph = Graph(initialization_object=matrix.astype(int), node_labels=node_labels)  
    return graph

def WL(ref, pred, sequence=None, kernel="WeisfeilerLehman",  n_iter=5):
    """pred, ref are 2D binary matrices"""
    
    node_labels = None
    if sequence is not None:
        node_labels = {i: s for i, s in enumerate(sequence)}
    
    pred_graph = mat2graph(pred, node_labels=node_labels)
    true_graph = mat2graph(ref, node_labels=node_labels)
    kernel = get_graph_kernel(kernel=kernel, n_iter=n_iter)

    kernel.fit_transform([true_graph])
    distance_score = kernel.transform([pred_graph])  

    return distance_score[0][0]

def triangular_metrics(ref, pred):
    """Compute F1, recall and precision from the upper triangular connection matrix. ref and pred are binary 2D numpy arrays"""
    assert type(ref) == np.ndarray, "ref must be a numpy array"
    assert type(pred) == np.ndarray, "pred must be a numpy array"
    assert ref.shape == pred.shape, "ref and pred must have the same shape"
    assert ref.shape[0] == ref.shape[1], "ref and pred must be square matrices"
    assert np.all(np.logical_or(ref == 0, ref == 1)), "ref must be binary"
    assert np.all(np.logical_or(pred == 0, pred == 1)), "pred must be binary"

    # get upper triangular matrix without diagonal
    ind = np.triu_indices(ref.shape[0], k=1) # k is the diagonal offset

    ref = ref[ind[0], ind[1]].ravel()
    pred = pred[ind[0], ind[1]].ravel()

    rec = recall_score(ref, pred)
    pre = precision_score(ref, pred)
    f1 = f1_score(ref, pred, zero_division=0)
    assert 2*rec*pre/(rec+pre) == f1, "F1 score is not consistent with precision and recall"
    return f1, rec, pre 

def get_tp(bp_x, bp_ref, strict):
    tp = 0
    for rbp in bp_x:
        cond = rbp in bp_ref
        if not strict:
            cond = cond or [rbp[0], rbp[1] - 1] in bp_ref or [rbp[0], rbp[1] + 1] in bp_ref or [rbp[0] + 1, rbp[1]] in bp_ref or [rbp[0] - 1, rbp[1]] in bp_ref
        if cond:
            tp += 1
    return tp   

# TODO requires validation, why using tp1 and tp2 in precision?
def bp_metrics(ref_bp, pre_bp, strict=True):
    """F1, recall and precision score from base pairs. Is the same as triangular but less efficient. strict=False takes into account the  neighbors for each base pair as correct"""
    assert all([type(bp) == list for bp in ref_bp]), "ref_bp must be a list of lists"
    assert all([type(bp) == list for bp in pre_bp]), "pre_bp must be a list of lists"

    # corner case when there are no positives
    if len(ref_bp) == 0 and len(pre_bp) == 0:
        return 1.0, 1.0, 1.0

    tp1 = get_tp(pre_bp, ref_bp, strict)
    tp2 = get_tp(ref_bp, pre_bp, strict)

    fn = len(ref_bp) - tp1
    fp = len(pre_bp) - tp1

    tpr = pre = f1 = 0.0
    if tp1 + fn > 0:
        tpr = tp1 / float(tp1 + fn)  # sensitivity (=recall =power)
    if tp1 + fp > 0:
        pre = tp2 / float(tp1 + fp)  # precision (=ppv)
    if tpr + pre > 0:
        f1 = 2 * pre * tpr / (pre + tpr)  # F1 score

    return f1, tpr, pre

def rnadistance(struct1, struct2):
    """Dotbracket structures to compare"""  
    assert len(struct1)==len(struct2), "structures must have the same lenght"
    
    echo_line = struct1 + "\n" + struct2
    # RNAdistance cannot handle pseudoknots (treated as unpaired) https://pubmed.ncbi.nlm.nih.gov/36077037/
    for b in MATCHING_BRACKETS[1:]:
        echo_line = echo_line.replace(b[0],'.').replace(b[1],'.')
    oo = subprocess.check_output([f"RNAdistance"], input=echo_line.encode('utf-8'))

    if oo == b"":
        raise ValueError("RNAdistance failed")
    
    return float(oo[2:])/len(struct1)
    