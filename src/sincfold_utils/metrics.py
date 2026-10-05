import shutil
import subprocess
from collections import Counter
import numpy as np
from scipy.sparse.csgraph import shortest_path
from sklearn.metrics import f1_score, precision_score, recall_score

from sincfold_utils.constants import MATCHING_BRACKETS
from sincfold_utils.structure import db2pairs, pairs2db, nested_projection


def _wl_histograms(adjs, labels, n_iter):
    """Weisfeiler-Lehman relabelling of several graphs with a shared label table, so labels are comparable
    between graphs. Returns the label histogram of each graph for the initial labels and each of the n_iter
    relabelling iterations."""
    neighbors = [[np.flatnonzero(row) for row in adj] for adj in adjs]
    table = {}
    labels = [np.array([table.setdefault(l, len(table)) for l in lab]) for lab in labels]
    hists = [[np.bincount(lab, minlength=len(table)) for lab in labels]]
    for _ in range(n_iter):
        table = {}
        labels = [np.array([table.setdefault((lab[v], tuple(np.sort(lab[nb[v]]))), len(table))
                            for v in range(len(lab))], dtype=int) for lab, nb in zip(labels, neighbors)]
        hists.append([np.bincount(lab, minlength=len(table)) for lab in labels])
    return hists


def _sp_histogram(adj, labels):
    """Counts of (label u, label v, shortest path length) over all connected ordered pairs of distinct nodes"""
    dist = shortest_path(adj, unweighted=True)
    u, v = np.nonzero(np.isfinite(dist) & ~np.eye(len(adj), dtype=bool))
    return Counter((labels[a], labels[b], dist[a, b]) for a, b in zip(u, v))


def WL(ref, pred, sequence=None, kernel="WeisfeilerLehman", n_iter=5):
    """Normalized graph kernel similarity between two structures. pred, ref are 2D binary matrices (adjacency
    matrices of the base pair graphs). Nodes are labelled with the sequence if given, or with their index otherwise.

    kernel: "WeisfeilerLehman" (WL subtree kernel), "WeisfeilerLehmanOptimalAssignment" or "ShortestPath".
    n_iter: number of WL relabelling iterations (n_iter + 1 histograms, including the initial labels).
    Same results as the corresponding grakel kernels (with normalize=True), without the dependency."""
    adjs = [ref.astype(int), pred.astype(int)]
    labels = list(sequence) if sequence is not None else [str(i) for i in range(ref.shape[0])]

    if kernel == "ShortestPath":
        h_ref, h_pred = [_sp_histogram(adj, labels) for adj in adjs]
        dot = lambda a, b: sum(a[k] * b[k] for k in a.keys() & b.keys())
        k_rp, k_rr, k_pp = dot(h_ref, h_pred), dot(h_ref, h_ref), dot(h_pred, h_pred)
    elif kernel in ("WeisfeilerLehman", "WeisfeilerLehmanOptimalAssignment"):
        # WL subtree kernel: dot product of the histograms; optimal assignment: histogram intersection
        base = np.dot if kernel == "WeisfeilerLehman" else lambda a, b: np.minimum(a, b).sum()
        hists = _wl_histograms(adjs, [labels, labels], n_iter)
        k_rp, k_rr, k_pp = (sum(base(h[i], h[j]) for h in hists) for i, j in ((0, 1), (0, 0), (1, 1)))
    else:
        raise ValueError(f"unknown kernel {kernel}")

    return k_rp / np.sqrt(k_rr * k_pp)

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

    # corner case when there are no positives, same convention as bp_metrics
    if ref.sum() == 0 and pred.sum() == 0:
        return 1.0, 1.0, 1.0

    rec = recall_score(ref, pred, zero_division=0)
    pre = precision_score(ref, pred, zero_division=0)
    if pre+rec>0:
        f1 = 2*rec*pre/(pre+rec)
    else:
        f1 = 0

    return f1, rec, pre 

def get_tp(bp_x, bp_ref, strict):
    """Number of base pairs in bp_x that are found in bp_ref (or a +-1 neighbor if not strict)"""
    bp_ref = set(map(tuple, bp_ref))
    tp = 0
    for i, j in bp_x:
        cond = (i, j) in bp_ref
        if not strict:
            cond = cond or (i, j - 1) in bp_ref or (i, j + 1) in bp_ref or (i + 1, j) in bp_ref or (i - 1, j) in bp_ref
        if cond:
            tp += 1
    return tp


def bp_metrics(ref_bp, pre_bp, strict=True):
    """F1, recall and precision score from base pairs. Is the same as triangular but less efficient. strict=False takes into account the  neighbors for each base pair as correct"""
    assert all([type(bp) == list for bp in ref_bp]), "ref_bp must be a list of lists"
    assert all([type(bp) == list for bp in pre_bp]), "pre_bp must be a list of lists"

    # corner case when there are no positives
    if len(ref_bp) == 0 and len(pre_bp) == 0:
        return 1.0, 1.0, 1.0

    tp_ref = get_tp(ref_bp, pre_bp, strict)  # reference pairs recovered by the prediction
    tp_pre = get_tp(pre_bp, ref_bp, strict)  # predicted pairs found in the reference

    tpr = pre = f1 = 0.0
    if len(ref_bp) > 0:
        tpr = tp_ref / float(len(ref_bp))  # sensitivity (=recall =power)
    if len(pre_bp) > 0:
        pre = tp_pre / float(len(pre_bp))  # precision (=ppv)
    if tpr + pre > 0:
        f1 = 2 * pre * tpr / (pre + tpr)  # F1 score

    return f1, tpr, pre


def bp_distance(ref_bp, pre_bp):
    """Base pair distance: number of base pairs in only one of the structures. Includes pseudoknots."""
    return len(set(map(tuple, ref_bp)) ^ set(map(tuple, pre_bp)))


def _nested_db(struct):
    """Dot-bracket with only '()' holding the canonical maximum nested subset of the base pairs"""
    return pairs2db(nested_projection(db2pairs(struct)), len(struct))


def _strip_db(struct):
    """Previous behavior: keep '()' and remove any other bracket type"""
    for b in MATCHING_BRACKETS[1:]:
        struct = struct.replace(b[0], '.').replace(b[1], '.')
    return struct


def rnadistance(struct1, struct2, projection="nested"):
    """Normalized RNAdistance (full tree edit distance / length) between two dot-bracket structures.

    RNAdistance cannot handle pseudoknots (https://pubmed.ncbi.nlm.nih.gov/36077037/), so both structures are
    projected to their maximum nested subset of base pairs (projection="nested"), which does not depend on the
    bracket types used to write the structures. projection="strip" keeps the previous behavior (only '()' pairs
    are kept), where the result depends on the tool that wrote the dot-bracket.
    Uses the RNAdistance binary if available, or the ViennaRNA python package otherwise."""
    assert len(struct1)==len(struct2), "structures must have the same lenght"

    project = _nested_db if projection == "nested" else _strip_db
    struct1, struct2 = project(struct1), project(struct2)

    if shutil.which("RNAdistance"):
        oo = subprocess.check_output(["RNAdistance"], input=(struct1 + "\n" + struct2).encode('utf-8'))
        if oo.strip() == b"":
            raise ValueError("RNAdistance failed")
        dist = float(oo.decode().split(":")[1])
    else:
        import RNA
        trees = [RNA.make_tree(RNA.expand_Full(s)) for s in (struct1, struct2)]
        dist = RNA.tree_edit_distance(*trees)

    return dist / len(struct1)
