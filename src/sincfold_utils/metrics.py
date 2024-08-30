from sklearn.metrics import f1_score
import numpy as np

def f1_triangular(ref, pred):
    """Compute F1 from the upper triangular connection matrix. ref and pred are binary 2D numpy arrays"""
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

    return f1_score(ref, pred, zero_division=0)

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
