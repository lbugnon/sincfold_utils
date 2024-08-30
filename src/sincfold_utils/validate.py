from sincfold_utils.constants import NT_DICT, VOCABULARY

def is_valid_sequence(seq):
    """Check if sequence is valid"""
    return set(seq.upper()) <= (set(NT_DICT.keys()).union(set(VOCABULARY)))