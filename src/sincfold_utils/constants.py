# Mapping of nucleotide symbols
# R	Guanine / Adenine (purine)
# Y	Cytosine / Uracil (pyrimidine)
# K	Guanine / Uracil
# M	Adenine / Cytosine
# S	Guanine / Cytosine
# W	Adenine / Uracil
# B	Guanine / Uracil / Cytosine
# D	Guanine / Adenine / Uracil
# H	Adenine / Cytosine / Uracil
# V	Guanine / Cytosine / Adenine
# N	Adenine / Guanine / Cytosine / Uracil
NT_DICT = {
    "R": ["G", "A"],
    "Y": ["C", "U"],
    "K": ["G", "U"],
    "M": ["A", "C"],
    "S": ["G", "C"],
    "W": ["A", "U"],
    "B": ["G", "U", "C"],
    "D": ["G", "A", "U"],
    "H": ["A", "C", "U"],
    "V": ["G", "C", "A"],
    "N": ["G", "A", "C", "U"],
}

VOCABULARY = ["A", "C", "G", "U"]

# Bracket types for pseudoknot pages, in the standard order (ViennaRNA, bpRNA, RNApolis):
# ()[]{}<> followed by Aa, Bb, ..., Zz
STANDARD_BRACKETS = [["(", ")"], ["[", "]"], ["{", "}"], ["<", ">"]] + [
    [chr(c), chr(c).lower()] for c in range(ord("A"), ord("Z") + 1)
]
# Order used by RNAstructure ct2dot: ( < { [ A B ...
RNASTRUCTURE_BRACKETS = [["(", ")"], ["<", ">"], ["{", "}"], ["[", "]"]] + STANDARD_BRACKETS[4:]

MATCHING_BRACKETS = STANDARD_BRACKETS
OPEN_TO_CLOSE = {op: cl for op, cl in STANDARD_BRACKETS}
CLOSE_TO_OPEN = {cl: op for op, cl in STANDARD_BRACKETS}