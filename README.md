# sincfold_utils
Simple package to work with RNA sequences and structures 


# TODO
- hide output messages from predictors if theres no errors.
- raise error if predictor not installed 

# Dot-bracket and pseudoknots
`sincfold_utils.structure` converts between base pairs (1-based `[[i, j], ...]`, i<j) and dot-bracket
in pure python, including pseudoknots:

```python
from sincfold_utils.structure import pairs2db, db2pairs, nested_projection, normalize_db
pairs2db([[1, 16], [2, 15], [9, 24], [10, 23]], 24)  # '((......[[....))......]]'
db2pairs("((<<))>>")          # RNAstructure (ct2dot) brackets are also parsed
normalize_db("((<<))>>")      # '(([[))]]'
```
Page 0 (`()`) is the maximum nested subset of pairs, page 1 (`[]`) the maximum nested subset of the
rest, and so on (`()[]{}<>AaBb...Zz`).
`io.bp2dot`, `io.mat2dot`, `io.dot2bp` and `io.ct2dot` use it (`ct2dot(f, brackets=RNASTRUCTURE_BRACKETS)` gives the
same string as RNAstructure ct2dot). The bundled binary is still available as `io.ct2dot_rnastructure` for comparisons.
