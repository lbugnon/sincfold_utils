"""Pure python/numpy conversions between base pairs and dot-bracket notation, with pseudoknot support.

Base pairs follow the package convention: a sorted list of ``[i, j]`` lists, 1-based, with ``i < j``.

Pseudoknots are encoded with nested "pages": page 0 is the maximum set of non-crossing base pairs (written
with ``()``), page 1 is the maximum non-crossing set among the remaining pairs (``[]``), and so on. The page
assignment depends only on the set of pairs, so the same structure always gets the same dot-bracket string,
no matter which tool produced it.
"""
import numpy as np

from sincfold_utils.constants import STANDARD_BRACKETS, OPEN_TO_CLOSE, CLOSE_TO_OPEN

UNPAIRED = "."


def validate_pairs(pairs, L=None):
    """Return pairs as a sorted list of [i, j] (1-based, i < j).

    Pairs given as (j, i) are flipped and exact duplicates are removed. Raises ValueError on
    self pairs, positions out of range, or nucleotides with more than one partner."""
    out = set()
    for bp in pairs:
        if len(bp) != 2:
            raise ValueError(f"base pair {bp} must have two positions")
        i, j = int(bp[0]), int(bp[1])
        if i == j:
            raise ValueError(f"self pair {bp}")
        if i > j:
            i, j = j, i
        if i < 1 or (L is not None and j > L):
            raise ValueError(f"base pair {bp} out of range 1..{L}")
        out.add((i, j))

    partner = {}
    for i, j in out:
        for x, y in ((i, j), (j, i)):
            if x in partner:
                raise ValueError(f"nucleotide {x} paired with {partner[x]} and {y} (multiplet)")
            partner[x] = y
    return [[i, j] for i, j in sorted(out)]


def pairs2ptable(pairs, L):
    """Partner array of length L: ptable[k-1] is the 1-based partner of nucleotide k, 0 if unpaired."""
    ptable = np.zeros(L, dtype=int)
    for i, j in validate_pairs(pairs, L):
        ptable[i - 1] = j
        ptable[j - 1] = i
    return ptable


def _stems(pairs):
    """Group sorted, valid pairs into stems of stacked pairs (i, j), (i+1, j-1), ...
    Two stems either cross completely or not at all, so pages can be computed on stems."""
    partner = {i: j for i, j in pairs}
    stems, seen = [], set()
    for i, j in pairs:
        if i in seen:
            continue
        stem = []
        while partner.get(i) == j:
            stem.append([i, j])
            seen.add(i)
            i, j = i + 1, j - 1
        stems.append(stem)
    return stems


def _max_nested_stems(stems):
    """Indices of the maximum-weight (number of pairs) subset of non-crossing stems.

    Interval DP over the compressed stem endpoints. On ties, the stem starting further 5' is kept."""
    if len(stems) == 0:
        return []
    ends = []  # (position, stem index)
    for k, s in enumerate(stems):
        ends.append((s[0][0], k))
        ends.append((s[0][1], k))
    ends.sort()
    m = len(ends)
    stem_of = [k for _, k in ends]
    other = [0] * m  # index of the other endpoint of the same stem
    first = {}
    for p, k in enumerate(stem_of):
        if k in first:
            other[first[k]] = p
            other[p] = first[k]
        else:
            first[k] = p
    weight = [len(s) for s in stems]

    # best[a][b]: best weight using stems fully inside compressed positions a..b
    best = [[0] * (m + 1) for _ in range(m + 2)]
    take = [[False] * (m + 1) for _ in range(m + 2)]
    for a in range(m - 1, -1, -1):
        row, row_next = best[a], best[a + 1]
        c = other[a]
        w = weight[stem_of[a]]
        for b in range(a, m):
            skip = row_next[b]
            if a < c <= b:
                inc = w + best[a + 1][c - 1] + (best[c + 1][b] if c + 1 <= b else 0)
                if inc >= skip:
                    row[b] = inc
                    take[a][b] = True
                    continue
            row[b] = skip

    # traceback
    chosen, stack = [], [(0, m - 1)]
    while stack:
        a, b = stack.pop()
        while a <= b:
            if take[a][b]:
                c = other[a]
                chosen.append(stem_of[a])
                if c + 1 <= b:
                    stack.append((c + 1, b))
                a, b = a + 1, c - 1
            else:
                a += 1
    return chosen


def pairs2pages(pairs):
    """Split base pairs into pages (lists of non-crossing pairs). Page 0 is the maximum nested subset."""
    pairs = validate_pairs(pairs)
    pages = []
    while pairs:
        stems = _stems(pairs)
        chosen = set(_max_nested_stems(stems))
        page = sorted(bp for k in chosen for bp in stems[k])
        pages.append(page)
        page_set = set(map(tuple, page))
        pairs = [bp for bp in pairs if tuple(bp) not in page_set]
    return pages


def pk_order(pairs):
    """Number of bracket types needed to write the structure (0 if no pairs, 1 if nested)."""
    return len(pairs2pages(pairs))


def nested_projection(pairs):
    """Maximum nested (pseudoknot-free) subset of the base pairs."""
    pages = pairs2pages(pairs)
    return pages[0] if pages else []


def crossing_pairs(pairs):
    """Base pairs that cross at least one other base pair."""
    pairs = validate_pairs(pairs)
    if len(pairs) == 0:
        return []
    bp = np.array(pairs)
    i, j = bp[:, 0][:, None], bp[:, 1][:, None]
    cross = ((i < i.T) & (i.T < j) & (j < j.T)) | ((i.T < i) & (i < j.T) & (j.T < j))
    return [pairs[k] for k in np.nonzero(cross.any(axis=1))[0]]


def pairs2db(pairs, L, brackets=STANDARD_BRACKETS):
    """Dot-bracket string of length L from 1-based base pairs. Pseudoknots use the next bracket types
    in ``brackets`` (default ()[]{}<>AaBb...Zz)."""
    pairs = validate_pairs(pairs, L)
    pages = pairs2pages(pairs)
    if len(pages) > len(brackets):
        raise ValueError(f"structure needs {len(pages)} bracket types, only {len(brackets)} available")
    db = [UNPAIRED] * L
    for (op, cl), page in zip(brackets, pages):
        for i, j in page:
            db[i - 1] = op
            db[j - 1] = cl
    return "".join(db)


def db2pairs(db, unpaired=UNPAIRED):
    """Sorted 1-based base pairs from a dot-bracket string. Accepts any bracket type in
    ()[]{}<>AaBb...Zz, so both the standard and the RNAstructure (ct2dot) orders are parsed.
    Raises ValueError on unbalanced brackets or unknown characters."""
    stacks = {op: [] for op in OPEN_TO_CLOSE}
    pairs = []
    for k, c in enumerate(db, start=1):
        if c in unpaired:
            continue
        if c in stacks:
            stacks[c].append(k)
        elif c in CLOSE_TO_OPEN:
            op = CLOSE_TO_OPEN[c]
            if not stacks[op]:
                raise ValueError(f"unbalanced '{c}' at position {k}")
            pairs.append([stacks[op].pop(), k])
        else:
            raise ValueError(f"unknown character '{c}' at position {k}")
    for op, stack in stacks.items():
        if stack:
            raise ValueError(f"unbalanced '{op}' at position {stack[-1]}")
    return sorted(pairs)


def normalize_db(db, brackets=STANDARD_BRACKETS):
    """Re-encode a dot-bracket string with canonical pages, e.g. RNAstructure '<>' becomes '[]'."""
    return pairs2db(db2pairs(db), len(db), brackets)
