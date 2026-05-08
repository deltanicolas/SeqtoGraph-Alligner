# Sequence-to-Graph Aligner

Global alignment of a nucleotide read **P** to a path in a **Directed Acyclic
Graph (DAG) **

**note** "This project was developed in accordance with the technical specifications provided for the Bioinformatics course at the University of Milano-Bicocca (UNIMIB). The algorithm logic and the specific scoring parameters (including the non-standard penalty matrix) strictly adhere to the academic assignment requirements."

---

## Algorithm overview

The algorithm extends **Needleman–Wunsch** to DAGs.  A fictitious epsilon node
`ε` is added as the universal source, so boundary conditions are identical to
the classic linear case.

**DP state**

```
D[i, v] = minimum cost to align prefix P[1..i] to any path ε ⇝ v
```

**Boundary conditions**

```
D[0, ε]  = 0
D[i, ε]  = i · δins                              (i ≥ 1)
D[0, v]  = min_{u ∈ Pred(v)} D[0, u] + δdel      (topological order)
```

**Recurrence** (i ≥ 1, v ≠ ε, processed in topological order)

```
D[i, v] = min {
    D[i−1, v]                         + δins,   # insertion
    min_{u∈Pred(v)} D[i−1, u] + Δ(P[i], L(v)), # match / mismatch
    min_{u∈Pred(v)} D[i,   u] + δdel,           # deletion  ← same row!
}
```

**Complexity:** O(k(n+m)) time, O(kn) space  
(O(n) space if backtracking is skipped)

---

## Default scoring matrix

|   | A | C | G | T |
|---|---|---|---|---|
| **A** | 0 | 4 | 4 | 4 |
| **C** | 4 | 1 | 4 | 4 |
| **G** | 4 | 4 | 1 | 4 |
| **T** | 4 | 4 | 4 | 1 |

Gap penalty: **δ = 4** (uniform for insertion and deletion).

---



## Installation

```bash
pip install textual          # only needed for the TUI (app.py)
```

Python 3.10+ required.  No other third-party dependencies.

---

## Usage

### CLI

```bash
# JSON graph + inline query 
python main.py -g examples/example.json -s AGT

# JSON graph + FASTA query, print full DP table
python main.py -g examples/example.json -q examples/query.fasta --table

# GFA graph, custom gap penalty
python main.py -g examples/example.gfa --gfa -s ACGTT --gap 3

# Custom scoring matrix, save JSON output
python main.py -g examples/example.json -s AGT \
    --matrix my_scores.tsv --json-out result.json

# Disable colour
python main.py -g examples/example.json -s AGT --no-color
```

**All CLI flags**

| Flag | Description |
|------|-------------|
| `-g FILE` | Graph file (required) |
| `--gfa` | Interpret graph as GFA v1 (default: JSON) |
| `-s SEQ` | Inline query sequence |
| `-q FILE` | FASTA query (first record used) |
| `--gap INT` | Uniform indel penalty δ (default: 4) |
| `--matrix FILE` | Custom 4×4 substitution matrix (TSV) |
| `--table` | Print the full DP table |
| `--json-out FILE` | Save result as JSON |
| `--no-color` | Disable ANSI colour codes |

### TUI (interactive)

```bash
python app.py
```

| Key | Action |
|-----|--------|
| `Ctrl+R` | Run alignment |
| `Ctrl+L` | Clear log |
| `Ctrl+Q` | Quit |

---

## Graph file formats

### JSON

```json
{
  "nodes": [
    {"id": "s",  "label": "A"},
    {"id": "v1", "label": "G"},
    {"id": "v3", "label": "T"},
    {"id": "v2", "label": "C"},
    {"id": "t",  "label": "T"}
  ],
  "edges": [
    {"from": "s",  "to": "v1"},
    {"from": "s",  "to": "v3"},
    {"from": "v1", "to": "v2"},
    {"from": "v2", "to": "t"},
    {"from": "v3", "to": "t"}
  ],
  "source": "s",
  "sink":   "t"
}
```

- `"source"` can be a string or a list of strings (multi-source graphs).
- All nodes listed under `"source"` are automatically connected to `ε`.

### GFA v1 (minimal subset)

Only `S` (Segment) and `L` (Link) lines are parsed.  
The label for each node is the **first nucleotide** of its sequence string.

```
H   VN:Z:1.0
S   s   ACGT
S   v1  GCGT
L   s   +   v1  +   0M
...
```

---

## Custom scoring matrix (TSV)

**Note:** This project uses a non-standard scoring matrix provided by the Bioinformatics UNIMIB course assignment. Unlike standard biological models where all matches equal 0, this matrix assigns specific costs to certain matches.

```
    A   C   G   T
A   0   4   4   4
C   4   1   4   4
G   4   4   1   4
T   4   4   4   1
```

Lines starting with `#` are ignored.

---

## Extending the aligner

- **Affine gap penalties** (Gotoh-style): add two extra DP states (`in-gap-ins`,
  `in-gap-del`).  Complexity stays O(k(n+m)).
- **Arbitrary scoring matrices**: replace `DEFAULT_SCORES` with BLOSUM62, PAM, etc.
- **Semi-global alignment**: modify boundary conditions so that the read can be
  placed anywhere in the graph (useful for short-read mapping).
- **Local alignment** (Smith–Waterman on DAG): clamp negative values to 0.
- **Multi-source graphs**: list all sources in `"source": [...]`; they are all
  connected to `ε` automatically.

---

## Example

```
Query P = AGT

Graph:  ε → s(A) → v1(G) → v2(C) → t(T)
                 ↘ v3(T) ↗

Expected DP table (Table 4):

   i\v |  ε |  s | v1 | v3 | v2 |  t
   ----+----+----+----+----+----+----
    0  |  0 |  4 |  8 |  8 | 12 | 12
    1  |  4 |  0 |  4 |  4 |  8 |  8
    2  |  8 |  4 |  1 |  4 |  5 |  8
    3  | 12 |  8 |  5 |  5 |  5 |  5

Optimal path π* = s → v3 → t  (labels: ATT)
Optimal cost = 5  (breakdown: 0 + 4 + 1)

  P :  A  G  T
  π*:  A  T  T
       |  .  |    (| match, . mismatch)

Alternative π₁ = s → v1 → v2 → t (labels: AGCT) costs 6 (Section 4.7).
```

## 📝 Author & Acknowledgments

Developed as a project for the Bioinformatics university course. The scoring matrix is based on the course materials provided for the 2026 session.