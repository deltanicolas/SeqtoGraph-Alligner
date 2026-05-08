
from __future__ import annotations

import re
from pathlib import Path
from typing import Dict, Generator, List, Optional, Tuple

from graph import EPSILON


# ---------------------------------------------------------------------------
# FASTA I/O
# ---------------------------------------------------------------------------

def parse_fasta(path: str) -> Generator[Tuple[str, str], None, None]:
    name: Optional[str] = None
    chunks: List[str] = []

    def _emit(n: str, seqs: List[str]) -> Tuple[str, str]:
        seq = "".join(seqs).upper().replace(" ", "").replace("\n", "")
        bad = re.sub(r"[ACGT]", "", seq)
        if bad:
            raise ValueError(
                f"Sequence '{n}' contains invalid characters: {set(bad)}"
            )
        return n, seq

    with open(path, encoding="utf-8") as fh:
        for line in fh:
            line = line.rstrip("\n")
            if line.startswith(">"):
                if name is not None:
                    yield _emit(name, chunks)
                name = line[1:].strip()
                chunks = []
            elif line and name is not None:
                chunks.append(line)

    if name is not None:
        yield _emit(name, chunks)


def first_sequence(path: str) -> str:
    for _, seq in parse_fasta(path):
        return seq
    raise ValueError(f"No sequences found in '{path}'.")


# ---------------------------------------------------------------------------
# Scoring matrix I/O
# ---------------------------------------------------------------------------

NUCLEOTIDES = ("A", "C", "G", "T")


def default_score_matrix() -> Dict[Tuple[str, str], int]:
    from aligner import DEFAULT_SCORES
    return dict(DEFAULT_SCORES)


def load_score_matrix(path: str) -> Dict[Tuple[str, str], int]:
    """
    Load a custom scoring matrix from a tab-separated file.

    Format (header row + 4 data rows)::

          A   C   G   T
        A 0   4   4   4
        C 4   1   4   4
        G 4   4   1   4
        T 4   4   4   1

    Returns a dict keyed by (row_nucleotide, col_nucleotide).
    """
    matrix: Dict[Tuple[str, str], int] = {}
    with open(path, encoding="utf-8") as fh:
        lines = [l.strip() for l in fh if l.strip() and not l.startswith("#")]
    if not lines:
        raise ValueError("Score matrix file is empty.")

    cols = lines[0].split()
    for row_line in lines[1:]:
        parts = row_line.split()
        row_nuc = parts[0].upper()
        for col_nuc, val in zip(cols, parts[1:]):
            matrix[(row_nuc, col_nuc.upper())] = int(val)
    for a in NUCLEOTIDES:
        for b in NUCLEOTIDES:
            if (a, b) not in matrix:
                raise ValueError(f"Score matrix missing entry ({a}, {b}).")
    return matrix


# ---------------------------------------------------------------------------
# Pretty-printers
# ---------------------------------------------------------------------------

RESET = "\033[0m"
BOLD  = "\033[1m"
GREEN = "\033[92m"
RED   = "\033[91m"
YELLOW = "\033[93m"
CYAN   = "\033[96m"
DIM    = "\033[2m"


def _colored(text: str, color: str, use_color: bool = True) -> str:
    return f"{color}{text}{RESET}" if use_color else text


def print_alignment(result, use_color: bool = True) -> None:
    """Print a full alignment report to stdout."""
    from aligner import AlignmentResult
    r: AlignmentResult = result

    sep = "─" * 60
    print(f"\n{_colored(sep, BOLD, use_color)}")
    print(f"{_colored('  Sequence-to-Graph Alignment Result', BOLD, use_color)}")
    print(f"{_colored(sep, BOLD, use_color)}")

    arrow = " → "
    path_str = arrow.join(r.path) if r.path else "(empty)"
    print(f"\n  Optimal path : {_colored(path_str, CYAN, use_color)}")
    print(f"  Path labels  : {_colored(r.path_labels, CYAN, use_color)}")
    print(f"  Optimal cost : {_colored(str(r.cost), GREEN, use_color)}")

    print(f"\n  {_colored('Alignment', BOLD, use_color)}\n")
    aq = r.aligned_query
    ag = r.aligned_graph
    mid = ""
    for op in r.operations:
        if op == "M":
            mid += _colored("|", GREEN, use_color)
        elif op == "X":
            mid += _colored(".", YELLOW, use_color)
        else:
            mid += " "

    print(f"  Query  {aq}")
    print(f"         {mid}")
    print(f"  Graph  {ag}")

    from aligner import DEFAULT_SCORES, DEFAULT_GAP
    costs: List[int] = []
    for op, qc, gc in zip(r.operations, aq, ag):
        if op in ("M", "X"):
            costs.append(DEFAULT_SCORES.get((qc, gc), 4))
        else:
            costs.append(DEFAULT_GAP)

    cost_str = " + ".join(str(c) for c in costs) + f" = {r.cost}"
    print(f"\n  Cost breakdown: {cost_str}")

    legend = "  Legend: " + "  ".join([
        _colored("| match", GREEN, use_color),
        _colored(". mismatch", YELLOW, use_color),
        "  insertion/deletion",
    ])
    print(f"\n{legend}")
    print(f"\n{_colored(sep, BOLD, use_color)}\n")


def format_dp_table(result, max_cols: int = 40) -> str:
    from aligner import AlignmentResult
    r: AlignmentResult = result

    topo = r.topo_order
    k = len(r.query)
    col_labels = ["ε" if v == EPSILON else v for v in topo]
    col_w = max(len(c) for c in col_labels) + 2
    row_label_w = max(len(str(k)) + 4, 6)

    if (len(topo) + 1) * col_w > max_cols * 2:
        return "(DP table too wide to display — use --no-table to suppress)"

    def _cell(i: int, v: str) -> str:
        val = r.dp_table.get((i, v))
        s = "∞" if val is None else str(val)
        return s.center(col_w)

    header_row = " " * row_label_w + "".join(lbl.center(col_w) for lbl in col_labels)
    sep_row = "-" * len(header_row)

    rows = [sep_row, header_row, sep_row]
    for i in range(k + 1):
        label = f"i={i}".ljust(row_label_w)
        cells = "".join(_cell(i, v) for v in topo)
        rows.append(label + cells)
    rows.append(sep_row)

    return "\n".join(rows)