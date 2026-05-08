
from __future__ import annotations

from dataclasses import dataclass, field
from typing import Dict, List, Optional, Tuple

from graph import Graph, EPSILON

INF = float("inf")

# ---------------------------------------------------------------------------
# Default scoring 
# ---------------------------------------------------------------------------

DEFAULT_GAP: int = 4  

DEFAULT_SCORES: Dict[Tuple[str, str], int] = {
    ("A", "A"): 0, ("A", "C"): 4, ("A", "G"): 4, ("A", "T"): 4,
    ("C", "A"): 4, ("C", "C"): 1, ("C", "G"): 4, ("C", "T"): 4,
    ("G", "A"): 4, ("G", "C"): 4, ("G", "G"): 1, ("G", "T"): 4,
    ("T", "A"): 4, ("T", "C"): 4, ("T", "G"): 4, ("T", "T"): 1,
}


# ---------------------------------------------------------------------------
# Result dataclass
# ---------------------------------------------------------------------------

@dataclass
class AlignmentResult:
    """Everything produced by a single alignment run."""
    cost: int                          # D[k, t] – optimal alignment cost
    path: List[str]                    # node ids on optimal path (excl. ε)
    path_labels: str                   # concatenated labels of path nodes
    query: str                         # original query string P
    operations: List[str]              # per-column: 'M' match, 'X' mismatch,
                                       #              'I' insertion, 'D' deletion
    aligned_query: str                 # P with gap characters inserted
    aligned_graph: str                 # graph path with gap characters inserted
    # Full DP table and topological order (for display)
    dp_table: Dict[Tuple[int, str], int] = field(repr=False)
    topo_order: List[str]              = field(repr=False)

    # ------------------------------------------------------------------

    def pretty_alignment(self, width: int = 60) -> str:
        """Return a multi-line string showing the alignment in BLAST style."""
        aq = self.aligned_query
        ag = self.aligned_graph
        mid = "".join(
            "|" if op == "M" else "." if op == "X" else " "
            for op in self.operations
        )
        lines = [f"Query  {aq}", f"       {mid}", f"Graph  {ag}"]
        return "\n".join(lines)

    def cost_breakdown(self) -> str:
        costs: List[int] = []
        for op, qc, gc in zip(self.operations, self.aligned_query, self.aligned_graph):
            if op == 'M':
                c = DEFAULT_SCORES.get((qc, gc), 0)
            elif op == 'X':
                c = DEFAULT_SCORES.get((qc, gc), 4)
            else:
                c = DEFAULT_GAP
            costs.append(c)
        parts = " + ".join(str(c) for c in costs)
        return f"{parts} = {self.cost}"


# ---------------------------------------------------------------------------
# Core alignment function
# ---------------------------------------------------------------------------

def align(
    query: str,
    graph: Graph,
    scores: Optional[Dict[Tuple[str, str], int]] = None,
    gap: int = DEFAULT_GAP,
) -> AlignmentResult:
    """
    Perform global sequence-to-graph alignment.

    Parameters
    ----------
    query  : nucleotide string P (A/C/G/T characters only).
    graph  : a Graph instance with epsilon node and defined sink.
    scores : substitution cost matrix; defaults to the paper's matrix.
    gap    : uniform indel cost δins = δdel.

    Returns
    -------
    AlignmentResult with cost, optimal path, aligned sequences and full DP table.
    """
    if scores is None:
        scores = DEFAULT_SCORES

    if graph.sink is None:
        raise ValueError("Graph has no sink node defined.")

    query = query.upper()
    k = len(query)
    topo = graph.topological_order()  
   
    D: Dict[Tuple[int, str], float] = {}
    B: Dict[Tuple[int, str], Optional[Tuple[int, str]]] = {}

    # ----------------------------------------------------------------
    # Boundary conditions
    # ----------------------------------------------------------------

    
    D[(0, EPSILON)] = 0
    B[(0, EPSILON)] = None


    for i in range(1, k + 1):
        D[(i, EPSILON)] = i * gap
        B[(i, EPSILON)] = (i - 1, EPSILON)


    for v in topo:
        if v == EPSILON:
            continue
        best: float = INF
        best_pred: Optional[str] = None
        for u in graph.predecessors[v]:
            cand = D.get((0, u), INF)
            if cand < INF:
                cand += gap
            if cand < best:
                best = cand
                best_pred = u
        D[(0, v)] = best
        B[(0, v)] = (0, best_pred) if best_pred is not None else None

    # ----------------------------------------------------------------
    # Fill 
    # ----------------------------------------------------------------

    def delta(c1: str, c2: str) -> int:
        return scores.get((c1, c2), 4)

    for i in range(1, k + 1):
        pi = query[i - 1]         
        for v in topo:
            if v == EPSILON:
                continue
            lv = graph.nodes[v]   


            best = D.get((i - 1, v), INF)
            if best < INF:
                best += gap
            best_from: Optional[Tuple[int, str]] = (i - 1, v) if best < INF else None

            for u in graph.predecessors[v]:
                prev = D.get((i - 1, u), INF)
                if prev < INF:
                    cand = prev + delta(pi, lv)  # type: ignore
                    if cand < best:
                        best = cand
                        best_from = (i - 1, u)


            for u in graph.predecessors[v]:
                prev = D.get((i, u), INF)
                if prev < INF:
                    cand = prev + gap
                    if cand < best:
                        best = cand
                        best_from = (i, u)

            D[(i, v)] = best
            B[(i, v)] = best_from

    # ----------------------------------------------------------------
    # Optimal cost
    # ----------------------------------------------------------------

    opt = D.get((k, graph.sink), INF)
    if opt == INF:
        raise ValueError(
            f"Sink '{graph.sink}' is unreachable from ε — check graph connectivity."
        )

    # ----------------------------------------------------------------
    # Backtracking
    # ----------------------------------------------------------------

    path_nodes: List[str] = []
    ops: List[str] = []
    aq: List[str] = []   
    ag: List[str] = [] 

    cur_i, cur_v = k, graph.sink

    while not (cur_i == 0 and cur_v == EPSILON):
        back = B.get((cur_i, cur_v))
        if back is None:
            # Reached a boundary cell that points to None (e.g. (0, ε))
            break
        prev_i, prev_v = back

        if prev_i == cur_i - 1 and prev_v == cur_v:
            ops.append("I")
            aq.append(query[cur_i - 1])
            ag.append("-")
            cur_i = prev_i
 

        elif prev_i == cur_i - 1:
            if cur_v != EPSILON:
                path_nodes.append(cur_v)
                lv = graph.nodes[cur_v]
                ops.append("M" if query[cur_i - 1] == lv else "X")
                aq.append(query[cur_i - 1])
                ag.append(lv)
            cur_i = prev_i
            cur_v = prev_v

        else:
            if cur_v != EPSILON:
                path_nodes.append(cur_v)
                lv = graph.nodes[cur_v]
                ops.append("D")
                aq.append("-")
                ag.append(lv)
            cur_v = prev_v

    path_nodes.reverse()
    ops.reverse()
    aq.reverse()
    ag.reverse()

    path_nodes = [n for n in path_nodes if n != EPSILON]
    path_labels = "".join(graph.nodes[n] for n in path_nodes) # type: ignore

    dp_int: Dict[Tuple[int, str], int] = {
        k: int(v) for k, v in D.items() if v < INF
    }

    return AlignmentResult(
        cost=int(opt),
        path=path_nodes,
        path_labels=path_labels,
        query=query,
        operations=ops,
        aligned_query="".join(aq),
        aligned_graph="".join(ag),
        dp_table=dp_int,
        topo_order=topo,
    )