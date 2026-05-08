
from __future__ import annotations

import json
from collections import deque
from typing import Dict, List, Optional, Tuple

EPSILON: str = "__ε__"
NUCLEOTIDES = frozenset("ACGT")


class Graph:
    def __init__(self) -> None:
        self.nodes: Dict[str, Optional[str]] = {EPSILON: None}
        self.successors: Dict[str, List[str]] = {EPSILON: []}
        self.predecessors: Dict[str, List[str]] = {EPSILON: []}
        self.sink: Optional[str] = None

    def add_node(self, node_id: str, label: str) -> None:
        label = label.upper()
        if label not in NUCLEOTIDES:
            raise ValueError(f"Invalid label '{label}' for node '{node_id}'. "
                             f"Expected one of {sorted(NUCLEOTIDES)}.")
        if node_id == EPSILON:
            raise ValueError("Cannot use the reserved epsilon identifier as a node id.")
        if node_id in self.nodes:
            return
        self.nodes[node_id] = label
        self.successors[node_id] = []
        self.predecessors[node_id] = []

    def add_edge(self, u: str, v: str) -> None:
        for n in (u, v):
            if n not in self.nodes:
                raise ValueError(f"Unknown node '{n}'.")
        if v not in self.successors[u]:
            self.successors[u].append(v)
        if u not in self.predecessors[v]:
            self.predecessors[v].append(u)

    def connect_epsilon(self, node_id: str) -> None:
        self.add_edge(EPSILON, node_id)


    def topological_order(self) -> List[str]:
        """
        Return a topological ordering of all nodes via Kahn's algorithm.
        Raises ValueError if the graph contains a cycle (not a DAG).
        Epsilon is always first in the output.
        """
        in_deg: Dict[str, int] = {n: len(self.predecessors[n]) for n in self.nodes}
        queue: deque[str] = deque(n for n, d in in_deg.items() if d == 0)
        order: List[str] = []

        while queue:
            node = min(queue, key=lambda x: ("" if x == EPSILON else x))
            queue.remove(node)
            order.append(node)
            for succ in self.successors[node]:
                in_deg[succ] -= 1
                if in_deg[succ] == 0:
                    queue.append(succ)

        if len(order) != len(self.nodes):
            raise ValueError("Cycle detected – the graph is not a DAG.")
        return order


    @property
    def real_nodes(self) -> List[str]:
        return [n for n in self.nodes if n != EPSILON]

    def __repr__(self) -> str:
        edges = [(u, v) for u in self.successors for v in self.successors[u]]
        return (f"Graph(nodes={list(self.nodes.keys())}, "
                f"edges={edges}, sink={self.sink})")

    @classmethod
    def from_dict(cls, data: dict) -> "Graph":
        g = cls()

        for n in data["nodes"]:
            g.add_node(n["id"], n["label"])

        for e in data["edges"]:
            g.add_edge(e["from"], e["to"])

        sources = data["source"]
        if isinstance(sources, str):
            sources = [sources]
        for src in sources:
            g.connect_epsilon(src)

        g.sink = data["sink"]
        return g

    @classmethod
    def from_json(cls, path: str) -> "Graph":
        with open(path, encoding="utf-8") as fh:
            return cls.from_dict(json.load(fh))

    @classmethod
    def from_gfa(cls, path: str) -> "Graph":
        """
        Parse a minimal GFA v1 file (S-lines and L-lines only).

        S-line: ``S  <name>  <sequence>``
          → label = first nucleotide of *sequence* (uppercased).
        L-line: ``L  <from>  <orient>  <to>  <orient>  <overlap>``
          → directed edge from → to (orientation is ignored).

        Source nodes (no incoming edges) are automatically connected to ε.
        If multiple sink nodes are found the first one is chosen.
        """
        g = cls()
        segments: List[str] = []
        link_pairs: List[Tuple[str, str]] = []

        with open(path, encoding="utf-8") as fh:
            for raw in fh:
                line = raw.strip()
                if not line or line.startswith("#"):
                    continue
                parts = line.split("\t")
                tag = parts[0]
                if tag == "S":
                    name = parts[1]
                    seq = parts[2].upper()
                    g.add_node(name, seq[0])
                    segments.append(name)
                elif tag == "L":
                    u, v = parts[1], parts[3]
                    link_pairs.append((u, v))

        for u, v in link_pairs:
            g.add_edge(u, v)

        link_targets = {v for _, v in link_pairs}
        link_sources = {u for u, _ in link_pairs}
        sources = [n for n in segments if n not in link_targets]
        sinks = [n for n in segments if n not in link_sources]

        if not sources:
            raise ValueError("GFA: no source node found (every node has incoming edges).")
        if not sinks:
            raise ValueError("GFA: no sink node found (every node has outgoing edges).")

        for src in sources:
            g.connect_epsilon(src)
        g.sink = sinks[0]
        return g