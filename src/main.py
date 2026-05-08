
from __future__ import annotations

import json
import sys


def main() -> int:
    from parser import parse_args
    from graph import Graph
    from aligner import align, DEFAULT_GAP, DEFAULT_SCORES
    from utils import (
        first_sequence, load_score_matrix,
        print_alignment, format_dp_table,
    )

    args = parse_args()
    try:
        if args.gfa:
            graph = Graph.from_gfa(args.graph)
        else:
            graph = Graph.from_json(args.graph)
    except (FileNotFoundError, ValueError, KeyError) as exc:
        print(f"[ERROR] Failed to load graph: {exc}", file=sys.stderr)
        return 1
    try:
        if args.fasta:
            query = first_sequence(args.fasta)
        else:
            query = args.sequence.upper().strip()
    except (FileNotFoundError, ValueError) as exc:
        print(f"[ERROR] Failed to load query: {exc}", file=sys.stderr)
        return 1

    if not query:
        print("[ERROR] Query sequence is empty.", file=sys.stderr)
        return 1
    try:
        scores = load_score_matrix(args.matrix) if args.matrix else None
    except (FileNotFoundError, ValueError) as exc:
        print(f"[ERROR] Failed to load scoring matrix: {exc}", file=sys.stderr)
        return 1
    use_color = not args.no_color
    try:
        result = align(query, graph, scores=scores, gap=args.gap)
    except ValueError as exc:
        print(f"[ERROR] Alignment failed: {exc}", file=sys.stderr)
        return 1
    print_alignment(result, use_color=use_color)

    if args.table:
        print("  DP Table\n")
        print(format_dp_table(result))
        print()

    if args.json_out:
        payload = {
            "query": result.query,
            "cost": result.cost,
            "path": result.path,
            "path_labels": result.path_labels,
            "aligned_query": result.aligned_query,
            "aligned_graph": result.aligned_graph,
            "operations": result.operations,
        }
        with open(args.json_out, "w", encoding="utf-8") as fh:
            json.dump(payload, fh, indent=2)
        print(f"  Result saved to '{args.json_out}'.")

    return 0


if __name__ == "__main__":
    sys.exit(main())