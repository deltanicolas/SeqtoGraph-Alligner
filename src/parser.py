
from __future__ import annotations

import argparse
import sys


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="seq2graph",
        description=(
            "Global Sequence-to-Graph aligner.\n"
            "Extends Needleman–Wunsch DP to directed acyclic graphs (DAGs).\n\n"
        ),
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=textwrap_dedent("""
Examples
--------
  # Align a FASTA query to a JSON graph
  python main.py -g examples/example.json -q examples/query.fasta

  # Provide the query directly on the command line
  python main.py -g examples/example.json -s AGT

  # Use a GFA graph file with custom gap penalty
  python main.py -g mygraph.gfa --gfa -s ACGT --gap 2

  # Launch the interactive TUI instead
  python app.py
        """),
    )
    graph_grp = parser.add_argument_group("Graph input")
    graph_grp.add_argument(
        "-g", "--graph",
        metavar="FILE",
        required=True,
        help="Path to the graph file (JSON or GFA).",
    )
    graph_grp.add_argument(
        "--gfa",
        action="store_true",
        default=False,
        help="Interpret the graph file as GFA (default: JSON).",
    )

    seq_grp = parser.add_argument_group("Query sequence")
    seq_excl = seq_grp.add_mutually_exclusive_group(required=True)
    seq_excl.add_argument(
        "-q", "--fasta",
        metavar="FILE",
        help="FASTA file – the first record is used as the query.",
    )
    seq_excl.add_argument(
        "-s", "--sequence",
        metavar="SEQ",
        help="Query sequence supplied directly (e.g. 'AGT').",
    )
    score_grp = parser.add_argument_group("Scoring")
    score_grp.add_argument(
        "--gap",
        metavar="INT",
        type=int,
        default=4,
        help="Uniform gap penalty δins = δdel (default: 4).",
    )
    score_grp.add_argument(
        "--matrix",
        metavar="FILE",
        default=None,
        help=(
            "Tab-separated substitution matrix file (4×4, nucleotides A/C/G/T). "
            "If omitted, the paper's default matrix is used."
        ),
    )
    out_grp = parser.add_argument_group("Output")
    out_grp.add_argument(
        "--no-color",
        action="store_true",
        default=False,
        help="Disable ANSI colour codes in the output.",
    )
    out_grp.add_argument(
        "--table",
        action="store_true",
        default=False,
        help="Print the full DP table after the alignment.",
    )
    out_grp.add_argument(
        "--json-out",
        metavar="FILE",
        default=None,
        help="Write the alignment result to a JSON file.",
    )

    return parser


def textwrap_dedent(text: str) -> str:
    import textwrap
    return textwrap.dedent(text)


def parse_args(argv=None) -> argparse.Namespace:
    return build_parser().parse_args(argv if argv is not None else sys.argv[1:])