#!/usr/bin/env python3
"""
app.py – Interactive TUI for the Sequence-to-Graph aligner.

Built with the Textual framework.  Layout:

  ┌─ Header ───────────────────────────────────────────────┐
  │  Sequence-to-Graph Aligner                             │
  ├─ Sidebar ──┬─ Main panel ───────────────────────────── │
  │  Graph     │  [Tabs: DP Table | Alignment | About]      │
  │  Sequence  │                                            │
  │  Gap / Mtx │                                            │
  │  [Align]   │                                            │
  ├────────────┴───────────────────────────────────────────┤
  │  Log                                                    │
  └────────────────────────────────────────────────────────┘

Run:
    python app.py
"""
from __future__ import annotations

import sys
import traceback
from pathlib import Path
from typing import Optional

# ── Textual imports ────────────────────────────────────────────────────────
try:
    from textual.app import App, ComposeResult
    from textual.binding import Binding
    from textual.containers import Container, Horizontal, Vertical, ScrollableContainer
    from textual.widgets import (
        Button, DataTable, Footer, Header, Input, Label,
        Log, Markdown, Static, TabbedContent, TabPane,
    )
    from textual.reactive import reactive
except ImportError:
    print(
        "[ERROR] Textual is not installed.\n"
        "  pip install textual\n"
        "Alternatively, use the CLI:\n"
        "  python main.py -g examples/example.json -s AGT --table",
        file=sys.stderr,
    )
    sys.exit(1)

# ── Project imports ────────────────────────────────────────────────────────
sys.path.insert(0, str(Path(__file__).parent))
from graph import Graph, EPSILON
from aligner import align, AlignmentResult, DEFAULT_GAP
from utils import first_sequence, format_dp_table


# ──────────────────────────────────────────────────────────────────────────
# CSS
# ──────────────────────────────────────────────────────────────────────────

APP_CSS = """
Screen {
    background: $surface;
}

Header {
    background: $accent;
    color: $text;
    text-style: bold;
}

/* ── Sidebar ── */
#sidebar {
    width: 28;
    min-width: 24;
    background: $panel;
    border-right: solid $primary;
    padding: 1 1;
}

#sidebar Label {
    color: $text-muted;
    text-style: bold;
    margin-top: 1;
}

#sidebar Input {
    margin-bottom: 1;
}

#btn-align {
    width: 100%;
    margin-top: 1;
    background: $success;
    color: $text;
    text-style: bold;
}

#btn-align:hover {
    background: $success-darken-1;
}

#btn-clear {
    width: 100%;
    margin-top: 0;
    background: $error;
    color: $text;
}

/* ── Main content ── */
#main-content {
    background: $surface;
    padding: 0 1;
}

/* ── Log pane ── */
#log-pane {
    height: 7;
    border-top: solid $primary;
    background: $panel;
}

#log-pane Label {
    background: $primary;
    color: $text;
    text-style: bold;
    width: 100%;
    padding: 0 1;
}

Log {
    height: 5;
    background: $panel;
}

/* ── DP Table tab ── */
#dp-table-container {
    padding: 1;
}

DataTable {
    height: auto;
    max-height: 30;
}

/* ── Alignment tab ── */
#alignment-display {
    padding: 1 2;
}

.result-section {
    margin-bottom: 1;
}

.result-label {
    text-style: bold;
    color: $accent;
}

.result-value {
    color: $text;
}

.align-match    { color: $success; }
.align-mismatch { color: $warning; }
.align-gap      { color: $text-muted; }

/* ── About tab ── */
#about-pane {
    padding: 2;
}
"""

# ──────────────────────────────────────────────────────────────────────────
# About page content (Markdown)
# ──────────────────────────────────────────────────────────────────────────

ABOUT_MD = """
# Sequence-to-Graph Aligner

Global alignment of a nucleotide sequence **P** to a path in a
**Directed Acyclic Graph (DAG)**, extending the classic Needleman–Wunsch
algorithm.

## Algorithm

The DP state **D[i, v]** stores the minimum cost to align prefix P[1..i]
to any path ε ⇝ v.

| Boundary / Recurrence | Formula |
|---|---|
| Origin | D[0, ε] = 0 |
| Column ε | D[i, ε] = i · δ_ins |
| Row 0 (deletions) | D[0, v] = min_u D[0, u] + δ_del |
| General (i ≥ 1, v ≠ ε) | min{ insertion, match/mismatch, deletion } |

## Default Scoring Matrix

| | A | C | G | T |
|---|---|---|---|---|
| **A** | 0 | 4 | 4 | 4 |
| **C** | 4 | 1 | 4 | 4 |
| **G** | 4 | 4 | 1 | 4 |
| **T** | 4 | 4 | 4 | 1 |

Gap penalty: **δ = 4** (insertion = deletion).

## File formats

- **Graph**: JSON (`-g file.json`) or GFA v1 (`-g file.gfa --gfa`)
- **Query**: FASTA file or direct string on the CLI (`-s AGT`)

## Example JSON graph

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
  "sink": "t"
}
```
"""


# ──────────────────────────────────────────────────────────────────────────
# Helper: render aligned sequences as Rich markup
# ──────────────────────────────────────────────────────────────────────────

def _render_alignment_md(r: AlignmentResult) -> str:
    aq = r.aligned_query
    ag = r.aligned_graph
    ops = r.operations

    mid = ""
    for op in ops:
        if op == "M":
            mid += "|"
        elif op == "X":
            mid += "."
        else:
            mid += " "

    from aligner import DEFAULT_SCORES, DEFAULT_GAP
    costs = []
    for op, qc, gc in zip(ops, aq, ag):
        if op in ("M", "X"):
            costs.append(str(DEFAULT_SCORES.get((qc, gc), 4)))
        else:
            costs.append(str(DEFAULT_GAP))

    arrow = " → ".join(r.path) if r.path else "(empty)"

    lines = [
        f"## Optimal path",
        f"**{arrow}**  (labels: `{r.path_labels}`)",
        f"",
        f"## Optimal cost:  **{r.cost}**",
        f"",
        f"## Alignment",
        f"```",
        f"Query  {aq}",
        f"       {mid}",
        f"Graph  {ag}",
        f"```",
        f"",
        f"**Cost breakdown:** `{' + '.join(costs)} = {r.cost}`",
        f"",
        f"*Legend:  `|` match  `.` mismatch  ` ` insertion/deletion*",
    ]
    return "\n".join(lines)


# ──────────────────────────────────────────────────────────────────────────
# Main App
# ──────────────────────────────────────────────────────────────────────────

class Seq2GraphApp(App):
    """Interactive TUI for the Sequence-to-Graph aligner."""

    TITLE = "Sequence-to-Graph Aligner"
    SUB_TITLE = "Chines (2026) · Needleman–Wunsch on DAGs"
    CSS = APP_CSS

    BINDINGS = [
        Binding("ctrl+q", "quit", "Quit"),
        Binding("ctrl+r", "run_alignment", "Run"),
        Binding("ctrl+l", "clear_log", "Clear log"),
    ]
    _result: Optional[AlignmentResult] = None

    # ── Layout ─────────────────────────────────────────────────────────────

    def compose(self) -> ComposeResult:
        yield Header()

        with Horizontal():
            # ── Sidebar ──────────────────────────────────────────────────
            with Vertical(id="sidebar"):
                yield Label("Graph file (JSON/GFA)")
                yield Input(
                    placeholder="examples/example.json",
                    id="input-graph",
                )
                yield Label("Format")
                yield Input(
                    placeholder="json  (or gfa)",
                    id="input-format",
                    value="json",
                )
                yield Label("Query sequence / FASTA")
                yield Input(
                    placeholder="AGT  or  path/to/query.fasta",
                    id="input-query",
                )
                yield Label("Gap penalty (δ)")
                yield Input(
                    placeholder="4",
                    id="input-gap",
                    value="4",
                )
                yield Label("Scoring matrix (optional)")
                yield Input(
                    placeholder="path/to/matrix.tsv",
                    id="input-matrix",
                )
                yield Button("▶  Align  (Ctrl+R)", id="btn-align", variant="success")
                yield Button("✕  Clear log", id="btn-clear", variant="error")

            # ── Main content ──────────────────────────────────────────────
            with Vertical(id="main-content"):
                with TabbedContent(id="tabs"):
                    with TabPane("DP Table", id="tab-dp"):
                        with ScrollableContainer(id="dp-table-container"):
                            yield DataTable(id="dp-table", zebra_stripes=True)

                    with TabPane("Alignment", id="tab-align"):
                        with ScrollableContainer(id="alignment-display"):
                            yield Markdown(
                                "Run an alignment to see results here.",
                                id="alignment-md",
                            )

                    with TabPane("About", id="tab-about"):
                        with ScrollableContainer(id="about-pane"):
                            yield Markdown(ABOUT_MD, id="about-md")

        # ── Log ──────────────────────────────────────────────────────────
        with Vertical(id="log-pane"):
            yield Label("  ■ Log")
            yield Log(id="log", auto_scroll=True)

        yield Footer()


    def on_mount(self) -> None:
        self._log("TUI ready. Fill in the fields and press ▶ Align.")
        self._log("Tip: use Ctrl+R to align, Ctrl+Q to quit.")
        example = Path(__file__).parent / "examples" / "example.json"
        if example.exists():
            self.query_one("#input-graph", Input).value = str(example)
            self.query_one("#input-query", Input).value = "AGT"
            self._log(f"Pre-loaded example: {example}")

    def on_button_pressed(self, event: Button.Pressed) -> None:
        if event.button.id == "btn-align":
            self.action_run_alignment()
        elif event.button.id == "btn-clear":
            self.action_clear_log()

    def action_run_alignment(self) -> None:
        graph_path = self.query_one("#input-graph", Input).value.strip()
        fmt        = self.query_one("#input-format", Input).value.strip().lower()
        raw_query  = self.query_one("#input-query", Input).value.strip()
        raw_gap    = self.query_one("#input-gap", Input).value.strip()
        mtx_path   = self.query_one("#input-matrix", Input).value.strip()

        if not graph_path:
            self._log("[ERROR] Graph file path is required.", error=True)
            return
        if not raw_query:
            self._log("[ERROR] Query sequence or FASTA path is required.", error=True)
            return
        try:
            gap = int(raw_gap) if raw_gap else DEFAULT_GAP
        except ValueError:
            self._log(f"[ERROR] Invalid gap penalty '{raw_gap}' – must be an integer.", error=True)
            return

        # Load graph
        try:
            if fmt == "gfa":
                graph = Graph.from_gfa(graph_path)
                self._log(f"Loaded GFA graph: {graph_path}")
            else:
                graph = Graph.from_json(graph_path)
                self._log(f"Loaded JSON graph: {graph_path}")
        except Exception as exc:
            self._log(f"[ERROR] Could not load graph: {exc}", error=True)
            return

        self._log(
            f"  Nodes: {len(graph.real_nodes)}  "
            f"Sink: {graph.sink}  "
            f"Topo order: {graph.topological_order()}"
        )

        # Load query
        try:
            if raw_query.upper().endswith(".FASTA") or raw_query.upper().endswith(".FA"):
                query = first_sequence(raw_query)
                self._log(f"Loaded FASTA query: {raw_query}  seq={query}")
            else:
                query = raw_query.upper()
                self._log(f"Using inline query: {query}")
        except Exception as exc:
            self._log(f"[ERROR] Could not load query: {exc}", error=True)
            return

        # Load scoring matrix
        scores = None
        if mtx_path:
            try:
                from utils import load_score_matrix
                scores = load_score_matrix(mtx_path)
                self._log(f"Loaded scoring matrix: {mtx_path}")
            except Exception as exc:
                self._log(f"[ERROR] Could not load matrix: {exc}", error=True)
                return

        # Run alignment
        try:
            self._log(f"Running alignment  P={query}  gap={gap} …")
            result = align(query, graph, scores=scores, gap=gap)
            self._result = result
            self._log(
                f"✓ Done!  cost={result.cost}  "
                f"path={' → '.join(result.path)}"
            )
        except Exception as exc:
            self._log(f"[ERROR] Alignment failed: {exc}", error=True)
            self._log(traceback.format_exc())
            return

        self._populate_dp_table(result)
        self._populate_alignment_tab(result)

    def action_clear_log(self) -> None:
        self.query_one("#log", Log).clear()

    # ── Private helpers ────────────────────────────────────────────────────────── 

    def _log(self, message: str, error: bool = False) -> None:
        log_widget = self.query_one("#log", Log)
        prefix = "✗ " if error else "  "
        log_widget.write_line(prefix + message)

    def _populate_dp_table(self, result: AlignmentResult) -> None:
        table = self.query_one("#dp-table", DataTable)
        table.clear(columns=True)

        topo = result.topo_order
        k = len(result.query)

        col_labels = ["i \\ v"] + ["ε" if v == EPSILON else v for v in topo]
        for lbl in col_labels:
            table.add_column(lbl, key=lbl)

        for i in range(k + 1):
            label = f"i={i}"
            row_data = [label]
            for v in topo:
                val = result.dp_table.get((i, v))
                cell = "∞" if val is None else str(val)
                if v in result.path or v == EPSILON:
                    row_data.append(f"[bold]{cell}[/bold]")
                else:
                    row_data.append(cell)
            table.add_row(*row_data)

        self._log(f"DP table populated  ({k+1} rows × {len(topo)+1} cols).")

    def _populate_alignment_tab(self, result: AlignmentResult) -> None:
        md_widget = self.query_one("#alignment-md", Markdown)
        md_widget.update(_render_alignment_md(result))
        tabs = self.query_one("#tabs", TabbedContent)
        tabs.active = "tab-align"

if __name__ == "__main__":
    app = Seq2GraphApp()
    app.run()