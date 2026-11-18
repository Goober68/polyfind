"""The row runner's shared checkpoint and the split of its rows (examples/rows_under_preset.py)."""

import importlib.util
import json
import multiprocessing
from pathlib import Path

SPEC = importlib.util.spec_from_file_location("rows_under_preset", Path(__file__).resolve().parents[1] / "examples" / "rows_under_preset.py")
RUP = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(RUP)


def test_the_parts_take_every_row_once():
    rows = [f"row{i}" for i in range(35)]
    parts = [RUP.part_of(rows, f"{k}/3") for k in (1, 2, 3)]
    assert sorted(sum(parts, [])) == sorted(rows)
    assert max(len(p) for p in parts) - min(len(p) for p in parts) <= 1
    assert RUP.part_of(rows, "") == rows


def _add(args):
    path, key = args
    RUP.stored(Path(path), lambda d: d.setdefault("rows", {}).update({key: {"cells": {}}}))


def test_runs_side_by_side_keep_each_others_cells(tmp_path):
    path = tmp_path / "checkpoint.json"
    with multiprocessing.get_context("fork").Pool(4) as pool:
        pool.map(_add, [(str(path), f"row{i}") for i in range(40)])
    assert sorted(json.loads(path.read_text())["rows"]) == sorted(f"row{i}" for i in range(40))
