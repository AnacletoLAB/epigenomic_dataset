"""Submodule relative to cell lines utilities."""
from typing import List


def get_available_metrics() -> List[str]:
    """Return list of the available metrics."""
    return ["mean", "max"]
