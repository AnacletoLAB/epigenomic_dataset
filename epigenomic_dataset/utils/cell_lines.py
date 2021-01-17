"""Submodule relative to cell lines utilities."""
from typing import List


def get_cell_lines() -> List[str]:
    """Return list of the available cell lines."""
    return [
        "A549",
        "GM12878",
        "H1",
        "HEK293",
        "HepG2",
        "K562",
        "MCF-7"
    ]
