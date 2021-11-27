"""Submodule relative to cell lines utilities."""
from typing import List


def get_window_sizes() -> List[int]:
    """Return list of the available window sizes."""
    return [
        64, 128, 256, 512, 1024
    ]
