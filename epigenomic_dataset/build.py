from .extract import extract
from .mine import mine
from .concatenate import concatenate
from typing import List


def build(
    bed_path: str,
    cell_lines: List[str],
    epigenomes_path: str = "epigenomes",
    targets_path: str = "targets",
    clear_download: bool = False,
    extraction_workers=-1,
    concatenation_workers=-1,
    mine_max: bool = True,
    mine_min: bool = False,
    mine_mean: bool = False,
    mine_median: bool = False,
    mine_variance: bool = False
):
    extract(
        bed_path,
        cell_lines,
        epigenomes_path=epigenomes_path,
        targets_path=targets_path,
        clear_download=clear_download,
        workers=extraction_workers
    )
    mine(
        targets_path,
        {
            "max": mine_max,
            "min": mine_min,
            "mean": mine_mean,
            "median": mine_median,
            "variance": mine_variance,
        },
        cell_lines
    )
    concatenate(
        targets_path,
        cell_lines,
        concatenation_workers
    )
