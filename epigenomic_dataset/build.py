from .extract import extract
from .mine import mine
from .concatenate import concatenate
from typing import List


def build(
    bed_path: str,
    cell_lines: List[str],
    assembly: str,
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
    """Build the dataset.

    Parameters
    --------------------------
    bed_path:str,
        Either path to the bed file containing the regions of interest.
    cell_lines:List[str],
        List of cell lines whose epigenomes are to retrieve.
    assembly: str,
        The genomic assembly of the data to be retrieved.
    epigenomes_path:str,
        Path where to store the epigenomic bigWig.
    targets_path:str,
        Path where to store the epigenomic bed.
    clear_download: bool = False,
        Whetever to delete the downloaded files or not.
        By default False.
    extraction_workers: int = -1,
        Number of workers to use.
        The default, -1, set the workers to the maximum available.
    concatenation_workers: int,
        Workers to use to parallelize the loading.
    mine_max: bool = True,
        Wether to mine max value in windows.
    mine_min: bool = False,
        Wether to mine min value in windows.
    mine_mean: bool = False,
        Wether to mine mean value in windows.
    mine_median: bool = False,
        Wether to mine median value in windows.
    mine_variance: bool = False
        Wether to mine variance value in windows.
    """
    extract(
        bed_path,
        cell_lines,
        assembly,
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
            "var": mine_variance,
        },
        cell_lines,
        assembly
    )
    concatenate(
        targets_path,
        cell_lines,
        concatenation_workers,
    )
