from glob import glob
from multiprocessing import Pool, cpu_count
import os
from typing import Dict, List
from tqdm.auto import tqdm
import numpy as np
import pandas as pd
from tabulate import tabulate
import csv
import gzip
import warnings
from .extract import load_epigenomes_table, load_accession_path


def compute_header(statistics: Dict[str, bool]) -> str:
    """Return the header for extracted epigenomes.

    Parameters
    ------------------
    statistics:Dict[str, bool]:
        The statistics to be computed

    Returns
    ------------------
    The header separated with tabs.
    """
    return "\t".join([
        "chrom", "chromStart", "chromEnd", "strand",
        *[s for s, enabled in statistics.items() if enabled]
    ])+'\n'


def get_callback(statistic: str):
    return {
        "mean": np.nanmean,
        "var": np.nanvar,
        "max": np.nanmax,
        "min": np.nanmin,
        "median": np.nanmedian
    }[statistic]


def get_target_path(root: str, cell_line: str, assembly: str, assay_term_name: str, target: str) -> str:
    """Return path where the target epigenomic data are to be stored.

    Parameters
    -----------------------
    root: str,
        Root from where to search the files
        to concatenate.
    cell_line: str,
        Cell line to consider.
    assembly: str,
        The genomic assembly of the data to be retrieved.
    assay_term_name: str,
        Name of the experiment.
    target: str,
        The name of the genomic target.

    Returns
    -----------------------
    The path where to store the epigenomic data.
    """
    return "{root}/{cell_line}/{file_name}.csv.gz".format(
        root=root,
        cell_line=cell_line,
        assembly=assembly,
        file_name=assay_term_name if target == "Unknown" else target
    )


def parse_extracted_epigenome(sources: List[str], target: str, statistics: Dict[str, bool]):
    """Parse the given source bed-like file.

    Parameters
    ----------------------------
    sources: List[str],
        Paths from where to load the sources.
    target: str,
        Epigenomic data target.
    statistics: Dict[str, bool]
        Statistics to be extracted.
    """
    header = compute_header(statistics)
    callbacks = [
        get_callback(s)
        for s, enabled in statistics.items()
        if enabled
    ]
    os.makedirs(os.path.dirname(target), exist_ok=True)

    source_files = [
        gzip.open(source, "rt")
        for source in sources
    ]

    readers = [
        csv.reader(source_file, delimiter='\t')
        for source_file in source_files
    ]
    try:
        with gzip.open(target, "wt") as t:
            # Starting by writing the head
            t.write(header)
            # And now we parse the lines one by one
            for rows in zip(*readers):
                # We extract the values
                chrom, chromStart, chromEnd, _, _, strand = rows[0][:6]
                # Convert the scores to float values
                scores = [
                    [
                        float(s) if s != "NA" else np.nan
                        for s in row[7:]
                    ]
                    for row in rows
                ]
                # Compute the averages in a fully defined way
                averaged_scores = [
                    np.nan
                    if len(sub_scores) == 0 or np.all(np.isnan(sub_scores))
                    else np.nanmean(sub_scores)
                    for sub_scores in scores
                ]
                # Compute the metrics
                metrics = [
                    str(np.nan)
                    if len(averaged_scores) == 0 or np.all(np.isnan(averaged_scores))
                    else cal(averaged_scores).astype(str)
                    for cal in callbacks
                ]
                # And write the results
                t.write("\t".join([
                    chrom, chromStart, chromEnd, strand,
                    *metrics
                ])+'\n')
    except EOFError:
        warnings.warn((
            "Unable to properly finish reading corrupted compressed files {}. "
            "I am now deleting these files. "
            "Just rerun the pipeline to retrieve them again."
        ).format(
            ", ".join(sources)
        ))
        for source in sources:
            os.remove(source)

    for source_file in source_files:
        source_file.close()


def _parse_extracted_epigenome(kwargs: Dict):
    return parse_extracted_epigenome(**kwargs)


def mine(
    root: str,
    statistics: Dict[str, bool],
    cell_lines: List[str],
    assembly: str
):
    """Extract and saves requested statistics from epigenomic files.

    Parameters
    -----------------------
    root: str,
        Root from where to search the files
        to concatenate.
    statistics: Dict[str, bool],
        Dictionary of statistics to extract from windows.
    cell_line: str,
        Cell line to consider.
    assembly: str,
        The genomic assembly of the data to be retrieved.
    """
    tasks = [
        {
            "sources": [
                load_accession_path(root, accession)
                for accession in group.accession
            ],
            "target": get_target_path(root, cell_line, assembly, assay_term_name, target),
            "statistics": statistics
        }
        for (cell_line, assay_term_name, target), group in load_epigenomes_table(cell_lines, assembly).groupby(["cell_line", "assay_term_name", "target"])
        if not os.path.exists(get_target_path(root, assembly, cell_line, assay_term_name, target))
    ]

    with Pool(cpu_count()) as p:
        list(tqdm(
            p.imap(
                _parse_extracted_epigenome,
                tasks
            ),
            desc="Parse extracted epigenomes",
            total=len(tasks)
        ))
        p.close()
        p.join()
