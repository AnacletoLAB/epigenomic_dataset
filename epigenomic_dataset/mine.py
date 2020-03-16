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
        "var": np.nanmean,
        "max": np.nanmax,
        "min": np.nanmin,
        "median": np.nanmedian
    }[statistic]


def get_target_path(root: str, cell_line: str, target: str):
    return "{root}/{cell_line}/{target}.csv.gz".format(
        root=root,
        cell_line=cell_line,
        target=target
    )


def parse_extracted_epigenome(sources: str, target: str, statistics: Dict[str, bool]):
    """Parse the given source bed-like file."""
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

    with gzip.open(target, "wt") as t:
        # Starting by writing the head
        t.write(header)
        # And now we parse the lines one by one
        for rows in zip(*readers):
            # We extract the values
            chrom, chromStart, chromEnd, _, _, strand = rows[0][:6]
            # Convert the scores to float values
            scores = np.nanmean([
                [
                    float(s) if s != "NA" else np.nan
                    for s in row[7:]
                ]
                for row in rows
            ], axis=0)
            metrics = [
                cal(scores).astype(str)
                for cal in callbacks
            ]
            # And write the results
            t.write("\t".join([
                chrom, chromStart, chromEnd, strand,
                *metrics
            ])+'\n')

    for source_file in source_files:
        source_file.close()


def _parse_extracted_epigenome(kwargs: Dict):
    return parse_extracted_epigenome(**kwargs)


def mine(root: str, statistics: Dict[str, bool], cell_lines: List[str]):
    tasks = [
        {
            "sources": [
                load_accession_path(root, accession)
                for accession in group.accession
            ],
            "target": get_target_path(root, cell_line, target),
            "statistics": statistics
        }
        for (cell_line, target), group in load_epigenomes_table(cell_lines).groupby(["cell_line", "target"])
        if not os.path.exists(get_target_path(root, cell_line, target))
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
