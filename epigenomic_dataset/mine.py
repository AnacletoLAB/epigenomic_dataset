from glob import glob
from multiprocessing import Pool, cpu_count
import os
from typing import Dict
from tqdm.auto import tqdm
import numpy as np
import pandas as pd
from tabulate import tabulate
import csv
import gzip


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
        "mean": np.mean,
        "var": np.mean,
        "max": np.max,
        "min": np.min,
        "median": np.median
    }[statistic]

def get_target_path(source:str):
    root, filename = os.path.split(source)
    return "{root}/parsed/{filename}".format(
        root=root,
        filename=filename
    )

def parse_extracted_epigenome(source: str, statistics: Dict[str, bool]):
    """Parse the given source bed-like file."""
    target = get_target_path(source)

    header = compute_header(statistics)
    callbacks = [
        get_callback(s)
        for s, enabled in statistics.items()
        if enabled
    ]
    nans = ["nan"]*len(callbacks)
    os.makedirs(os.path.dirname(target), exist_ok=True)
    with gzip.open(source, "rt") as s:
        with gzip.open(target, "wt") as t:
            # Opening file reader
            reader = csv.reader(s, delimiter='\t')
            # Starting by writing the head
            t.write(header)
            # And now we parse the lines one by one
            for row in reader:
                # We extract the values
                chrom, chromStart, chromEnd, _, _, strand = row[:6]
                # Convert the scores to float values
                scores = np.array([
                    float(s)
                    for s in row[7:]
                    if s != "NA"
                ])
                if scores.size != 0:
                    metrics = [
                        cal(scores).astype(str)
                        for cal in callbacks
                    ]
                else:
                    metrics = nans
                # And write the results
                t.write("\t".join([
                    chrom, chromStart, chromEnd, strand,
                    *metrics
                ])+'\n')


def _parse_extracted_epigenome(kwargs: Dict):
    return parse_extracted_epigenome(**kwargs)


def mine(root: str, statistics: Dict[str, bool]):
    tasks = [
        {
            "source": source,
            "statistics": statistics
        }
        for source in glob(f"{root}/*.bed.gz")
        if not os.path.exists(get_target_path(source))
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
