import pandas as pd
import os
from tqdm.auto import tqdm
from typing import List
from glob import glob
from multiprocessing import Pool, cpu_count
from .extract import load_epigenomes_table
from .mine import get_target_path


def to_dict(path: str):
    feature_name = path.split(os.sep)[-1].split(".")[0]
    df = pd.read_csv(path, sep="\t", index_col=[0, 1, 2, 3]).round(2)
    df.columns = pd.MultiIndex.from_product([[feature_name], df.columns])
    return df


def concatenate(
    root: str,
    cell_lines: List[str],
    workers: int
):
    """Concatenate the targets into a single file.

    Parameters
    ----------------------------------
    root: str,
        Root from where to search the files
        to concatenate.
    cell_lines: List[str],
        Cell lines to consider.
    workers: int,
        Workers to use to parallelize the loading.
    """
    if workers == -1:
        workers = cpu_count()
    with Pool(min(cpu_count(), workers)) as p:
        for cell_line in tqdm([
            directory_name
            for directory_name in os.listdir(root)
            if os.path.isdir(f"{root}/{directory_name}")
        ], leave=False, desc="Concatenating cell lines"):
            path = "{root}/{cell_line}.csv.xz".format(
                root=root,
                cell_line=cell_line
            )
            paths = glob("{root}/{cell_line}/*.csv.gz".format(
                root=root,
                cell_line=cell_line
            ))

            if os.path.exists(path):
                continue

            pd.concat(list(tqdm(
                p.imap(to_dict, paths),
                desc="Concatenating files",
                total=len(paths),
                leave=False
            )), axis=1).reset_index().to_csv(path, index=False)

        p.close()
        p.join()
