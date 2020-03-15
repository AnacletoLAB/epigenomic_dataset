import pandas as pd
import os
from tqdm.auto import tqdm
from typing import List
from multiprocessing import Pool, cpu_count
from .extract import load_epigenomes_table
from .mine import get_target_path


def to_dict(path: str, target: str):
    df = pd.read_csv(path, sep="\t", index_col=[0, 1, 2, 3]).round(2)
    df.columns = pd.MultiIndex.from_product([[target], df.columns])
    return df


def _to_dict(kwargs):
    return to_dict(**kwargs)


def concatenate(root: str, cell_lines: List[str], workers: int):
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
    table = load_epigenomes_table(cell_lines)
    if workers == -1:
        workers = cpu_count()
    with Pool(min(cpu_count(), workers)) as p:
        for cell_line, group in tqdm(table.groupby("cell_line"), leave=False, desc="Concatenating cell lines"):
            path = "{root}/{cell_line}.csv.gz".format(
                root=root,
                cell_line=cell_line
            )
            if os.path.exists(path):
                continue
            paths = [
                {
                    "path": get_target_path(root, cell_line, target),
                    "target": target
                }
                for target in group.target.unique()
            ]

            pd.concat(list(tqdm(
                p.imap(_to_dict, paths),
                desc="Concatenating files",
                total=len(paths),
                leave=False
            )), axis=1).reset_index().to_csv(path, index=False)

        p.close()
        p.join()
