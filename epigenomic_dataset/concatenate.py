import pandas as pd
import os
from tqdm.auto import tqdm
from typing import List
from multiprocessing import Pool, cpu_count


def to_dict(accession: str, target:str, root: str):
    df = pd.read_csv(
        "{root}/parsed/{accession}.bed.gz".format(
            root=root,
            accession=accession
        ), sep="\t", index_col=[0, 1, 2, 3]).round(2)
    df.columns = pd.MultiIndex.from_product([[target], df.columns])
    return df


def _to_dict(kwargs):
    return to_dict(**kwargs)


def concatenate(root: str, cell_lines: List[str], workers: int):
    data_table = pd.read_csv(
        "{}/epigenomes.csv".format(os.path.dirname(os.path.abspath(__file__))))
    if workers == -1:
        workers = cpu_count()
    filtered = data_table[data_table.cell_line.isin(cell_lines)]
    for cell_line, group in tqdm(filtered.groupby("cell_line"), leave=False, desc="Concatenating cell lines"):
        path = f"{root}/{cell_line}.csv.gz"
        if os.path.exists(path):
            continue
        tasks = [
            {
                "accession": row.accession,
                "target": row.target,
                "root": root
            }
            for _, row in group.iterrows()
        ]
        with Pool(min(cpu_count(), workers)) as p:
            concatenation = pd.concat(list(tqdm(
                p.imap(
                    _to_dict,
                    tasks
                ),
                desc="Concatenating files",
                total=len(tasks),
                leave=False
            )), axis=1)
            p.close()
            p.join()

        concatenation.reset_index().to_csv(path, index=False)
