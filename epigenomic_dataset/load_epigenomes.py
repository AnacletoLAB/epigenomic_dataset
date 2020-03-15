from typing import Tuple
import os
from encodeproject.utils import download
import pandas as pd


def load_epigenomes(
    cell_line: str = "K562",
    dataset: str = "fantom",
    regions: str = "promoters",
    window_size: int = 200,
    root: str = "datasets"
) -> Tuple[pd.DataFrame, pd.DataFrame]:
    """Return epigenomic data and labels for given parameters.
    
    Parameters
    ----------------------------------------
    cell_line: str = "K562",
        Cell line to consider. By default K562.
        Currently available cell lines are
        listed in the repository README file.
    dataset: str = "fantom",
        Dataset to consider. By default fantom.
        Currently available datasets are
        listed in the repository README file.
    regions: str = "promoters",
        Regions to consider. By default promoters.
        Currently available regions are
        listed in the repository README file.
    window_size: int = 200,
        Window size to consider. By default 200.
        Currently available window sizes are
        listed in the repository README file.
    root: str = "datasets"
        Where to store the downloaded data.

    Returns
    ----------------------------------------
    Return tuple with input and output DataFrames.
    """
    repository = "https://github.com/LucaCappelletti94/epigenomic_dataset/blob/master/preprocessed"
    get_parameter = "?raw=true"
    data_path_placeholder = "{{root}}/{dataset}/{window_size}/{regions}/{cell_line}.csv.gz".format(
        root=root,
        dataset=dataset,
        window_size=window_size,
        regions=regions,
        cell_line=cell_line
    )
    data_path = data_path_placeholder.format(root=root)
    label_path_placeholder = "{{root}}/{dataset}/{window_size}/{regions}.bed".format(
        root=root,
        dataset=dataset,
        window_size=window_size,
        regions=regions
    )
    label_path = label_path_placeholder.format(root=root)

    if not os.path.exists(data_path):
        download(
            url=data_path_placeholder.format(root=repository)+get_parameter,
            path=data_path
        )
    if not os.path.exists(label_path):
        download(
            url=label_path_placeholder.format(root=repository)+get_parameter,
            path=label_path
        )
    return pd.read_csv(data_path, low_memory=False), pd.read_csv(label_path, sep="\t")[[
        "chrom", "chromStart", "chromEnd", "strand", cell_line.replace("-", "")
    ]]