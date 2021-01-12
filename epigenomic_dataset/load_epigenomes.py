from typing import Tuple
import os
from downloaders import BaseDownloader
import pandas as pd


def load_epigenomes(
    cell_line: str = "K562",
    assembly: str = "hg38",
    dataset: str = "fantom",
    region: str = "promoters",
    metric: str = "mean",
    window_size: int = 256,
    root: str = "datasets",
    verbose: int = 2
) -> Tuple[pd.DataFrame, pd.DataFrame]:
    """Return epigenomic data and labels for given parameters.

    Parameters
    ----------------------------------------
    cell_line: str = "K562",
        Cell line to consider. By default K562.
        Currently available cell lines are
        listed in the repository README file.
    assembly: str,
        The genomic assembly of the data to be retrieved.
    dataset: str = "fantom",
        Dataset to consider. By default fantom.
        Currently available datasets are
        listed in the repository README file.
    region: str = "promoters",
        Region to consider. By default promoters.
        Currently available region are
        listed in the repository README file.
    metric: str = "mean",
        The metric to load.
    window_size: int = 200,
        Window size to consider. By default 200.
        Currently available window sizes are
        listed in the repository README file.
    root: str = "datasets"
        Where to store the downloaded data.
    verbose: int = 2,
        Verbosity level.

    Returns
    ----------------------------------------
    Return tuple with input and output DataFrames.
    """
    repository = "https://github.com/LucaCappelletti94/epigenomic_dataset/blob/master/preprocessed"
    get_parameter = "?raw=true"
    data_path_placeholder = "{{root}}/{dataset}/{assembly}/{window_size}/{region}/{cell_line}.csv.xz".format(
        root=root,
        dataset=dataset,
        assembly=assembly,
        window_size=window_size,
        region=region,
        cell_line=cell_line
    )
    data_path = data_path_placeholder.format(root=root)
    label_path_placeholder = "{{root}}/{dataset}/{assembly}/{window_size}/{region}.bed.xz".format(
        root=root,
        dataset=dataset,
        assembly=assembly,
        window_size=window_size,
        region=region
    )
    label_path = label_path_placeholder.format(root=root)

    downloader = BaseDownloader(target_directory=root, verbose=verbose)

    downloader.download(
        urls=data_path_placeholder.format(root=repository)+get_parameter,
        paths=data_path
    )
    downloader.download(
        urls=label_path_placeholder.format(root=repository)+get_parameter,
        paths=label_path
    )

    dtypes = {
        "chrom": "str",
        "chromStart": "int",
        "chromEnd": "int",
        "strand": "str"
    }

    X = pd.read_csv(
        data_path,
        index_col=[0, 1, 2, 3],
        header=[0, 1],
        low_memory=False,
        dtype=dtypes
    )

    X.index.rename(
        list(dtypes.keys()),
        inplace=True
    )

    X = X[[
        col
        for col in X.columns
        if metric in col
    ]]

    X = X.droplevel(1, axis=1)

    y = pd.read_csv(
        label_path,
        index_col=[0, 1, 2, 3],
        sep="\t",
        low_memory=False,
        dtype=dtypes
    ).astype(int)

    normalized_cell_line = cell_line.replace("-", "").upper()

    if normalized_cell_line not in y.columns:
        raise ValueError(
            (
                "The requested cell line {} is not present in the labels. "
                "The available cell lines are {}"
            ).format(normalized_cell_line, ", ".join(y.columns))
        )

    y = y[[normalized_cell_line]]

    return X, y