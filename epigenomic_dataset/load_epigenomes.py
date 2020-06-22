from typing import Tuple
import os
from encodeproject.utils import download
import pandas as pd


def load_epigenomes(
    cell_line: str = "K562",
    dataset: str = "fantom",
    regions: str = "promoters",
    window_size: int = 200,
    root: str = "datasets",
    drop_unique_group_by: bool = True
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
    drop_unique_group_by: bool = True,
        Whetever to drop the group by column layer
        when it is unique, eg when only a max group by
        is available in the dataset.

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
    label_path_placeholder = "{{root}}/{dataset}/{window_size}/{regions}.bed.gz".format(
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

    X.columns.names = (None,)*len(X.columns.names)

    if X.columns.levels[1].size == 1 and drop_unique_group_by:
        X = X.droplevel(1, axis=1)

    y = pd.read_csv(
        label_path,
        index_col=[0, 1, 2, 3],
        sep="\t",
        dtype=dtypes
    ).astype(int)

    y = y[[cell_line.replace("-", "").upper()]]

    return X, y
