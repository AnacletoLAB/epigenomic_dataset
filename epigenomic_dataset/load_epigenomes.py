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
        dtype=dtypes
    ).astype(int)

    y = y[[cell_line.replace("-", "").upper()]]

    return X, y


def load_task(
    cell_line: str = "K562",
    assembly: str = "hg38",
    dataset: str = "fantom",
    metric: str = "mean",
    window_size: int = 256,
    root: str = "datasets",
    verbose: int = 2,
    only_active: bool = False,
    only_inactive: bool = False,
    load_promoters: bool = False,
    load_enhancers: bool = False
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
    only_active: bool = False,
        Wether to filter for only active.
    only_inactive: bool = False,
        Wether to filter for only inactive.
    load_promoters: bool = False,
        Wether to load promoters.
    load_enhancers: bool = False
        Wether to load enhancers.

    Returns
    ----------------------------------------
    Return tuple with input and output DataFrames.
    """
    if only_active and only_inactive:
        raise ValueError(
            "It is not possible to require both only active and inactive."
        )
    if not (load_promoters or load_enhancers):
        raise ValueError(
            "You need to load either the promoters or the enhancers."
        )
    if (only_active or only_inactive) and not (load_promoters and load_enhancers):
        raise ValueError(
            "When requiring to filter only active or inactive regions, "
            "you must load both enhancers and promoters."
        )

    (promoters_epi, promoters_labels), (enhancers_epi, enhancers_labels) = [
        load_epigenomes(
            cell_line=cell_line,
            assembly=assembly,
            dataset=dataset,
            region=region,
            metric=metric,
            window_size=window_size,
            root=root,
            verbose=verbose
        ) if enabled else (None, None)
        for region, enabled in (
            ("promoters", load_promoters),
            ("enhancers", load_enhancers),
        )
    ]
    if only_active:
        promoters_epi = promoters_epi[promoters_labels.numpy() == 1]
        enhancers_epi = enhancers_epi[enhancers_labels.numpy() == 1]
        promoters_labels = promoters_labels[promoters_labels.numpy() == 1]
        enhancers_labels = enhancers_labels[promoters_labels.numpy() == 1]
        enhancers_labels[enhancers_labels.columns[0]] = 0
    elif only_inactive:
        promoters_epi = promoters_epi[promoters_labels.numpy() == 0]
        enhancers_epi = enhancers_epi[enhancers_labels.numpy() == 0]
        promoters_labels = promoters_labels[promoters_labels.numpy() == 0]
        promoters_labels[promoters_labels.columns[0]] = 1
        enhancers_labels = enhancers_labels[promoters_labels.numpy() == 0]

    if only_active or only_inactive:
        epigenomic_data = pd.concat([
            promoters_epi,
            enhancers_epi
        ])
        labels = pd.concat([
            promoters_labels,
            enhancers_labels
        ])
    else:
        epigenomic_data = pd.concat([
            region
            for region in (promoters_epi, enhancers_epi)
            if region is not None
        ])
        labels = pd.concat([
            region
            for region in (promoters_labels, enhancers_labels)
            if region is not None
        ])

    return epigenomic_data, labels
