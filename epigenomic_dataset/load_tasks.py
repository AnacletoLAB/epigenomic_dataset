from typing import Tuple
import pandas as pd
from tqdm.auto import tqdm
from .load_epigenomes import load_epigenomes


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
        promoters_epi = promoters_epi[promoters_labels.to_numpy() == 1]
        enhancers_epi = enhancers_epi[enhancers_labels.to_numpy() == 1]
        promoters_labels = promoters_labels[promoters_labels.to_numpy() == 1]
        enhancers_labels = enhancers_labels[enhancers_labels.to_numpy() == 1]
        enhancers_labels[enhancers_labels.columns[0]] = 0
    elif only_inactive:
        promoters_epi = promoters_epi[promoters_labels.to_numpy() == 0]
        enhancers_epi = enhancers_epi[enhancers_labels.to_numpy() == 0]
        promoters_labels = promoters_labels[promoters_labels.to_numpy() == 0]
        promoters_labels[promoters_labels.columns[0]] = 1
        enhancers_labels = enhancers_labels[enhancers_labels.to_numpy() == 0]

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


def active_promoters_vs_inactive_promoters(
    cell_line: str = "K562",
    assembly: str = "hg38",
    dataset: str = "fantom",
    metric: str = "mean",
    window_size: int = 256,
    root: str = "datasets",
    verbose: int = 2,
):
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

    Returns
    ----------------------------------------
    Return tuple with input and output DataFrames.
    """
    return load_task(
        cell_line=cell_line,
        assembly=assembly,
        dataset=dataset,
        metric=metric,
        window_size=window_size,
        root=root,
        verbose=verbose,
        only_active=False,
        only_inactive=False,
        load_enhancers=False,
        load_promoters=True,
    )


def active_enhancers_vs_inactive_enhancers(
    cell_line: str = "K562",
    assembly: str = "hg38",
    dataset: str = "fantom",
    metric: str = "mean",
    window_size: int = 256,
    root: str = "datasets",
    verbose: int = 2,
):
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

    Returns
    ----------------------------------------
    Return tuple with input and output DataFrames.
    """
    return load_task(
        cell_line=cell_line,
        assembly=assembly,
        dataset=dataset,
        metric=metric,
        window_size=window_size,
        root=root,
        verbose=verbose,
        only_active=False,
        only_inactive=False,
        load_enhancers=True,
        load_promoters=False,
    )


def active_vs_inactive(
    cell_line: str = "K562",
    assembly: str = "hg38",
    dataset: str = "fantom",
    metric: str = "mean",
    window_size: int = 256,
    root: str = "datasets",
    verbose: int = 2,
):
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

    Returns
    ----------------------------------------
    Return tuple with input and output DataFrames.
    """
    return load_task(
        cell_line=cell_line,
        assembly=assembly,
        dataset=dataset,
        metric=metric,
        window_size=window_size,
        root=root,
        verbose=verbose,
        only_active=False,
        only_inactive=False,
        load_enhancers=True,
        load_promoters=True,
    )


def active_enhancers_vs_active_promoters(
    cell_line: str = "K562",
    assembly: str = "hg38",
    dataset: str = "fantom",
    metric: str = "mean",
    window_size: int = 256,
    root: str = "datasets",
    verbose: int = 2,
):
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

    Returns
    ----------------------------------------
    Return tuple with input and output DataFrames.
    """
    return load_task(
        cell_line=cell_line,
        assembly=assembly,
        dataset=dataset,
        metric=metric,
        window_size=window_size,
        root=root,
        verbose=verbose,
        only_active=True,
        only_inactive=False,
        load_enhancers=True,
        load_promoters=True,
    )


def inactive_enhancers_vs_inactive_promoters(
    cell_line: str = "K562",
    assembly: str = "hg38",
    dataset: str = "fantom",
    metric: str = "mean",
    window_size: int = 256,
    root: str = "datasets",
    verbose: int = 2,
):
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

    Returns
    ----------------------------------------
    Return tuple with input and output DataFrames.
    """
    return load_task(
        cell_line=cell_line,
        assembly=assembly,
        dataset=dataset,
        metric=metric,
        window_size=window_size,
        root=root,
        verbose=verbose,
        only_active=False,
        only_inactive=True,
        load_enhancers=True,
        load_promoters=True,
    )


def load_all_tasks(
    cell_line: str = "K562",
    assembly: str = "hg38",
    dataset: str = "fantom",
    metric: str = "mean",
    window_size: int = 256,
    root: str = "datasets",
    verbose: int = 2,
    leave: bool = False
):
    """Return generator with all the tasks.

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
    leave: bool = False,
        Wether to leave the loading bar.

    Returns
    ----------------------------------------
    Return tuple with input and output DataFrames.
    """
    return (
        (
            task(
                cell_line=cell_line,
                assembly=assembly,
                dataset=dataset,
                metric=metric,
                window_size=window_size,
                root=root,
                verbose=verbose,
            ),
            task.__name__
        )
        for task in tqdm(
            (
                active_enhancers_vs_inactive_enhancers,
                active_promoters_vs_inactive_promoters,
                active_enhancers_vs_active_promoters,
                inactive_enhancers_vs_inactive_promoters,
                active_vs_inactive
            ),
            desc="Executing CRR prediction tasks",
            disable=verbose==0,
            leave=leave
        )
    )
