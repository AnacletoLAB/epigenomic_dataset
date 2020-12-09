import pandas as pd
import os
from typing import List
from encodeproject import download
from pybwtool import extract as extract_bigwig
from tqdm.auto import tqdm
from multiprocessing import Pool, cpu_count
import warnings


def extraction_job(
    bed_path: str,
    epigenome_path: str,
    target_path: str,
    url: str,
    clear_download: bool
):
    # Download file if it does not already exist
    if not os.path.exists(epigenome_path):
        download(url, epigenome_path)

    # Extract the features
    extract_bigwig(
        bed_path=bed_path,
        bigwig_path=epigenome_path,
        target=target_path
    )

    # Remove the bigwig file if required
    if clear_download:
        os.remove(epigenome_path)


def _extraction_job(kwargs):
    return extraction_job(**kwargs)


def load_epigenomes_table(
    cell_lines: List[str],
    assembly: str,
) -> pd.DataFrame:
    """Return epigenomic data table.

    Parameters
    ----------------------------
    cell_lines:List[str],
        List of cell lines whose epigenomes are to retrieve.
    assembly: str,
        The genomic assembly of the data to be retrieved.

    Returns
    ----------------------------
    Pandas DataFrame with epigenomes meta data.
    """
    # Loading the epigenomes metadata
    epigenomes = pd.read_csv(
        "{pwd}/epigenomes_metadata/{assembly}.csv".format(
            pwd=os.path.dirname(os.path.abspath(__file__)),
            assembly=assembly
        )
    )
    # Filtering epigenomes for required cell lines
    return epigenomes[epigenomes.cell_line.isin(cell_lines)]


def load_accession_path(root: str, accession: str) -> str:
    return "{root}/{accession}.bed.gz".format(
        root=root,
        accession=accession
    )


def build_extraction_tasks(
    bed_path: str,
    cell_lines: List[str],
    assembly: str,
    epigenomes_path: str,
    target_path: str,
    clear_download: bool
) -> List:
    """Download bigwigs from ENCODE and extract regions specified in given bed file.

    Parameters
    --------------------------
    bed_path:str,
        Either path to the bed file containing the regions of interest.
    cell_lines:List[str],
        List of cell lines whose epigenomes are to retrieve.
    assembly: str,
        The genomic assembly of the data to be retrieved.
    epigenomes_path:str,
        Path where to store the epigenomic bigWig.
    target_path:str,
        Path where to store the epigenomic bed.
    clear_download: bool = False,
        Whetever to delete the downloaded files or not.
        By default False.

    Returns
    --------------------------
    Returns list of tasks to be executed.
    """
    # Loading the epigenomes metadata
    epigenomes = load_epigenomes_table(cell_lines, assembly)
    # Build the tasks
    return [
        {
            "bed_path": bed_path,
            # Where to store the downloaded bigWig file
            "epigenome_path": "{epigenomes_path}/{accession}.{file_format}".format(
                epigenomes_path=epigenomes_path,
                **epigenome.to_dict()
            ),
            # Where to store the extracted regions
            "target_path": load_accession_path(target_path, epigenome.accession),
            "url": epigenome.url,
            "clear_download": clear_download
        }
        for _, epigenome in epigenomes.iterrows()
        if not os.path.exists("{target_path}/{accession}.bed.gz".format(target_path=target_path, **epigenome.to_dict()))
    ]


def extract(
    bed_path: str,
    cell_lines: List[str],
    assembly: str,
    epigenomes_path: str = "epigenomes",
    targets_path: str = "targets",
    clear_download: bool = False,
    workers: int = -1
):
    """Download bigwigs from ENCODE and extract regions specified in given bed file.

    Parameters
    --------------------------
    bed_path:str,
        Either path to the bed file containing the regions of interest.
    cell_lines:List[str],
        List of cell lines whose epigenomes are to retrieve.
    assembly: str,
        The genomic assembly of the data to be retrieved.
    epigenomes_path:str,
        Path where to store the epigenomic bigWig.
    targets_path:str,
        Path where to store the epigenomic bed.
    clear_download: bool = False,
        Whetever to delete the downloaded files or not.
        By default False.
    workers: int = -1,
        Number of workers to use.
        The default, -1, set the workers to the maximum available.

    Raises
    --------------------------
    ValueError,
        If workers number is neiter -1 nor a strictly positive integer.
    ValueError,
        If given nan threshold is not a float value between 0 and 1.
    """
    # Creating target directory if doesn't exist already
    os.makedirs(epigenomes_path, exist_ok=True)
    os.makedirs(targets_path, exist_ok=True)
    # Create the building tasks list
    tasks = build_extraction_tasks(
        bed_path,
        cell_lines,
        assembly,
        epigenomes_path,
        targets_path,
        clear_download
    )
    # Set workers number
    if workers == -1:
        workers = cpu_count()
    workers = min(workers, cpu_count(), len(tasks))
    if len(tasks) == 0:
        warnings.warn("No epigenome still to be parsed, moving on.",
                      category=UserWarning)
        return
    if workers < 1:
        raise ValueError(
            "Given workers number {} is neither -1 or a strictly positive integer.".format(workers))
    # Downloading and elaborating data
    with Pool(workers) as p:
        list(tqdm(
            p.imap(_extraction_job, tasks),
            total=len(tasks),
            desc="Parsing epigenomes"
        ))
        p.close()
        p.join()
