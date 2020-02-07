import pandas as pd
import os
from typing import List
from encodeproject import download
from pybwtool import extract
from tqdm.auto import tqdm


def build(
    bed_path: str,
    cell_lines: List[str],
    nan_threshold: float = 0.9,
    path: str = "epigenomes",
    clear_download: bool = False
):
    """Download bigwigs from ENCODE and extract regions specified in given bed file.

    Parameters
    --------------------------
    bed_path:str,
        Either path to the bed file containing the regions of interest.
    cell_lines:List[str],
        List of cell lines whose epigenomes are to retrieve.
    nan_threshold: float = 0.9,
        Percentage of NaNs to allow in every region.
    path:str="epigenomes",
        Path where to store the epigenomic data
    clear_download: bool = False,
        Whetever to delete the downloaded files or not.
        By default False.
    """
    # Loading the epigenomes metadata
    epigenomes = pd.read_csv(
        "{}/epigenomes.csv".format(os.path.dirname(os.path.abspath(__file__)))
    )
    # Filtering epigenomes for required cell lines
    epigenomes = epigenomes[epigenomes.cell_line.isin(cell_lines)]
    # Creating target directory if doesn't exist already
    os.makedirs(path, exist_ok=True)
    # Downloading and elaborating data
    for _, row in tqdm(epigenomes.iterrows(), total=len(epigenomes), desc="Downloading epigenomes"):
        # Where to store the downloaded bigWig file
        epigenome_path = "{path}/{accession}.{file_format}".format(
            path=path,
            **row.to_dict()
        )
        # Whete to store the extracted regions
        target_path = "{path}/{accession}.bed".format(
            path=path,
            **row.to_dict()
        )
        # If the file was already parsed
        if not os.path.exists(target_path):
            # Download file
            download(row.url, epigenome_path)
            # Extract the features
            bed, scores = extract(
                bed_path=bed_path,
                bigwig_path=epigenome_path,
                nan_threshold=nan_threshold
            )
            # Save the obtained features
            pd.concat([bed, scores], axis=1).to_csv(target_path, sep="\t")
        # Remove the bigwig file
        if clear_download:
            os.remove(epigenome_path)
