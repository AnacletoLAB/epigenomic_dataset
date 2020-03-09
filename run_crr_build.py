from crr_labels import fantom, roadmap
from epigenomic_dataset import build
import pandas as pd
from typing import List
import os
from notipy_me import Notipy

def get_bed_path(root, target):
    return "{root}/{target}/{target}.bed".format(
        root=root,
        target=target
    )

def bed_files_exist(root):
    return (
        os.path.exists(get_bed_path(root, "promoters")) and
        os.path.exists(get_bed_path(root, "enhancers"))
    )

def build_bed(
    bed: pd.DataFrame,
    root: str,
    target: str,
    cell_lines: List[str],
    workers: int = 12
):
    path = get_bed_path(root, target)
    os.makedirs(os.path.dirname(path), exist_ok=True)
    bed.to_csv(path, sep="\t", index=False)
    regions_path = "{root}/{target}/regions.bed".format(
        root=root,
        target=target
    )
    bed["name"] = [
        "{chrom}.{chromStart}.{chromEnd}".format(
            **row.to_dict()
        )
        for _, row in bed.iterrows()
    ]

    bed["score"] = 0

    bed[["chrom", "chromStart", "chromEnd", "name", "score", "strand"]].to_csv(
        regions_path,
        sep="\t",
        header=False,
        index=False
    )
    build(
        bed_path=regions_path,
        cell_lines=cell_lines,
        epigenomes_path="epigenomes",
        targets_path="{root}/{target}".format(
            root=root,
            target=target
        ),
        workers=workers
    )


with Notipy():
    cell_lines = ["A549", "GM12878", "H1", "HEK293", "HepG2", "K562"]
    cell_lines_encode = cell_lines + ["MCF-7"]
    cell_lines_fantom = cell_lines + ["MCF7"]
    cell_lines_roadmap = ["A549", "GM12878", "H1", "HepG2", "K562"]
    windows_size = 1000

    if not bed_files_exist("fantom"):
        print("Retrieving FANTOM labels")
        enhancers, promoters = fantom(
            cell_lines=cell_lines_fantom,  # list of cell lines to be considered.
            # window size to use for the various regions.
            window_size=windows_size,
            # whetever to drop the rows where no activation is detected for every rows.
            drop_always_inactive_rows=False
        )
    else:
        print("Loading FANTOM labels")
        enhancers = pd.read_csv(get_bed_path("fantom", "enhancers"), sep="\t")
        promoters = pd.read_csv(get_bed_path("fantom", "promoters"), sep="\t")

    assert (enhancers.chromEnd - enhancers.chromStart == windows_size).all()
    assert (promoters.chromEnd - promoters.chromStart == windows_size).all()

    build_bed(
        enhancers,
        root="fantom",
        target="enhancers",
        cell_lines=cell_lines_encode
    )

    build_bed(
        promoters,
        root="fantom",
        target="promoters",
        cell_lines=cell_lines_encode
    )

    
    if not bed_files_exist("roadmap"):
        print("Retrieving ROADMAP labels")
        enhancers, promoters = roadmap(
            cell_lines=cell_lines_roadmap,  # List of cell lines to be considered.
            # Window size to use for the various regions.
            window_size=windows_size,
        )
    else:
        print("Loading ROADMAP labels")
        enhancers = pd.read_csv(get_bed_path("roadmap", "enhancers"), sep="\t")
        promoters = pd.read_csv(get_bed_path("roadmap", "promoters"), sep="\t")

    assert (enhancers.chromEnd - enhancers.chromStart == windows_size).all()
    assert (promoters.chromEnd - promoters.chromStart == windows_size).all()

    build_bed(
        enhancers,
        root="roadmap",
        target="enhancers",
        cell_lines=cell_lines_roadmap
    )

    build_bed(
        promoters,
        root="roadmap",
        target="promoters",
        cell_lines=cell_lines_roadmap
    )
