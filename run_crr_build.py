from crr_labels import fantom, roadmap
from epigenomic_dataset import build
import pandas as pd
from typing import List
import os
from tqdm.auto import tqdm
from notipy_me import Notipy


def get_bed_path(root, region, window_size):
    return "{root}/{window_size}/{region}.bed".format(
        root=root,
        window_size=window_size,
        region=region
    )


def bed_files_exist(root, window_size):
    return (
        os.path.exists(get_bed_path(root, "promoters", window_size)) and
        os.path.exists(get_bed_path(root, "enhancers", window_size))
    )


def run_pipeline(
    bed: pd.DataFrame,
    root: str,
    region: str,
    windows_size: int,
    cell_lines: List[str]
):
    path = get_bed_path(root, region, windows_size)
    os.makedirs(os.path.dirname(path), exist_ok=True)
    bed.to_csv(path, sep="\t", index=False)
    regions_path = "{root}/{windows_size}/{region}/regions.bed".format(
        root=root,
        windows_size=windows_size,
        region=region
    )
    os.makedirs(os.path.dirname(regions_path), exist_ok=True)

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
        targets_path="{root}/{windows_size}/{region}".format(
            root=root,
            windows_size=windows_size,
            region=region
        ),
        extraction_workers=12,
        concatenation_workers=20
    )


if __name__ == "__main__":
    with Notipy() as r:
        cell_lines = ["A549", "GM12878", "H1", "HEK293", "HepG2", "K562"]
        cell_lines_encode = cell_lines + ["MCF-7"]
        cell_lines_fantom = cell_lines + ["MCF7"]
        cell_lines_roadmap = ["A549", "GM12878", "H1", "HepG2", "K562"]
        windows_sizes = (1000, 500, 300, 200, 100)

        for windows_size in tqdm(windows_sizes, desc="Parsing window sizes"):
            ####################################################
            # HERE WE BUILD FANTOM                             #
            ####################################################

            if not bed_files_exist("fantom", windows_size):
                print("Retrieving FANTOM labels")
                enhancers, promoters = fantom(
                    # list of cell lines to be considered.
                    cell_lines=cell_lines_fantom,
                    # window size to use for the various regions.
                    window_size=windows_size,
                    # whetever to drop the rows where no activation is detected for every rows.
                    drop_always_inactive_rows=False
                )
            else:
                print("Loading FANTOM labels")
                enhancers = pd.read_csv(get_bed_path("fantom", "enhancers", windows_size), sep="\t")
                promoters = pd.read_csv(get_bed_path("fantom", "promoters", windows_size), sep="\t")

            run_pipeline(
                enhancers,
                root="fantom",
                region="enhancers",
                windows_size=windows_size,
                cell_lines=cell_lines_encode
            )
            run_pipeline(
                promoters,
                root="fantom",
                region="promoters",
                windows_size=windows_size,
                cell_lines=cell_lines_encode
            )
            
            r.add_report{
                "window_size":windows_size,
                "dataset":"fantom"
            })

            ####################################################
            # HERE WE BUILD ROADMAP                            #
            ####################################################

            if not bed_files_exist("roadmap", windows_size):
                print("Retrieving ROADMAP labels")
                enhancers, promoters = roadmap(
                    # List of cell lines to be considered.
                    cell_lines=cell_lines_roadmap,
                    # Window size to use for the various regions.
                    window_size=windows_size,
                )
            else:
                print("Loading ROADMAP labels")
                enhancers = pd.read_csv(get_bed_path("roadmap", "enhancers", windows_size), sep="\t")
                promoters = pd.read_csv(get_bed_path("roadmap", "promoters", windows_size), sep="\t")
            
            run_pipeline(
                enhancers,
                root="roadmap",
                region="enhancers",
                windows_size=windows_size,
                cell_lines=cell_lines_roadmap
            )
            run_pipeline(
                promoters,
                root="roadmap",
                region="promoters",
                windows_size=windows_size,
                cell_lines=cell_lines_roadmap
            )

            r.add_report({
                "window_size":windows_size,
                "dataset":"roadmap"
            })