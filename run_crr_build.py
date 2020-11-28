from crr_labels import fantom, roadmap
from epigenomic_dataset import build
import pandas as pd
from typing import List
import os
from tqdm.auto import tqdm
from notipy_me import Notipy


def get_bed_path(root: str, assembly: str, region: str, window_size: int) -> str:
    return "{root}/{assembly}/{window_size}/{region}.bed".format(
        root=root,
        assembly=assembly,
        window_size=window_size,
        region=region
    )


def bed_files_exist(root: str, assembly: str, window_size: int):
    return (
        os.path.exists(get_bed_path(root, assembly, "promoters", window_size)) and
        os.path.exists(get_bed_path(root, assembly, "enhancers", window_size))
    )


def run_pipeline(
    bed: pd.DataFrame,
    root: str,
    assembly: str,
    region: str,
    windows_size: int,
    cell_lines: List[str]
):
    path = get_bed_path(root, assembly, region, windows_size)
    os.makedirs(os.path.dirname(path), exist_ok=True)
    bed.to_csv(path, sep="\t", index=False)
    regions_path = "{root}/{assembly}/{windows_size}/{region}/regions.bed".format(
        root=root,
        assembly=assembly,
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
        assembly=assembly,
        targets_path="{root}/{assembly}/{windows_size}/{region}".format(
            root=root,
            assembly=assembly,
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
        windows_sizes = (1024, 512, 256, 128, 64)
        assembly = "hg38"
        # We are not computing RoadMap right now
        # because we are still choosing the states from the model to be used.
        build_roadmap = False

        for windows_size in tqdm(windows_sizes, desc="Parsing window sizes"):
            ####################################################
            # HERE WE BUILD FANTOM                             #
            ####################################################

            if not bed_files_exist("fantom", assembly, windows_size):
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
                enhancers = pd.read_csv(get_bed_path(
                    "fantom", assembly, "enhancers", windows_size), sep="\t")
                promoters = pd.read_csv(get_bed_path(
                    "fantom", assembly, "promoters", windows_size), sep="\t")

            run_pipeline(
                enhancers,
                root="fantom",
                assembly=assembly,
                region="enhancers",
                windows_size=windows_size,
                cell_lines=cell_lines_encode
            )
            run_pipeline(
                promoters,
                root="fantom",
                assembly=assembly,
                region="promoters",
                windows_size=windows_size,
                cell_lines=cell_lines_encode
            )

            r.add_report({
                "window_size": windows_size,
                "dataset": "fantom"
            })

            if build_roadmap:
                ####################################################
                # HERE WE BUILD ROADMAP                            #
                ####################################################

                if not bed_files_exist("roadmap", assembly, windows_size):
                    print("Retrieving ROADMAP labels")
                    enhancers, promoters = roadmap(
                        # List of cell lines to be considered.
                        cell_lines=cell_lines_roadmap,
                        # Window size to use for the various regions.
                        window_size=windows_size,
                    )
                else:
                    print("Loading ROADMAP labels")
                    enhancers = pd.read_csv(get_bed_path(
                        "roadmap", assembly, "enhancers", windows_size), sep="\t")
                    promoters = pd.read_csv(get_bed_path(
                        "roadmap", assembly, "promoters", windows_size), sep="\t")

                run_pipeline(
                    enhancers,
                    root="roadmap",
                    assembly=assembly,
                    region="enhancers",
                    windows_size=windows_size,
                    cell_lines=cell_lines_roadmap
                )
                run_pipeline(
                    promoters,
                    root="roadmap",
                    assembly=assembly,
                    region="promoters",
                    windows_size=windows_size,
                    cell_lines=cell_lines_roadmap
                )

                r.add_report({
                    "window_size": windows_size,
                    "dataset": "roadmap"
                })
