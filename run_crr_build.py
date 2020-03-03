from crr_labels import fantom, roadmap
from epigenomic_dataset import build
from notipy_me import Notipy


with Notipy():
    cell_lines = ["A549", "GM12878", "H1", "HEK293", "HepG2", "K562"]
    cell_lines_encode = cell_lines + ["MCF-7"]
    cell_lines_fantom = cell_lines + ["MCF7"]
    cell_lines_roadmap = ["A549", "GM12878", "H1", "HepG2", "K562"]
    windows_size = 1000

    # print("Downloading labels from FANTOM")
    # enhancers, promoters = fantom(
    #     cell_lines=cell_lines_fantom, # list of cell lines to be considered.
    #     window_size=windows_size, # window size to use for the various regions.
    #     genome = "hg19", # considered genome version. Currently supported only "hg19".
    #     center_enhancers = "peak", # how to center the enhancer window, either around "peak" or the "center" of the region.
    #     enhancers_threshold = 0, # activation threshold for the enhancers.
    #     promoters_threshold = 5, # activation threshold for the promoters.
    #     drop_always_inactive_rows = False, # whetever to drop the rows where no activation is detected for every rows.
    #     binarize = True, # Whetever to return the data binary-encoded, zero for inactive, one for active.
    #     nrows = None # the number of rows to read, usefull when testing pipelines for creating smaller datasets.
    # )

    # enhancers.to_csv("fantom_enhancers.bed", sep="\t")
    # promoters.to_csv("fantom_promoters.bed", sep="\t")

    # enhancers[["chrom","chromStart","chromEnd"]].to_csv(
    #     "enhancers_regions.bed",
    #     sep="\t",
    #     header=False,
    #     index=False
    # )
    # print("Parsing enhancers epigenomes")
    # build(
    #     bed_path="enhancers_regions.bed",
    #     cell_lines=cell_lines_encode,
    #     nan_threshold=1,
    #     epigenomes_path="epigenomes",
    #     targets_path="fantom_enhancers",
    #     workers=4
    # )

    # promoters[["chrom","chromStart","chromEnd"]].to_csv(
    #     "promoters_regions.bed",
    #     sep="\t",
    #     header=False,
    #     index=False
    # )

    # print("Parsing promoters epigenomes")
    # build(
    #     bed_path="promoters_regions.bed",
    #     cell_lines=cell_lines_encode,
    #     nan_threshold=1,
    #     epigenomes_path="epigenomes",
    #     targets_path="fantom_promoters",
    #     workers=4
    # )

    print("Downloading labels from ROADMAP")
    enhancers, promoters = roadmap(
        cell_lines = cell_lines_roadmap, # List of cell lines to be considered.
        window_size = windows_size, # Window size to use for the various regions.
        genome = "hg19", # Considered genome version. Currently supported only "hg19".
        states = 18, # Number of the states of the model to consider. Currently supported only "15" and "18".
        enhancers_labels = ("7_Enh", "9_EnhA1", "10_EnhA2"), # Labels to encode as active enhancers.
        promoters_labels = ("1_TssA",), # Labels to enode as active promoters.
    )

    enhancers.to_csv("roadmap_enhancers.bed", sep="\t")
    promoters.to_csv("roadmap_promoters.bed", sep="\t")

    enhancers[["chrom","chromStart","chromEnd"]].to_csv(
        "roadmap_enhancers_regions.bed",
        sep="\t",
        header=False,
        index=False
    )
    print("Parsing enhancers epigenomes")
    build(
        bed_path="roadmap_enhancers_regions.bed",
        cell_lines=cell_lines_encode,
        nan_threshold=1,
        epigenomes_path="epigenomes",
        targets_path="roadmap_enhancers",
        workers=1
    )

    promoters[["chrom","chromStart","chromEnd"]].to_csv(
        "roadmap_promoters_regions.bed",
        sep="\t",
        header=False,
        index=False
    )

    print("Parsing promoters epigenomes")
    build(
        bed_path="roadmap_promoters_regions.bed",
        cell_lines=cell_lines_encode,
        nan_threshold=1,
        epigenomes_path="epigenomes",
        targets_path="roadmap_promoters",
        workers=1
    )