from crr_labels import fantom
from epigenomic_dataset import build
from notipy_me import Notipy


with Notipy():
    cell_lines_encode = [
        "A549",
        "GM12878",
        "H1",
        "HEK293",
        "HepG2",
        "K562",
        "MCF-7"
    ]
    cell_lines_fantom = [
        "A549",
        "GM12878",
        "H1",
        "HEK293",
        "HepG2",
        "K562",
        "MCF7"
    ]
    windows_size = 1000
    path_to_bed = "regions.bed"

    enhancers, promoters = fantom(
        cell_lines=cell_lines_fantom, # list of cell lines to be considered.
        window_size=windows_size, # window size to use for the various regions.
        genome = "hg19", # considered genome version. Currently supported only "hg19".
        center_enhancers = "peak", # how to center the enhancer window, either around "peak" or the "center" of the region.
        enhancers_threshold = 0, # activation threshold for the enhancers.
        promoters_threshold = 5, # activation threshold for the promoters.
        drop_always_inactive_rows = False, # whetever to drop the rows where no activation is detected for every rows.
        binarize = True, # Whetever to return the data binary-encoded, zero for inactive, one for active.
        nrows = None # the number of rows to read, usefull when testing pipelines for creating smaller datasets.
    )

    enhancers.to_csv("enhancers.bed", sep="\t")
    promoters.to_csv("promoters.bed", sep="\t")

    enhancers[[
        "chrom",
        "chromStart",
        "chromEnd"
    ]].to_csv("regions.bed", sep="\t")

    build(
        bed_path=path_to_bed,
        cell_lines=cell_lines_encode
    )