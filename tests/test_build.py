import os
from epigenomic_dataset import build
import shutil
import pytest


def test_build():
    for assembly in ("hg19", "hg38"):
        build("tests/test.bed", ["GM12892"], assembly)
        with pytest.warns(UserWarning):
            build("tests/test.bed", ["GM12892"], assembly)
        shutil.rmtree("targets")
