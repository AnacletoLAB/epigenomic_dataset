import os
from epigenomic_dataset import build
import shutil

def test_build():
    build("tests/test.bed", ["GM12892"])
    build("tests/test.bed", ["GM12892"])
    shutil.rmtree("targets")