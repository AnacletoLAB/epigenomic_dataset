import os
from epigenomic_dataset import build
import shutil
import pytest

def test_build():
    build("tests/test.bed", ["GM12892"])
    with pytest.warns(UserWarning): 
        build("tests/test.bed", ["GM12892"])
    shutil.rmtree("targets")