# author: Roy Kid
# contact: lijichen365@126.com
# date: 2023-01-29
# version: 0.0.1

import os
import subprocess
from pathlib import Path

import molq
import pytest

_REPO_URL = "https://github.com/molcrafts/chemfile-testcases.git"
_DEFAULT_DIR = Path(__file__).parent / "chemfile-testcases"


@pytest.fixture(scope="session", name="TEST_DATA_DIR")
def find_test_data() -> Path:
    """
    Ensure the chemfile-testcases repository is present and up-to-date.

    * Set env var **CHEMFILE_DATA_DIR** to override the clone location.
    * If the directory already contains a `.git` folder → `git pull`.
      Otherwise clone afresh.
    """
    data_dir = Path(os.getenv("CHEMFILE_DATA_DIR", _DEFAULT_DIR)).expanduser()
    local_submitor = molq.LocalSubmitor("download-test-data")
    if (data_dir / ".git").exists():
        local_submitor.local_submit(
            cmd=["git", "pull", "--ff-only"],
            cwd=data_dir,
            job_name="update-test-data",
            block=True,
            quiet=True,
        )
    else:
        data_dir.parent.mkdir(parents=True, exist_ok=True)
        local_submitor.local_submit(
            cmd=["git", "clone", "--depth", "1", _REPO_URL, str(data_dir)],
            cwd=data_dir.parent,
            job_name="download-test-data",
            block=True,
            quiet=True,
        )
    return data_dir
