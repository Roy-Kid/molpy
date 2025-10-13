# author: Roy Kid
# contact: lijichen365@126.com
# date: 2023-01-29
# version: 0.0.1

import os
from pathlib import Path

import molq
import pytest

_REPO_URL = "https://github.com/molcrafts/tests-data.git"
_DEFAULT_DIR = Path(__file__).parent / "tests-data"


@pytest.fixture(scope="session", name="TEST_DATA_DIR")
def find_test_data() -> Path:
    """
    Ensure the tests-data repository is present and up-to-date.
    * If the directory already contains a `.git` folder → `git pull`.
      Otherwise clone afresh.
    """
    local_submitor = molq.LocalSubmitor("download-test-data")
    if (_DEFAULT_DIR / ".git").exists():
        local_submitor.local_submit(
            cmd=["git", "pull", "--ff-only"],
            cwd=_DEFAULT_DIR,
            job_name="update-test-data",
            block=True,
            quiet=True,
            cleanup_temp_files=True,
        )
    else:
        _DEFAULT_DIR.parent.mkdir(parents=True, exist_ok=True)
        local_submitor.local_submit(
            cmd=["git", "clone", "--depth", "1", _REPO_URL, str(_DEFAULT_DIR)],
            cwd=_DEFAULT_DIR.parent,
            job_name="download-test-data",
            block=True,
            quiet=True,
        )
    return _DEFAULT_DIR
