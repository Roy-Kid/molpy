# author: Roy Kid
# contact: lijichen365@126.com
# date: 2023-01-29
# version: 0.0.1

from pathlib import Path
import pytest

@pytest.fixture(name="test_data_path", scope="session")
def find_test_data() -> Path:

    data_path = Path(__file__).parent / "chemfile-testcases"
    if not data_path.exists():
        pytest.skip("test data repository unavailable")
    return data_path

