import csv
from collections import defaultdict
from io import StringIO
from pathlib import Path
from typing import Any

from .frame import Block


def to_dict_of_list(list_of_dict: list[dict[str, Any]]) -> dict[str, list[Any]]:
    result = defaultdict(list)
    for item in list_of_dict:
        for key, value in item.items():
            result[key].append(value)
    return dict(result)


def to_list_of_dict(list_of_dict: dict[str, list[Any]]) -> list[dict[str, Any]]:
    keys = list(list_of_dict.keys())
    length = len(next(iter(list_of_dict.values())))
    return [{key: list_of_dict[key][i] for key in keys} for i in range(length)]


def read_csv(file: Path | StringIO, delimiter: str = ",") -> Block:
    """
    Read a CSV file or StringIO object and return a Block object.

    Args:
        file: Path to the CSV file or a StringIO object containing CSV data.
        delimiter: Delimiter used in the CSV file (default is comma).

    Returns:
        Block: A Block object containing the data from the CSV.
    """
    if isinstance(file, StringIO):
        file.seek(0)  # Ensure we read from the start
        reader = csv.DictReader(file, delimiter=delimiter)
        data = [row for row in reader]
    else:
        file = Path(file)
        if not file.exists():
            raise FileNotFoundError(f"File {file} does not exist.")

        with open(file, "r", newline="") as csvfile:
            reader = csv.DictReader(csvfile, delimiter=delimiter)
            data = [row for row in reader]

    # Convert list of dicts to dict of lists for Block
    dict_of_lists = to_dict_of_list(data)

    # Handle empty data case - create empty arrays for each column
    if not dict_of_lists and isinstance(file, StringIO):
        file.seek(0)
        reader = csv.DictReader(file, delimiter=delimiter)
        fieldnames = reader.fieldnames or []
        dict_of_lists = {field: [] for field in fieldnames}
    elif not dict_of_lists and not isinstance(file, StringIO):
        with open(file, "r", newline="") as csvfile:
            reader = csv.DictReader(csvfile, delimiter=delimiter)
            fieldnames = reader.fieldnames or []
            dict_of_lists = {field: [] for field in fieldnames}

    return Block(dict_of_lists)
