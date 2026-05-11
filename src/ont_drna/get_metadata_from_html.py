#!/usr/bin/env python

# Standard library imports
import argparse
import csv
import json
from configparser import ConfigParser
from pathlib import Path
from re import sub
from typing import Any, Dict, List

# Third-party imports
from bs4 import BeautifulSoup
from numpy import nan


def collect_html_files(input_files: List[Path]) -> List[Path]:
    """Collect all HTML files from given input paths.

    Args:
        input_files: List of file or directory paths.

    Returns:
        List of HTML file paths.
    """
    files = []

    for input_file in input_files:
        path = Path(input_file)

        if path.is_file() and path.suffix.lower() == ".html":
            files.append(path)
        elif path.is_dir():
            files.extend(path.rglob("*.html"))

    return files


def extract_json_from_html(file_path: Path, keys_of_interest: List[str]) -> Dict[str, Any]:
    """Extract selected JSON fields from an HTML report.

    Args:
        file_path: Path to the HTML file.
        keys_of_interest: List of keys to extract from the JSON data.

    Returns:
        Dictionary containing extracted key-value pairs.
    """
    tmp_output_dict = {}
    soup = BeautifulSoup(file_path.read_text(encoding="utf-8"), "lxml")
    script_content = soup.find_all("script")[0]
    text = json.loads(script_content.string.strip().replace("const reportData=", ""))

    for key in keys_of_interest:
        if key in text and (isinstance(text[key], list)):
            for item in text[key]:
                if isinstance(item, dict) and "title" in item and "value" in item:
                    tmp_output_dict[item["title"]] = item["value"]
                else:
                    tmp_output_dict[key] = text[key]
        elif key in text and isinstance(text[key], dict):
            tmp_output_dict.update(text[key])
        elif key in text:
            tmp_output_dict[key] = text[key]

    return tmp_output_dict


def filter_dict(input_dict: Dict[str, Any], exclude_keys: List[str]) -> Dict[str, Any]:
    """Filter out keys from a dictionary.

    Args:
        data: Dictionary to filter.
        exclude_keys: List of keys to exclude.

    Returns:
        Filtered dictionary.
    """
    output_dict = {}
    for key, value in input_dict.items():
        # Remove non-alphanumeric characters, trim whitespace, and convert to lowercase with underscores
        cleaned_key = sub(r"^[^\w]+", "", key).strip()
        cleaned_key = cleaned_key.replace(" ", "_").lower()
        if isinstance(value, str):
            cleaned_value = sub(r"^[^\w]+", "", value).strip()
        else:
            cleaned_value = value
        if key == "sample_id":
            cleaned_key = "sample_name"
        if key in exclude_keys or cleaned_key in exclude_keys:
            continue
        output_dict[cleaned_key] = cleaned_value

    return output_dict


def sort_output_and_write_csv(list_col_values: List[Dict[str, Any]], output_file: Path, config: ConfigParser) -> None:
    """Sort the extracted data and write to a CSV file.

    Args:
        list_col_values: List of dictionaries containing row data.
        output_file: Path to the output CSV file.
        config: Configuration parser.
    """
    if not list_col_values:
        return

    # Collect unique column names from the data
    unique_col_names = set()
    for row_dict in list_col_values:
        unique_col_names.update(row_dict.keys())

    # Get ordered keys from config and ensure they are at the beginning of the column names
    ordered_keys_sections = config.options("ordered_keys")
    ordered_keys = []
    for section in ordered_keys_sections:
        ordered_keys.extend(json.loads(config.get("ordered_keys", section)))

    # Ensure ordered keys are at the beginning of the column names
    column_names = ordered_keys + list(unique_col_names - set(ordered_keys))
    try:
        with output_file.open(mode="w", newline="", encoding="utf-8") as file:
            writer = csv.DictWriter(file, fieldnames=column_names, restval=nan)
            writer.writeheader()
            writer.writerows(list_col_values)
    except ValueError as e:
        print(f"Error writing CSV: {e}. Used keys: {column_names}")
        raise


def parse_args() -> argparse.Namespace:
    """Parse command-line arguments.

    Returns:
        Parsed arguments namespace.
    """
    parser = argparse.ArgumentParser(description="Parse MinKNOW HTML reports.")

    parser.add_argument("config_file", type=Path, help="Path to config file")
    parser.add_argument("inputs", nargs="+", type=Path, help="HTML files or directories")
    parser.add_argument(
        "-o",
        "--output",
        type=Path,
        default=Path("output.csv"),
        help="Output CSV file",
    )

    return parser.parse_args()


def main() -> None:
    """Main entry point for the script."""
    args = parse_args()

    config = ConfigParser()
    config.read(args.config_file)

    keys_of_interest = json.loads(config.get("settings", "keys_of_interest"))
    exclude_keys = json.loads(config.get("settings", "exclude_keys"))

    html_files = collect_html_files(args.inputs)

    extracted = [extract_json_from_html(path, keys_of_interest) for path in html_files]
    filtered = [filter_dict(row, exclude_keys) for row in extracted]

    sort_output_and_write_csv(filtered, args.output, config)

    print(json.dumps(filtered, indent=2))


if __name__ == "__main__":
    main()
