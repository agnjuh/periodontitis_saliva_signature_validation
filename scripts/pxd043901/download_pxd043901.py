#!/usr/bin/env python3

from pathlib import Path
import json
import subprocess
import sys
from urllib.parse import urlparse

import requests


ACCESSION = "PXD043901"
OUTPUT_DIR = Path("data/raw") / ACCESSION

TARGET_NAMES = {
    "proteinGroups.txt",
    "peptides.txt",
    "evidence.txt",
    "parameters.txt",
    "summary.txt",
}

API_URLS = [
    f"https://www.ebi.ac.uk/pride/ws/archive/v2/projects/{ACCESSION}/files?page=0&pageSize=1000",
    f"https://www.ebi.ac.uk/pride/ws/archive/v2/projects/{ACCESSION}",
]


def collect_file_records(obj):
    records = []

    if isinstance(obj, dict):
        filename = (
            obj.get("fileName")
            or obj.get("filename")
            or obj.get("name")
        )

        urls = []

        for key in (
            "downloadLink",
            "downloadUrl",
            "fileDownloadLink",
            "ftpUrl",
            "url",
        ):
            value = obj.get(key)
            if isinstance(value, str) and value.startswith(
                ("http://", "https://", "ftp://")
            ):
                urls.append(value)

        locations = obj.get("publicFileLocations")
        if isinstance(locations, list):
            for location in locations:
                if isinstance(location, dict):
                    for value in location.values():
                        if isinstance(value, str) and value.startswith(
                            ("http://", "https://", "ftp://")
                        ):
                            urls.append(value)
                elif isinstance(location, str):
                    urls.append(location)

        if filename and urls:
            records.append(
                {
                    "filename": Path(str(filename)).name,
                    "url": urls[0],
                }
            )

        for value in obj.values():
            records.extend(collect_file_records(value))

    elif isinstance(obj, list):
        for item in obj:
            records.extend(collect_file_records(item))

    return records


def download_file(url: str, destination: Path) -> None:
    destination.parent.mkdir(parents=True, exist_ok=True)

    print(f"Downloading: {destination.name}")
    print(f"From: {url}")

    command = [
        "curl",
        "-L",
        "--fail",
        "--retry",
        "3",
        "--continue-at",
        "-",
        "-o",
        str(destination),
        url,
    ]

    subprocess.run(command, check=True)


def main() -> None:
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

    response_json = None
    successful_url = None

    for api_url in API_URLS:
        print(f"Querying PRIDE API: {api_url}")

        try:
            response = requests.get(api_url, timeout=60)
            response.raise_for_status()
            response_json = response.json()
            successful_url = api_url
            break
        except Exception as error:
            print(f"API query failed: {error}")

    if response_json is None:
        raise RuntimeError(
            "Could not retrieve the PRIDE project file list."
        )

    with open(
        OUTPUT_DIR / "pride_api_response.json",
        "w",
        encoding="utf-8",
    ) as handle:
        json.dump(response_json, handle, indent=2)

    records = collect_file_records(response_json)

    unique_records = {}
    for record in records:
        unique_records[record["filename"]] = record["url"]

    print(f"\nPRIDE endpoint used: {successful_url}")
    print(f"File records found: {len(unique_records)}\n")

    for filename in sorted(unique_records):
        print(filename)

    with open(
        OUTPUT_DIR / "available_files.tsv",
        "w",
        encoding="utf-8",
    ) as handle:
        handle.write("filename\turl\n")
        for filename, url in sorted(unique_records.items()):
            handle.write(f"{filename}\t{url}\n")

    selected = {
        filename: url
        for filename, url in unique_records.items()
        if filename in TARGET_NAMES
    }

    if not selected:
        print(
            "\nNo exact MaxQuant target filenames were detected."
        )
        print(
            "Inspect data/raw/PXD043901/available_files.tsv "
            "for the actual archive names."
        )
        sys.exit(1)

    print("\nSelected processed files:")
    for filename in sorted(selected):
        print(f"  {filename}")

    for filename, url in selected.items():
        destination = OUTPUT_DIR / filename

        if destination.exists() and destination.stat().st_size > 0:
            print(f"Already present: {destination}")
            continue

        download_file(url, destination)

    print("\nDownload complete.")


if __name__ == "__main__":
    main()
