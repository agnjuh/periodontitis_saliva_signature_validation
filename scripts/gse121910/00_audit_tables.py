#!/usr/bin/env python3

from pathlib import Path
import pandas as pd


INPUTS = {
    "expression": Path("data/processed/miRNA_expression.tsv"),
    "case_control": Path("data/processed/known_miRNA_case_vs_control.tsv"),
}

OUTPUT = Path("results/tables/dataset_audit.tsv")


def audit_table(name: str, path: Path) -> dict:
    if not path.exists():
        raise FileNotFoundError(f"Missing input file: {path}")

    df = pd.read_csv(path, sep="\t")

    print("=" * 80)
    print(f"Table: {name}")
    print(f"Path: {path}")
    print(f"Shape: {df.shape[0]} rows × {df.shape[1]} columns")
    print(f"Duplicate complete rows: {df.duplicated().sum()}")
    print("\nColumns:")
    for column in df.columns:
        print(f"  - {column}")

    print("\nMissing values by column:")
    print(df.isna().sum().to_string())

    print("\nFirst five rows:")
    print(df.head().to_string(index=False))
    print()

    return {
        "table": name,
        "path": str(path),
        "rows": df.shape[0],
        "columns": df.shape[1],
        "duplicate_complete_rows": int(df.duplicated().sum()),
        "total_missing_values": int(df.isna().sum().sum()),
    }


def main() -> None:
    OUTPUT.parent.mkdir(parents=True, exist_ok=True)

    records = [
        audit_table(name, path)
        for name, path in INPUTS.items()
    ]

    pd.DataFrame(records).to_csv(OUTPUT, sep="\t", index=False)
    print(f"Audit summary written to: {OUTPUT}")


if __name__ == "__main__":
    main()
