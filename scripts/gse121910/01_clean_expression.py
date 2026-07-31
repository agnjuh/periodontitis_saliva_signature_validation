#!/usr/bin/env python3

from pathlib import Path
import pandas as pd

INPUT = Path("data/processed/miRNA_expression.tsv")

OUT_DIR = Path("results/tables")
OUT_DIR.mkdir(parents=True, exist_ok=True)

df = pd.read_csv(INPUT, sep="\t")

print(f"Original rows: {len(df)}")

# teljesen azonos sorok eltávolítása
df = df.drop_duplicates()

print(f"After removing identical duplicates: {len(df)}")

# ugyanazon miRNA összevonása
summary = (
    df.groupby("miRNA-ID", as_index=False)
      .agg({
          "Control Raw Reads Count":"sum",
          "Case Raw Reads Count":"sum",
          "RPM-Control":"sum",
          "RPM-Case":"sum"
      })
)

print(f"Unique miRNAs: {len(summary)}")

summary.to_csv(
    OUT_DIR / "miRNA_expression_collapsed.tsv",
    sep="\t",
    index=False
)

print(summary.head())
