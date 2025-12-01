#!/usr/bin/env python3

import pandas as pd
from pathlib import Path

# Input: JLP05-style data (CTG subset with model, I, R, C, etc.)
INPUT_FILE = "data/JLP05.csv"

# Output: same schema as simulated data, but with proper q (df) inferred
OUTPUT_FILE = "data/CTG_clean.csv"

df = pd.read_csv(INPUT_FILE)

# Drop any unnamed index column (e.g., the leading "" column)
df = df.loc[:, ~df.columns.str.match(r"^Unnamed")]

# ----------------------------------------------------------------------
# Infer q (k, df) from model / I / R / C, matching the R infer_k logic
# ----------------------------------------------------------------------
def infer_q(row):
    model = row.get("model", None)

    # Normalize model to string safely
    model = str(model) if not pd.isna(model) else ""

    I = row.get("I", pd.NA)
    R = row.get("R", pd.NA)
    C = row.get("C", pd.NA)

    # t-test / regression / rank tests / Cox: q = 1
    if model in {
        "t-test",
        "linear_regression",
        "logistic_regression",
        "cox",
        "wilcoxon",
        "mann_whitney",
    }:
        return 1

    # ANOVA-like: q = I - 1, if I present
    if model in {"anova", "kruskal_wallis", "repeated_measures"}:
        if pd.isna(I):
            return pd.NA
        return int(I) - 1

    # Chi-squared: q = (R - 1) * (C - 1), if R and C present
    if model == "chi_squared":
        if pd.isna(R) or pd.isna(C):
            return pd.NA
        return (int(R) - 1) * (int(C) - 1)

    # Unknown model → cannot infer
    return pd.NA


df["q"] = df.apply(infer_q, axis=1)

# Optional: use pandas nullable integer dtype for q
try:
    df["q"] = df["q"].astype("Int64")
except TypeError:
    # If some values are non-integer / NA, leave as-is
    pass

# ----------------------------------------------------------------------
# Build a DataFrame with the same columns as the simulated data:
# p_value, ejab_value, is_null, n, q, effect_size
# ----------------------------------------------------------------------
clean_df = pd.DataFrame(
    {
        "p_value": df["pValue"],
        "ejab_value": df["JAB"],
        # Ground truth is unknown for real CTG data → leave as missing
        "is_null": pd.NA,
        # Sample size
        "n": df["N"],
        # Parameter dimension q inferred above
        "q": df["q"],
        # No effect size in CTG dataset → leave as missing
        "effect_size": pd.NA,
    }
)

output_path = Path(OUTPUT_FILE)
output_path.parent.mkdir(parents=True, exist_ok=True)
clean_df.to_csv(output_path, index=False)