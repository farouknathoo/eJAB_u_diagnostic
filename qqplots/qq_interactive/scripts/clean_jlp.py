#!/usr/bin/env python3

import pandas as pd
from pathlib import Path

# Input: CTG-style data (example snippet you provided)
INPUT_FILE = "data/JLP05.csv"

# Output: same schema as simulated data
OUTPUT_FILE = "data/CTG_clean.csv"

df = pd.read_csv(INPUT_FILE)

# Drop any unnamed index column (e.g., the leading "" column)
df = df.loc[:, ~df.columns.str.match(r"^Unnamed")]

# Build a DataFrame with the same columns as the simulated data:
# p_value, ejab_value, is_null, n, q, effect_size
clean_df = pd.DataFrame(
    {
        "p_value": df["pValue"],
        "ejab_value": df["JAB"],
        # Ground truth is unknown for real CTG data → leave as missing
        "is_null": pd.NA,
        # Sample size
        "n": df["N"],
        # Parameter dimension q is not given here → leave as missing
        "q": pd.NA,
        # No effect size in CTG dataset → leave as missing
        "effect_size": pd.NA,
    }
)

output_path = Path(OUTPUT_FILE)
output_path.parent.mkdir(parents=True, exist_ok=True)
clean_df.to_csv(output_path, index=False)
