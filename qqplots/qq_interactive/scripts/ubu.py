#!/usr/bin/env python3

import argparse
from pathlib import Path

import numpy as np
import pandas as pd
import scipy.stats as st


def compute_u_B_for_row(row, alpha):
    """
    Compute u and B(u) for a single row, using the formulas from the slides:

    L(n, q) = 1 - F_{chi^2_q}( (n^{1/q} / (n^{1/q} - 1)) * ln n )

    u = (p - L) / (alpha - L)

    B(u) = -e * u * log(u), if u < e^{-1}
           1, otherwise
    """
    p = row["p_value"]
    n = row["n"]

    # q: parameter dimension; if missing/invalid, default to 1
    q = row.get("q", np.nan)
    if pd.isna(q) or q <= 0:
        q = 1.0
    else:
        q = float(q)

    # Guard against non-positive n
    if n <= 0:
        return pd.Series({"u": np.nan, "B_u": np.nan})

    # Compute L(n, q)
    t = n ** (1.0 / q)
    lam = (t - 1.0) / t  # same lambda as in eJAB definition
    if lam <= 0:
        return pd.Series({"u": np.nan, "B_u": np.nan})

    arg = (1.0 / lam) * np.log(n)
    L = 1.0 - st.chi2.cdf(arg, df=q)

    denom = alpha - L
    if denom <= 0:
        # Outside the theoretical region where this transformation makes sense
        return pd.Series({"u": np.nan, "B_u": np.nan})

    u = (p - L) / denom

    # Numerical safety: clip u into (0,1)
    u_clipped = np.clip(u, 1e-15, 1 - 1e-15)

    # Compute B(u)
    e_inv = 1.0 / np.e
    if u_clipped < e_inv:
        B_u = -np.e * u_clipped * np.log(u_clipped)
    else:
        B_u = 1.0

    return pd.Series({"u": u, "B_u": B_u})


def add_u_B_columns(input_path: Path, output_path: Path, alpha: float) -> None:
    """
    Read a flagged CSV (simulation or CTG) with columns at least:
        p_value, ejab_value, n, q
    and append columns:
        u, B_u
    using the formulas from the slides.
    """
    df = pd.read_csv(input_path)

    required = {"p_value", "ejab_value", "n"}
    missing = required - set(df.columns)
    if missing:
        raise ValueError(
            f"Input file {input_path} is missing required columns: {missing}"
        )

    # Apply row-wise computation of u and B(u)
    ub = df.apply(lambda r: compute_u_B_for_row(r, alpha), axis=1)
    df["u"] = ub["u"]
    df["B_u"] = ub["B_u"]

    output_path.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(output_path, index=False)
    print(f"Wrote file with u and B(u) to {output_path}")


def main():
    parser = argparse.ArgumentParser(
        description=(
            "Add u and B(u) columns to flagged eJAB data, "
            "for either simulation or CTG datasets."
        )
    )

    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument(
        "--simulation",
        action="store_true",
        help="Use simulation flagged data: data/ejab_simulation_flagged.csv",
    )
    group.add_argument(
        "--ctg",
        action="store_true",
        help="Use CTG flagged data: data/CTG_clean_flagged.csv",
    )

    parser.add_argument(
        "--alpha",
        type=float,
        default=0.05,
        help="Significance level alpha used in defining candidates (default: 0.05).",
    )
    parser.add_argument(
        "--sim-file",
        type=str,
        default="data/ejab_simulation_flagged.csv",
        help="Path to flagged simulation CSV (default: data/ejab_simulation_flagged.csv).",
    )
    parser.add_argument(
        "--ctg-file",
        type=str,
        default="data/CTG_clean_flagged.csv",
        help="Path to flagged CTG CSV (default: data/CTG_clean_flagged.csv).",
    )
    parser.add_argument(
        "--output",
        type=str,
        default=None,
        help=(
            "Output CSV path. If not provided, defaults to:\n"
            "  data/ejab_simulation_flagged_with_u_B.csv  (for --simulation)\n"
            "  data/CTG_clean_flagged_with_u_B.csv        (for --ctg)"
        ),
    )

    args = parser.parse_args()

    if args.simulation:
        input_path = Path(args.sim_file)
        default_output = Path("data/ejab_simulation_flagged_with_u_B.csv")
    else:  # --ctg
        input_path = Path(args.ctg_file)
        default_output = Path("data/CTG_clean_flagged_with_u_B.csv")

    output_path = Path(args.output) if args.output is not None else default_output

    add_u_B_columns(input_path=input_path, output_path=output_path, alpha=args.alpha)


if __name__ == "__main__":
    main()
