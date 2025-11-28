#!/usr/bin/env python3

import argparse
from pathlib import Path

import pandas as pd


def flag_candidates(
    input_path: Path,
    output_path: Path,
    alpha: float,
) -> None:
    """
    Read a CSV with columns:
        p_value, ejab_value, is_null, n, q, effect_size
    and write out rows where:
        p_value <= alpha  and  ejab_value >= 1
    """
    df = pd.read_csv(input_path)

    if "p_value" not in df.columns or "ejab_value" not in df.columns:
        raise ValueError(
            f"Input file {input_path} must contain 'p_value' and 'ejab_value' columns."
        )

    mask = (df["p_value"] <= alpha) & (df["ejab_value"] > 1.0)
    flagged = df.loc[mask].copy()

    output_path.parent.mkdir(parents=True, exist_ok=True)
    flagged.to_csv(output_path, index=False)
    print(f"Wrote {flagged.shape[0]} flagged rows to {output_path}")


def main():
    parser = argparse.ArgumentParser(
        description=(
            "Flag candidate Type I errors where p_value <= alpha "
            "and ejab_value >= 1, for either simulation or CTG data."
        )
    )

    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument(
        "--simulation",
        action="store_true",
        help="Use simulation data: data/ejab_simulation_data.csv",
    )
    group.add_argument(
        "--ctg",
        action="store_true",
        help="Use cleaned CTG data: data/CTG_clean.csv",
    )

    parser.add_argument(
        "--alpha",
        type=float,
        default=0.05,
        help="Significance level cutoff (default: 0.05).",
    )
    parser.add_argument(
        "--sim-file",
        type=str,
        default="data/ejab_simulation_data.csv",
        help="Path to simulation input CSV (default: data/ejab_simulation_data.csv).",
    )
    parser.add_argument(
        "--ctg-file",
        type=str,
        default="data/CTG_clean.csv",
        help="Path to cleaned CTG input CSV (default: data/CTG_clean.csv).",
    )
    parser.add_argument(
        "--output",
        type=str,
        default=None,
        help=(
            "Output CSV path. If not provided, defaults to:\n"
            "  data/ejab_simulation_flagged.csv  (for --simulation)\n"
            "  data/CTG_clean_flagged.csv        (for --ctg)"
        ),
    )

    args = parser.parse_args()

    if args.simulation:
        input_path = Path(args.sim_file)
        default_output = Path("data/ejab_simulation_flagged.csv")
    else:  # --ctg
        input_path = Path(args.ctg_file)
        default_output = Path("data/CTG_clean_flagged.csv")

    output_path = Path(args.output) if args.output is not None else default_output

    flag_candidates(input_path=input_path, output_path=output_path, alpha=args.alpha)


if __name__ == "__main__":
    main()
