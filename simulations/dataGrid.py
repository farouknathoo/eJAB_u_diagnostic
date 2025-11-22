"""
Grid-based simulation study for eJAB diagnostic performance analysis.

This script generates simulations across a systematic grid of:
- π₀ (null proportion): controls the mixture of true nulls vs alternatives
- n (sample size): affects power and diagnostic performance
- d (effect size): Cohen's d for alternatives

Output is used for publication-quality diagnostic and ROC plots.
"""

import sys
import argparse
import numpy as np
import pandas as pd
import scipy.stats as st
from pathlib import Path
from itertools import product

# Reproducible RNG
RNG_SEED = 277


# Edited to have consistent num_sums across all n
def get_num_sims(n: int) -> int:
    """Reduce simulation count for large n to manage memory."""
    if n <= 200:
        return 10_000
    elif n <= 1000:
        return 10_000
    elif n <= 5000:
        return 10_000
    else:
        return 10_000


def ejab_statistic(n: int, q_dim: int, p_values: np.ndarray) -> np.ndarray:
    """
    Compute eJAB_01 statistic from p-values.
    
    eJAB_01 = sqrt(n) * exp{ -0.5 * ((n^(1/q) - 1) / n^(1/q)) * Q_{chi^2_q}(1 - p) }
    where Q_{chi^2_q} is the chi-square quantile with df = q_dim.
    
    Parameters
    ----------
    n : int
        Sample size per study
    q_dim : int
        Degrees of freedom for chi-square distribution
    p_values : np.ndarray
        Array of p-values
        
    Returns
    -------
    np.ndarray
        Array of eJAB statistics
    """
    # Guard against extreme p leading to 1 - p outside (0,1)
    one_minus_p = 1.0 - np.clip(p_values, 1e-300, 1.0 - 1e-16)
    chi2_q = st.chi2.ppf(one_minus_p, df=q_dim)
    lam = (n**(1.0 / q_dim) - 1.0) / (n**(1.0 / q_dim))
    return np.sqrt(n) * np.exp(-0.5 * lam * chi2_q)


def simulate_grid_cell(
    n: int,
    pi0: float,
    effect_size: float,
    q_dim: int,
    num_sims: int,
    rng: np.random.Generator
) -> pd.DataFrame:
    """
    Generate simulations for a single grid cell (n, π₀, d combination).
    
    Parameters
    ----------
    n : int
        Sample size per study
    pi0 : float
        Proportion of true nulls (0-1)
    effect_size : float
        Cohen's d for alternative hypothesis
    q_dim : int
        Degrees of freedom for chi-square distribution
    num_sims : int
        Number of simulated studies
    rng : np.random.Generator
        Random number generator
        
    Returns
    -------
    pd.DataFrame
        Simulation results with columns: p_value, ejab_value, is_null, n, q, pi0, effect_size
    """
    # Process in chunks to manage memory
    chunk_size = min(num_sims, 10_000)
    num_chunks = int(np.ceil(num_sims / chunk_size))
    
    results_list = []
    
    for i in range(num_chunks):
        start = i * chunk_size
        end = min((i + 1) * chunk_size, num_sims)
        chunk_n = end - start
        
        # Assign truth based on π₀
        is_null_chunk = rng.random(chunk_n) < pi0
        true_effects = np.where(is_null_chunk, 0.0, effect_size)
        
        # Each study: n observations from N(true_effect, 1)
        data = rng.normal(loc=true_effects[:, None], scale=1.0, size=(chunk_n, n))
        
        # Two-sided one-sample t-test against 0 for each study
        p_chunk = st.ttest_1samp(data, 0.0, axis=1, alternative="two-sided").pvalue
        
        # eJAB_01 statistic
        ejab_chunk = ejab_statistic(n, q_dim, p_chunk)
        
        # Store results
        df_chunk = pd.DataFrame({
            "p_value": p_chunk,
            "ejab_value": ejab_chunk,
            "is_null": is_null_chunk.astype(int),
            "n": n,
            "q": q_dim,
            "pi0": pi0,
            "effect_size": effect_size
        })
        results_list.append(df_chunk)
        
        # Free memory
        del data
    
    return pd.concat(results_list, ignore_index=True)


def run_grid_simulations(
    pi0_values: list[float],
    n_values: list[int],
    effect_sizes: list[float],
    q_dim: int = 1,
    output_file: str = "data/ejab_grid_simulation.csv",
    seed: int = RNG_SEED,
    scale: float = 1.0,
) -> pd.DataFrame:
    """
    Run simulations across a grid of (π₀, n, d) combinations.
    Uses adaptive simulation counts based on sample size (see get_num_sims).
    
    Parameters
    ----------
    pi0_values : list[float]
        List of null proportions to test
    n_values : list[int]
        List of sample sizes to test
    effect_sizes : list[float]
        List of effect sizes to test
    q_dim : int, default=1
        Degrees of freedom for chi-square distribution
    output_file : str, default="data/ejab_grid_simulation.csv"
        Path to save output CSV file
    seed : int, default=RNG_SEED
        Random seed for reproducibility
        
    Returns
    -------
    pd.DataFrame
        Combined dataframe with all simulation results
    """
    rng = np.random.default_rng(seed)
    output_path = Path(output_file)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    
    # Generate all grid combinations
    grid = list(product(pi0_values, n_values, effect_sizes))
    total_cells = len(grid)
    
    # Calculate total simulations using adaptive counts
    total_sims = int(sum(get_num_sims(n) for _, n, _ in grid) * scale)
    
    print(f"Running {total_cells} grid cells...")
    print(f"  π₀ values: {pi0_values}")
    print(f"  n values: {n_values}")
    # Clean formatting for effect sizes to avoid verbose numpy float repr
    es_fmt = [f"{es:.2f}" for es in effect_sizes]
    print(f"  effect sizes: {es_fmt}")
    print(f"  Simulations per cell: ADAPTIVE (100k for n≤200, 50k for n≤1000, 20k for n≤5000, 10k for n>5000) × scale={scale}")
    print(f"  Total simulations: {total_sims:,}")
    print()
    
    all_results = []
    
    for idx, (pi0, n, d) in enumerate(grid, 1):
        num_sims = max(1, int(get_num_sims(n) * scale))
        print(f"[{idx}/{total_cells}] π₀={pi0:.2f}, n={n}, d={d:.2f} ({num_sims:,} sims)...", end=" ")
        
        df_cell = simulate_grid_cell(
            n=n,
            pi0=pi0,
            effect_size=d,
            q_dim=q_dim,
            num_sims=num_sims,
            rng=rng
        )
        
        all_results.append(df_cell)
        print("✓")
    
    # Combine all results
    combined_df = pd.concat(all_results, ignore_index=True)
    combined_df.to_csv(output_path, index=False)
    
    print(f"\n{'='*70}")
    print(f"All results saved to: {output_path}")
    print(f"Total rows: {len(combined_df):,}")
    print(f"{'='*70}")
    
    return combined_df


if __name__ == "__main__":
    # Parse command line arguments
    parser = argparse.ArgumentParser(description="Run grid-based eJAB simulation study")
    parser.add_argument("--clean", action="store_true", 
                       help="Delete previously generated data before running")
    parser.add_argument(
        "--dense-n",
        action="store_true",
        help="Use a dense sample size grid 40..2000 (step 40) and save to *_dense_n.csv",
    )
    parser.add_argument(
        "--n-step",
        type=int,
        default=40,
        help="Step size for dense-n grid (default 40)",
    )
    parser.add_argument(
        "--dense-es",
        action="store_true",
        help="Use a dense effect size grid [0.0, 0.1, ..., 1.0]",
    )
    parser.add_argument(
        "--es",
        type=str,
        default=None,
        help="Comma-separated list of effect sizes to use (e.g., '0.2,0.5'). Overrides other ES settings.",
    )
    parser.add_argument(
        "--scale",
        type=float,
        default=1.0,
        help="Scale factor for simulations per cell (e.g., 0.1 = 10% of default)",
    )
    args = parser.parse_args()
    
    # Define parameter grid
    PI0_VALUES = [0.50, 0.75, 0.90]  # 50%, 75%, 90% nulls
    # Default (sparse) sample sizes retained from original study design
    base_n_values = [40, 80, 120, 200, 500, 1000, 2000, 5000, 10000, 20000]
    if args.dense_n:
        # Dense n grid with configurable step
        step = max(1, int(args.n_step))
        N_VALUES = list(range(40, 2000 + step, step))
    else:
        N_VALUES = base_n_values
    # Effect sizes: standardized effect size (μ/σ) where σ=1.0
    if args.es:
        try:
            EFFECT_SIZES = [float(x) for x in args.es.split(',') if x.strip() != '']
        except ValueError:
            print("Error: --es must be a comma-separated list of numbers, e.g., --es 0.2,0.5")
            sys.exit(2)
    elif args.dense_es:
        EFFECT_SIZES = [round(x, 2) for x in np.linspace(0.0, 1.0, 11)]
    else:
        EFFECT_SIZES = [0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.50, 0.80]
    # Output file naming based on density flags
    if args.es is not None and args.dense_n:
        OUTPUT_FILE = "data/ejab_grid_simulation_dense_n_es_custom.csv"
    elif args.dense_n and args.dense_es:
        OUTPUT_FILE = "data/ejab_grid_simulation_dense_n_es.csv"
    elif args.dense_n:
        OUTPUT_FILE = "data/ejab_grid_simulation_dense_n.csv"
    elif args.dense_es:
        OUTPUT_FILE = "data/ejab_grid_simulation_dense_es.csv"
    else:
        OUTPUT_FILE = "data/ejab_grid_simulation.csv"
    Q_DIM = 1
    
    # Clean old data if requested
    if args.clean:
        output_path = Path(__file__).parent / OUTPUT_FILE
        if output_path.exists():
            print(f"Cleaning: Removing {output_path}")
            output_path.unlink()
        else:
            print(f"No existing data file to clean: {output_path}")
        print()
    
    print("=" * 70)
    print("eJAB GRID-BASED SIMULATION STUDY")
    print("=" * 70)
    print(f"Grid dimensions:")
    print(f"  π₀ (null proportion): {len(PI0_VALUES)} values")
    n_label = f" (dense, step={args.n_step})" if args.dense_n else ""
    print(f"  n (sample size): {len(N_VALUES)} values" + n_label)
    d_label = " (dense)" if args.dense_es else ""
    print(f"  d (effect size): {len(EFFECT_SIZES)} values" + d_label)
    print(f"  Total grid cells: {len(PI0_VALUES) * len(N_VALUES) * len(EFFECT_SIZES)}")
    print(f"q_dim: {Q_DIM}")
    print(f"scale: {args.scale}")
    print("=" * 70)
    print()
    
    # Run the simulations
    results = run_grid_simulations(
        pi0_values=PI0_VALUES,
        n_values=N_VALUES,
        effect_sizes=EFFECT_SIZES,
        q_dim=Q_DIM,
        output_file=OUTPUT_FILE,
        scale=args.scale,
    )
    
    print("\nSimulation complete!")
    print(f"\nGrid cell summary:")
    print(results.groupby(['pi0', 'n', 'effect_size']).size().head(20))
