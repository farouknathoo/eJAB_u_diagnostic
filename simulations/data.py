import numpy as np
import pandas as pd
import scipy.stats as st
from pathlib import Path

# Settings
RNG_SEED = 277
PI0 = 0.75          # proportion of true nulls
NUM_SIMS = 70_000
SAMPLE_SIZES = [15, 30, 50, 100, 200]
EFFECT_SIZES = [0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.50, 0.80]
# EFFECT_SIZES = [0.60, 0.70, 0.80]
Q_DIM = 1
OUTPUT_FILE = "data/ejab_simulation_data.csv"

rng = np.random.default_rng(RNG_SEED)

p_values = []
ejab_values = []
is_nulls = []
n_values = []
effect_values = []


def ejab_statistic(p, n, q=1):
    # clip p to avoid 0 or 1
    p = np.clip(p, 1e-300, 1.0 - 1e-16)
    chi2_q = st.chi2.ppf(1.0 - p, df=q)
    lam = (n**(1.0 / q) - 1.0) / (n**(1.0 / q))
    return np.sqrt(n) * np.exp(-0.5 * lam * chi2_q)


for _ in range(NUM_SIMS):
    n = rng.choice(SAMPLE_SIZES)
    effect = rng.choice(EFFECT_SIZES)

    # choose null vs alternative
    is_null = rng.random() < PI0
    mu = 0.0 if is_null else effect

    # generate N(mu, 1) data
    data = rng.normal(loc=mu, scale=1.0, size=n)

    # two-sided one-sample t-test against 0
    p = st.ttest_1samp(data, 0.0, alternative="two-sided").pvalue

    # eJAB_01 with q = 1
    ejab = ejab_statistic(p, n, q=Q_DIM)

    p_values.append(p)
    ejab_values.append(ejab)
    is_nulls.append(int(is_null))
    n_values.append(n)
    effect_values.append(effect)

df = pd.DataFrame({
    "p_value": p_values,
    "ejab_value": ejab_values,
    "is_null": is_nulls,
    "n": n_values,
    "q": Q_DIM,
    "effect_size": effect_values,
})

output_path = Path(OUTPUT_FILE)
output_path.parent.mkdir(parents=True, exist_ok=True)
df.to_csv(output_path, index=False)
