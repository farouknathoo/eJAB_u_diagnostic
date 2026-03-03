# eJAB Type I Error Analysis (Adaptive C*(alpha))

## Settings

- Alpha: 0.05
- Up: 0.1
- Studies: 500
- C grid: [0, 1] with 200 points

## C*(alpha) Estimation

- C*(alpha = 0.05): 0.0101
- Method: pointwise grid search minimising (proportion - alpha/up)^2

## Detection Results

- Predicted T1Es: 53 / 500 (10.6%)

## Ground Truth Comparison

- True T1Es (null & p < alpha): 16 / 500 (3.2%)
- Correctly predicted: 16 / 16 (100.0%)
- False positives: 37
- Sensitivity: 1.000
- Specificity: 0.924

## Plots

- `calibration_new_output.pdf` page 1 -- Calibration curve using adaptive C*(alpha)
- `calibration_new_output.pdf` page 2 -- C*(alpha) vs alpha
- `calibration_new_output.pdf` page 3 -- Diagnostic QQ-plot with best-fit line

## Candidates

See `candidates.csv` (53 candidates)
