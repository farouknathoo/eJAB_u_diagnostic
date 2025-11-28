# Interactive QQ plot

The procedure so far:
- `clean_jlp.py` reformats the CTG data (JLP05.csv) to be the same format as the simulation data
- `flag.py` outputs .csv of flagged data for CTG (run w flag `--ctg`) or simulated (`--simulation`) data (must use a flag else does not run)
- `ubu.py` computes u and B(u) for either flagged sets, outputting .csv files containing corresponding rows
- `qqunif.r`


Running in this order will recreate all data files. Optionally, the .zip files have been included for ease.