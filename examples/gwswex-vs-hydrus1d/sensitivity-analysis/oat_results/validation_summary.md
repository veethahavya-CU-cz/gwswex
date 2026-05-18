# Validation summary — OAT-tuned comparison notebooks

## Full-horizon RMSE and NSE

| setup | soil | rc | runtime [s] | impl RMSE [cm] | impl NSE | expl RMSE [cm] | expl NSE |
|---|---|---|---|---|---|---|---|
| basic | loam | 0 | 15.4 | 0.67 | 1.00 | 1.56 | 0.99 |
| basic | sand | 0 | 18.0 | 0.66 | 1.00 | 0.29 | 1.00 |
| basic | clay | 0 | 15.8 | 0.86 | 0.99 | 3.96 | 0.89 |
| basic | sand-loam | 0 | 13.3 | 0.13 | 1.00 | 0.58 | 1.00 |
| basic | sand-clay | 0 | 18.1 | 0.55 | 1.00 | 0.35 | 1.00 |
| basic | loam-clay | 0 | 14.9 | 0.33 | 1.00 | 1.39 | 0.97 |
| intensive | loam | 0 | 8.9 | 16.45 | 0.80 | 16.00 | 0.81 |
| intensive | sand | 0 | 10.9 | 14.11 | 0.91 | 12.66 | 0.93 |
| intensive | clay | 0 | 12.3 | 37.59 | 0.25 | 41.65 | 0.08 |
| intensive | sand-loam | 0 | 7.8 | 19.00 | 0.84 | 22.77 | 0.77 |
| intensive | sand-clay | 0 | 9.0 | 11.51 | 0.92 | 11.18 | 0.92 |
| intensive | loam-clay | 0 | 9.1 | 6.14 | 0.96 | 5.70 | 0.97 |

## Phase-resolved RMSE and NSE

| setup | soil | solver | wet RMSE [cm] | wet NSE | dry RMSE [cm] | dry NSE |
|---|---|---|---|---|---|---|
| basic | loam | impl | 0.68 | 1.00 | 0.71 | 0.94 |
| basic | loam | expl | 1.53 | 0.98 | 1.71 | 0.65 |
| basic | sand | impl | 0.38 | 1.00 | 0.89 | 0.86 |
| basic | sand | expl | 0.13 | 1.00 | 0.41 | 0.97 |
| basic | clay | impl | 0.67 | 0.99 | 1.06 | 0.99 |
| basic | clay | expl | 3.21 | 0.88 | 4.88 | 0.87 |
| basic | sand-loam | impl | 0.16 | 1.00 | 0.12 | 0.99 |
| basic | sand-loam | expl | 0.79 | 1.00 | 0.31 | 0.95 |
| basic | sand-clay | impl | 0.48 | 0.99 | 0.64 | 0.94 |
| basic | sand-clay | expl | 0.43 | 0.99 | 0.27 | 0.99 |
| basic | loam-clay | impl | 0.24 | 0.99 | 0.42 | 0.98 |
| basic | loam-clay | expl | 0.42 | 0.96 | 2.00 | 0.64 |
| intensive | loam | impl | 1.21 | 1.00 | 26.10 | -1.46 |
| intensive | loam | expl | 4.27 | 0.99 | 24.82 | -1.23 |
| intensive | sand | impl | 1.91 | 1.00 | FAIL | FAIL |
| intensive | sand | expl | 0.98 | 1.00 | FAIL | FAIL |
| intensive | clay | impl | 10.66 | 0.95 | 51.68 | -4.08 |
| intensive | clay | expl | 18.59 | 0.84 | 64.75 | -6.98 |
| intensive | sand-loam | impl | 5.58 | 0.98 | FAIL | FAIL |
| intensive | sand-loam | expl | 8.90 | 0.96 | FAIL | FAIL |
| intensive | sand-clay | impl | 1.50 | 1.00 | 17.39 | -3.13 |
| intensive | sand-clay | expl | 9.09 | 0.96 | 15.70 | -2.37 |
| intensive | loam-clay | impl | 3.67 | 0.99 | 9.25 | 0.56 |
| intensive | loam-clay | expl | 4.45 | 0.99 | 6.74 | 0.77 |
