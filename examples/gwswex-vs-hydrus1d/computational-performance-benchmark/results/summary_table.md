# Computational performance — intensive-sand-loam

- CPU: Apple M1 Pro (10p / 10l)
- RAM: 16.0 GiB
- OS: Darwin 25.4.0 (64bit)
- OMP threads (GWSWEX): 8
- HYDRUS pool size (workers): 8
- HYDRUS setup path: phydrus.Model rebuilt per simulation

| Model | n_e | wall [s] (med) | wall [s] (min..max) | user CPU [s] | sys CPU [s] | max RSS [MiB] | disk write [MiB] |
|---|---:|---:|---|---:|---:|---:|---:|
| GWSWEX-explicit | 1 | 0.980 | 0.978..1.031 | 0.957 | 0.085 | 67.2 | 0.02 |
| GWSWEX-explicit | 10 | 2.448 | 2.429..3.663 | 10.204 | 0.156 | 68.4 | 0.53 |
| GWSWEX-explicit | 100 | 17.422 | 16.580..30.053 | 99.607 | 0.551 | 68.4 | 68.26 |
| GWSWEX-explicit | 1000 | 167.132 | 160.800..188.722 | 983.417 | 4.027 | 169.4 | 1299.14 |
| GWSWEX-implicit | 1 | 0.489 | 0.484..0.493 | 0.468 | 0.086 | 170.4 | 0.00 |
| GWSWEX-implicit | 10 | 1.228 | 1.226..1.274 | 4.431 | 0.108 | 170.4 | 0.04 |
| GWSWEX-implicit | 100 | 9.041 | 9.039..9.180 | 44.296 | 0.221 | 170.4 | 3.52 |
| GWSWEX-implicit | 1000 | 87.508 | 86.730..87.583 | 442.158 | 0.963 | 170.4 | 34.36 |
| HYDRUS-1D | 1 | 1.757 | 1.740..1.764 | 0.918 | 0.035 | 192.3 | 13.64 |
| HYDRUS-1D | 10 | 11.277 | 10.231..11.497 | 67.365 | 2.960 | 192.3 | 143.70 |
| HYDRUS-1D | 100 | 102.260 | 101.330..103.556 | 748.559 | 24.351 | 192.3 | 2223.17 |
| HYDRUS-1D | 1000 | 1054.112 | 987.558..1141.080 | 7591.133 | 260.452 | 192.3 | 18783.57 |
