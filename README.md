# Drug-release-modeling
Contains scrips in Python and Matlab to model: 1) In vitro drug release and 2) Cell growth in bioreactors
Designed for researchers, students and scientists in pharmaceutical and bioprocess engineering
1) kinetic_fitting.py
Fits in vitro drug release data to multiple models (Zero-order, First-order, Higuchi, Korsmeyer-Peppas), computes R2 and RMSE to assess model performance, plots original data and fitted curves
How to run the script: Python (packages: numpy, scipy, matplotlib, pandas), Define your data arrays x (time) and y (cumulative drug release)
3) cell_expansion.m
Simulates viable and total cell growth in batch bioreactor, includes glucose-limited Monod kinetics
How to run the script: Ready for ODE solvers (Matlab), can be extended to perfusion bioreactors by adjusting flow parameters, Define bioreactor volume
