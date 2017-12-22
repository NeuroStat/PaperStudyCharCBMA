Scripts
===========

For some set sizes (K10, K20 and K35), there is a separate folder containing a file called *master_fixran_HPC.sh*. This file contains the main body of the design of the study (figure 2 in the paper).
This is run on the UGent HPC infrastructure.
We obtain with this:
  - Group analyses of evaluation condition (N=200, 400 or 700).
  - Analyses of small studies using the 3 group level models (FE, ME and OLS) at second phase GLM
  - Fixed and random effects coordinate-based meta-analyses using *R* based on results from individual studies.

The analyses for Activation Likelihood Estimation were performed using MATLAB scripts obtained from prof. dr. Simon Eickhoff.

The other set sizes (K = 12, 14, 16, 18 and 30) are located in the ContK folder. The same scripts/procedures are used as with K = 10, 20 or 35. It is a bit more efficient though.
