This folder contains the R code for the paper:
“Joint Modeling of Semi-continuous Longitudinal Exposures and a Survival Outcome”

Contents

`share.func.R`: Shared helper functions used by the simulation scripts.

`sim-GVA-multi-2.R`: Main simulation script. Running this script reproduces the simulation results reported in Table 1 of the paper, and Tables C1-C4 in the supplementary materials.

`sim-GVA-multi-5.R`: Main simulation script. Running this script reproduces the simulation results reported in Table 2 of the paper.

`sim-GVA-multi-10.R`: Main simulation script. Running this script reproduces the simulation results reported in Table C5 in the supplementary materials.

How to run

1. Open R or RStudio and set the working directory to this folder.
2. Run the simulation: Create a folder named `results`. Source or run `sim-GVA-multi-2.R` to reproduce Table 1.

Notes

The scripts assume that `share.func.R` is available in the same directory (or sourced with the correct relative path).
Output files (if any) will be written to the working directory unless otherwise specified within the scripts.

Contact

For questions about the code, please check the github page: https://github.com/k285li/joint-model, or contact the corresponding author of the manuscript.
