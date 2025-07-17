# flexIC Simulation Code and Results

This repository contains the simulation materials associated with the manuscript:

> **flexIC: A Tunable Extension of the Iman–Conover Procedure for Marginal-Preserving Correlated Simulation in R**  
> Submitted to *Behavior Research Methods* (BR-Org-25-581)

---

## Overview

The `flexIC` method extends the classic Iman–Conover rank-correlation algorithm with tunable precision (`ε`) and deterministic marginal preservation. This repository includes:

- R code for running the full simulation
- Output table comparing `flexIC` to Iman–Conover across conditions
- Runtime and bias reduction summaries for each cell

All simulations use only base R and the `MASS` package (for `mvrnorm()`).

---

## Contents

| File | Description |
|------|-------------|
| `FINALGOOD1.R` | Complete script to reproduce the simulation |
| `flexIC_vs_IC_bias_results.csv` | Summary output of bias, runtime, and redraws |
| `LICENSE.md` | MIT License |
| `README.md` | This file |

---

## Reproducibility

To rerun the simulation:

```r
source("FINALGOOD1.R")
