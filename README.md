# catalytic_dengue_java
Accompanying codes for the paper: Estimating dengue force of infection from age-stratified surveillance data in Java, Indonesia

Folders contents:
- data: data used for analysis (only for Jakarta province, which is open access from: https://surveilans-dinkes.jakarta.go.id/sarsbaru/dashboard.php
- codes: analysis codes:
  - stan_run: scripts to run stan models, model 1: poisson hierarchical, model 2: nb hierarchical, model 3: poisson non-hierarchical & model 4: nb non-hierarchical; province 1: jakarta & province 2: west java
- model: stan model codes
- output: model fit and figures

You can run the Jakarta data and replicate the Jakarta results by:
1. run all scripts in stan_run folder
2. run processing_output.R
3. run visualisations_script.R
