### *In Review*
**Please note that this work is currently in review at North American Journal of Fisheries Management. It is pre-decisional, and therefore should not be cited without consulting first with the author(s).**


# grasscarp

von Bertalanffy growth models for triploid Grass Carp (*Ctenopharyngodon idella*) in Lake Gaston, North Carolina and Virginia. The growth model uses a non-centered parameterization of multivariate normal VBGF parameters. One model incorporates year-specific offsets using a hierarchical prior, and the other incorporates a continuous effect of hydrilla abundance on parameters L <sub>$\infty$</sub> and K.

Models are fit using Hamiltonian Monte Carlo sampling in Stan using the No-U-Turn Sampler (NUTS). To run the analysis, the user will need a recent version of R with `rstan` installed. Note that this may require installation of other software (Rtools v 3.5) or R packages. Help for installing and getting started with `rstan` and other (R)Stan resources can be found <a href="https://mc-stan.org/users/interfaces/rstan">here</a>. 


## Files

Files are organized into subdirectories under the assumption that all workflow occurs through the use of an R project that relies on version history. Likewise, all file paths in R source files are coded as relative directories to make file read/write more streamlined.

`grasscarp_models_stan.R`: R script used for package loading, function declarations, data read and manipulation, and fitting VBGF models in Stan.This is the file to use for re-running the analysis.

`data/`: Subdirectory containing all data files necessary to conduct the analysis.

`models/`: Subdirectory containing Stan model(s) for estimating VBGF parameters.

`results/`: Output directory for results. Also contains plotting code for result of `vonbert_hydrilla_mv.stan`.
