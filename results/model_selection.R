# Package load ----
library(rstan)
library(loo)

# Results read ----
load("results/vonbert_cohort-x-year.rda")
yearfit <- fit

load("results/vonbert_cohort-x-ha.rda")
hafit <- fit

load("results/vonbert_cohort-x-nha.rda")
nhafit <- fit

load("results/vonbert_cohort-x-kgha.rda")
kghafit <- fit

yearpars <- extract(yearfit)
hapars <- extract(hafit)
nhapars <- extract(nhafit)
kghapars <- extract(kghafit)

# Model selection and cross-validation -----
yearloo <- loo(yearpars$log_lik, r_eff=relative_eff(yearpars$log_lik, chain_id = rep(1:3, 5000))) 
haloo <- loo(hapars$log_lik, r_eff=relative_eff(hapars$log_lik, chain_id = rep(1:3, 5000))) 
nhaloo <- loo(nhapars$log_lik, r_eff=relative_eff(nhapars$log_lik, chain_id = rep(1:3, 5000))) 
kghaloo <- loo(kghapars$log_lik, r_eff=relative_eff(kghapars$log_lik, chain_id = rep(1:3, 5000))) 


selection_table <- loo_compare(x = list(yearloo, haloo, nhaloo, kghaloo))
selection_table

write.table(selection_table, 'results/model_selection_table.csv', row.names = FALSE,
            quote = FALSE, sep = ",")


# Statistical significance testing ----
# Year - both significant
year_linf_beta <- quantile(yearpars$ba_linf, c(0.50, 0.025, 0.975))
year_k_beta <- quantile(yearpars$ba_k, c(0.50, 0.025, 0.975))

# Hectares of hydrilla - both significant
ha_linf_beta <- quantile(hapars$ba_linf, c(0.50, 0.025, 0.975))
ha_k_beta <- quantile(hapars$ba_k, c(0.50, 0.025, 0.975))

# No. Fish per ha hydrilla - Linf significant
nha_linf_beta <- quantile(nhapars$ba_linf, c(0.50, 0.025, 0.975))
nha_k_beta <- quantile(nhapars$ba_k, c(0.50, 0.025, 0.975))

# Biomass of fish per ha hydrilla - Linf significant
kgha_linf_beta <- quantile(kghapars$ba_linf, c(0.50, 0.025, 0.975))
kgha_k_beta <- quantile(kghapars$ba_k, c(0.50, 0.025, 0.975))

# Gather coefficients for covariate effects on linf
linf_coeffs <- apply(
  rbind(year_linf_beta, ha_linf_beta, nha_linf_beta, kgha_linf_beta),
  2, round, 3
)
  
# Gather coefficients for covariate effects on k
k_coeffs <- apply(
  rbind(year_k_beta, ha_k_beta, nha_k_beta, kgha_k_beta),
  2, round, 3
  )

# Gather columns with means and 95% CRIs
linfcol <- paste0(linf_coeffs[,1], " (", linf_coeffs[,2], " - ", linf_coeffs[,3], ")")
kcol <- paste0(k_coeffs[,1], " (", k_coeffs[,2], " - ", k_coeffs[,3], ")")


coeffs_out <- data.frame(
  model = row.names(linf_coeffs),
  Linf = linfcol,
  K = kcol
)

write.table(coeffs_out, 'results/coefficient_table.csv', row.names = FALSE,
            quote = FALSE, sep = ",")






