# In this file we are setting up the functions for
# - the dynamics of the resources
# - the contribution of the resources to the encounter rate
# - the senescence mortality
# and provide a function for setting up Asta's model with
# benthos and algae.
#
# This file is sourced by the file run.R that runs the 
# climate change scenarios

# library(tidyverse)
# library(mizerExperimental)
# library(mizerStarvation)
library(mizer)

# see also mizerEncounter
background_encounter <- function(params, n, n_pp, n_other, ...) {
  n_pps_params <- params@other_params[["n_pps"]]
  idx_sp <- (length(params@w_full) - length(params@w) + 1):length(params@w_full)
  prey <- n_pps_params$interaction_resources %*% n_other$n_pps
  prey[, idx_sp] <- prey[, idx_sp] + params@interaction %*% n
  prey <- sweep(prey, 2, params@w_full * params@dw_full, "*")
  avail_energy <- Re(base::t(mvfft(base::t(params@ft_pred_kernel_e) * 
                                     mvfft(base::t(prey)),
                                   inverse = TRUE))) / length(params@w_full)
  avail_energy <- avail_energy[, idx_sp, drop = FALSE]
  avail_energy[avail_energy < 1e-18] <- 0
  
  params@search_vol * avail_energy
}

# see also mizerPredMort, resource_semichemostat

#### the following needs to be fixed to align with parametrisation and algorithm in resource_semichemostat
background_semichemostats <- function(params, n_other, rates, dt, component,
                                     ...) {
  c <- params@other_params[[component]]
  # name of interaction parameter for this component in species_params
  # interaction_component <- paste0("interaction_", component)
  # interaction <- params@species_params[[interaction_component]]
  mort <- base::t(c$interaction_resources) %*% rates$pred_rate
  mur <- c$rate + mort
  n_pps_steady <- c$rate * c$capacity / mur 
  return(n_pps_steady + (n_other[[component]] - n_pps_steady) * exp(-mur * dt))
}


newMultiResourceParams <- function(sp, ..., nResourceSpectra = 1, resource_sigma=0) {
  
  
  params <- newMultispeciesParams(sp, ...)
  
  #  params <- newMultispeciesParams(NS_species_params)
  
  initial_n_pps <- matrix(params@initial_n_pp/nResourceSpectra,
                          ncol = length(params@initial_n_pp),
                          nrow = nResourceSpectra )
  
  component_params <- params@resource_params
  component_params$kappa <- 
    component_params$kappa/nResourceSpectra
  
  S <- length(params@species_params$interaction_resource)
  
  component_params$interaction_resources <-
    matrix(exp(rnorm(S*nResourceSpectra,-resource_sigma^2/2,resource_sigma)),
           nrow=S, ncol=nResourceSpectra)
  
  component_params$capacity <- 
    params@cc_pp / nResourceSpectra
  
  component_params$rate <-
    params@rr_pp #/ nResourceSpectra ## ??? really need to divide?
  
  params <- setComponent(params = params, component = "n_pps",
                         initial_value = initial_n_pps,
                         dynamics_fun =  "background_semichemostats",
                         component_params = component_params)
  
  
  
  # Include these extra resources in the encounter rate
  params <- setRateFunction(params, "Encounter", "background_encounter")
  
  params
}

params <- newMultispeciesParams(NS_species_params,inter)
out1 <- project(params)


params_npps <- newMultiResourceParams(NS_species_params,inter,nResourceSpectra = 2)
out2 <- project(params_npps)

plot(out2)
