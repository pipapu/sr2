# Implementation and use example of a large number of resource size spectra in mizer
library(mizer)

###### Implementation: ######

# see also mizerEncounter
background_encounter <- function(params, n, n_pp, n_other, ...) {
  n_pps_params <- params@other_params[["n_pps"]]
  idx_sp <- (length(params@w_full) - length(params@w) + 1):length(params@w_full)
  prey <- n_pps_params$interaction %*% n_other$n_pps
  prey[, idx_sp] <- prey[, idx_sp] + params@interaction %*% n
  prey <- sweep(prey, 2, params@w_full * params@dw_full, "*")
  avail_energy <- Re(base::t(mvfft(base::t(params@ft_pred_kernel_e) * 
                                     mvfft(base::t(prey)),
                                   inverse = TRUE))) / length(params@w_full)
  avail_energy <- avail_energy[, idx_sp, drop = FALSE]
  avail_energy[avail_energy < 1e-18] <- 0
  
  encounter <- params@search_vol * avail_energy
  return(encounter)
}

### define this to extract diets of fish:
background_encounter_component <- function(params, n, n_pp, n_other, component, ...) {
  ### fix me: extract encounter with only one component of background spectrum
  # n_pps_params <- params@other_params[["n_pps"]]
  # idx_sp <- (length(params@w_full) - length(params@w) + 1):length(params@w_full)
  # prey <- n_pps_params$interaction %*% n_other$n_pps
  # prey[, idx_sp] <- prey[, idx_sp] + params@interaction %*% n
  # prey <- sweep(prey, 2, params@w_full * params@dw_full, "*")
  # avail_energy <- Re(base::t(mvfft(base::t(params@ft_pred_kernel_e) * 
  #                                    mvfft(base::t(prey)),
  #                                  inverse = TRUE))) / length(params@w_full)
  # avail_energy <- avail_energy[, idx_sp, drop = FALSE]
  # avail_energy[avail_energy < 1e-18] <- 0
  # 
  # encounter <- params@search_vol * avail_energy
  # return(encounter)
}

# see also mizerPredMort, resource_semichemostat

#### the following needs to be fixed to align with parametrisation and algorithm in resource_semichemostat
background_semichemostats <- function(params, n_other, rates, dt, component,
                                      ...) {
  c <- params@other_params[[component]]
  # name of interaction parameter for this component in species_params
  # interaction_component <- paste0("interaction_", component)
  # interaction <- params@species_params[[interaction_component]]
  mort <- base::t(c$interaction) %*% rates$pred_rate
  mur <- c$rate + mort
  n_pps_steady <- c$rate * c$capacity / mur
  new_n_pps <- n_pps_steady + (n_other[[component]] - n_pps_steady) * exp(-mur * dt)
  return(new_n_pps)
}


newMultiResourceParams <- function(sp, ..., nResourceSpectra = 1, resource_sigma=0) {
  
  
  params <- newMultispeciesParams(sp, ...)
  
  #  params <- newMultispeciesParams(NS_species_params)
  
  initial_n_pps <- matrix(params@initial_n_pp/nResourceSpectra,
                          ncol = length(params@initial_n_pp),
                          nrow = nResourceSpectra,
                          byrow = T)

  component_params <- params@resource_params
  component_params$kappa <- 
    component_params$kappa/nResourceSpectra
  
  S <- length(params@species_params$interaction_resource)
  
  component_params$interaction <-
    matrix(exp(rnorm(S*nResourceSpectra,-resource_sigma^2/2,resource_sigma)),
           nrow=S, ncol=nResourceSpectra)
  
  component_params$capacity <- 
    t(matrix(params@cc_pp / nResourceSpectra,
             ncol = nResourceSpectra,
             nrow = length(params@initial_n_pp)))
  
  component_params$rate <-
    t(matrix(params@rr_pp,
             ncol = nResourceSpectra,
             nrow = length(params@initial_n_pp)))
  
  params <- setComponent(params = params, component = "n_pps",
                         initial_value = initial_n_pps,
                         dynamics_fun =  "background_semichemostats",
                         component_params = component_params)
  
  
  
  # Include these extra resources in the encounter rate
  params <- setRateFunction(params, "Encounter", "background_encounter")
  
  params
}


###### Use examples: ######

### Compare implementation with original code (results for fish should be identical)

params <- newMultispeciesParams(NS_species_params,inter)
out1 <- project(params)
plotBiomass(out1)


nR <- 5

params_npps <- newMultiResourceParams(NS_species_params,inter,nResourceSpectra = nR)
out2 <- project(params_npps)
plotBiomass(out2)

### Now an example with a large number or resource spectra and fish species

nS <- 100  # you can make this arbitrarily large
nR <- 2*nS
nu <- 0.6  # target diet partitioning exponent
interaction_spread_sigma <- sqrt(2*log(nR+nS))/nu

initial_inter <- 
  matrix(exp(rnorm(nS*nS,-interaction_spread_sigma^2/2,interaction_spread_sigma)),
         nrow = nS, ncol = nS)

# Make all interactions < 1, fixing gamma below to compensate
inter_row_max <- apply(initial_inter, MARGIN=1, FUN = max)
initial_inter <- diag(1/inter_row_max) %*% initial_inter

# Initial set of species (no all will survivie to the end!)
initial_species_parms <- 
  data.frame(
    species=as.character(seq(nS)),
    w_inf=sort(10*exp(log(30)/log(3)*rexp(nS)))
    )
initial_species_parms$k_vb <-
  exp(0.6-(1/4)*log(initial_species_parms$w_inf))

# Generate MizerParams object:
params_npps <- 
  newMultiResourceParams(initial_species_parms,
                         initial_inter,
                         nResourceSpectra = nR, 
                         resource_sigma = interaction_spread_sigma )
with(species_params(params_npps),plot(gamma/w_inf^(0.1)~w_inf,log="xy"))

# Fix gamma as promised above:
species_params(params_npps)$gamma <-
  species_params(params_npps)$gamma * inter_row_max * 100
with(species_params(params_npps),plot(gamma/w_inf^(0.1)~w_inf,log="xy"))

# The species size spectrum
plot(log(rev(seq(nrow(species_params(params_npps)))))~log(sort(species_params(params_npps)$w_inf)))
lm(log(sort(species_params(params_npps)$w_inf))~log(rev(seq(nrow(species_params(params_npps))))))

# First simulation
out3 <- project(params_npps,t_max = 1000)

# Remove species that went extinct
species_to_remove <- tail(getBiomass(out3),n = 1) < 1e9
mean(species_to_remove)
params_npps2 <- removeSpecies(params_npps,species_to_remove)
params_npps2@other_params[["n_pps"]]$interaction <-
  params_npps@other_params[["n_pps"]]$interaction[!species_to_remove,]

# Second simulation
out4 <- project(params_npps2,t_max = 1000)
plotBiomass(out4)

# Remove species that went extinct
species_to_remove2 <- tail(getBiomass(out4),n = 1) < 1e9
mean(species_to_remove2)
params_npps3 <- removeSpecies(params_npps2,species_to_remove2)
params_npps3@other_params[["n_pps"]]$interaction <-
  params_npps2@other_params[["n_pps"]]$interaction[!species_to_remove2,]

# Third simulation
out5 <- project(params_npps3,t_max = 1000)
plotBiomass(out5)

# Remove species that went extinct
species_to_remove3 <- tail(getBiomass(out5),n = 1) < 1e9
mean(species_to_remove3)
params_npps4 <- removeSpecies(params_npps3,species_to_remove3)
params_npps4@other_params[["n_pps"]]$interaction <-
  params_npps4@other_params[["n_pps"]]$interaction[!species_to_remove3,]

# Fourth simulation (there should be no more extinctions)
out6 <- project(params_npps4,t_max = 1000)
plotBiomass(out6)

# Compare initial and final species size spectrum:
lm(log(sort(species_params(params_npps4)$w_inf))~log(rev(seq(nrow(species_params(params_npps4))))))
plot(log(rev(seq(nrow(species_params(params_npps)))))~log(sort(species_params(params_npps)$w_inf)))
points(log(rev(seq(nrow(species_params(params_npps4)))))~log(sort(species_params(params_npps4)$w_inf)),col="red")

# Some more analysis of result:
final_species_parms <- species_params(params_npps3)

BB <- getBiomass(out6)

# Transform this data in %
BB_percentage <- apply(t(BB[,-1]), 2, function(x){x*100/sum(x,na.rm=T)})

S <- nrow(BB_percentage)

library(RColorBrewer)
coul <- colorRampPalette(brewer.pal(9, "Set1"))(ncol(BB_percentage))[sample.int(S,replace = F)]

# Make a stacked barplot --> it will be in %!
barplot(BB_percentage[,seq(1,ncol(BB_percentage),length.out = 30)],col = coul, border="white", xlab="group")
# Compare this e.g. with Fig. 10 here: http://axel.rossberg.net/paper/Shephard2012ICES-JMS_Fishing_drives_LSI.pdf

## Animation of spectra takes a while, but result is interesting!  Note that, contrary to what one might expect, variability at small sizes is larger than at intermediate sizes. 
animateSpectra(out6,total = T)

# To compute diet partitioning functions (http://axel.rossberg.net/paper/Rossberg2011PRSB_Diet_Partitioning.pdf), we need to make getDiet work by registering params@other_encounter for all background components.  This requires some function closer trickery....
