# Implementation and use example of a large number of resource size spectra in mizer
library(mizer)
source("getRecruitment.R")

## but fix:
setInitialValues <- function (params, sim) 
{
  no_t <- dim(sim@n)[1]
  if (!identical(dim(sim@n)[2:3], dim(params@initial_n))) {
    stop("The consumer size spectrum of the simulation in `sim` has a ", 
         "different size from that in `params`.")
  }
  if (!identical(length(sim@n_pp[no_t, ]), length(params@initial_n_pp))) {
    stop("The resource size spectrum of the simulation in `sim` has a ", 
         "different size from that in `params`.")
  }
  if (!identical(length(sim@n_other[no_t, ]), length(params@initial_n_other))) {
    stop("The number of other components in the simulation in `sim` is ", 
         "different from that in `params`.")
  }
  if (!identical(length(sim@effort[no_t, ]), length(params@initial_effort))) {
    stop("The number of gears in the simulation in `sim` is ", 
         "different from that in `params`.")
  }
  if (!identical(dimnames(sim@effort)[[2]], names(params@initial_effort))) {
    stop("The gears in the simulation in `sim` have different names ", 
         "from those in `params`.")
  }
  params@initial_n[] <- sim@n[no_t, , ]
  params@initial_n_pp[] <- sim@n_pp[no_t, ]
  params@initial_n_other[] <- sim@n_other[no_t, ]
  params@initial_effort[] <- sim@effort[no_t, ]  ## fixed here
  params
}

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
  if(substr(component,1,4) != "n_pp") stop("This is not an n_pp component.");
  i <- as.integer(substr(component,5,nchar(component)))
  if(i <= 0) stop("Problem with component indexing.")
  n_pps_params <- params@other_params[["n_pps"]]
  idx_sp <- (length(params@w_full) - length(params@w) + 1):length(params@w_full)
  prey <- outer(n_pps_params$interaction[,i], n_other$n_pps[i,])
  prey <- sweep(prey, 2, params@w_full * params@dw_full, "*")
  avail_energy <- Re(base::t(mvfft(base::t(params@ft_pred_kernel_e) *
                                     mvfft(base::t(prey)),
                                   inverse = TRUE))) / length(params@w_full)
  avail_energy <- avail_energy[, idx_sp, drop = FALSE]
  avail_energy[avail_energy < 1e-18] <- 0

  encounter <- params@search_vol * avail_energy

  return(encounter)
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


dummy_background_dynamics <- function(...){
  return(0)
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
  
  for(i in 1:nResourceSpectra){
    params <- setComponent(params = params, component = paste0("n_pp",i),
                           initial_value = 0,
                           dynamics_fun = "dummy_background_dynamics",
                           encounter_fun =  "background_encounter_component")
  }
  
  
  # Include these extra resources in the encounter rate
  params <- setRateFunction(params, "Encounter", "background_encounter")
  
  params
}


###### Use examples: ######

### Compare implementation with original code (results for fish should be identical)

params <- newMultispeciesParams(NS_species_params,inter)
out1 <- project(params)
plotBiomass(out1)

## Generate diet partitioning function:
params <- setInitialValues(params,out1)
diets <- apply(getDiet(params,proportion = F),MARGIN = c(1,3),sum)


nR <- 2

params_npps <- newMultiResourceParams(NS_species_params,inter,nResourceSpectra = nR)
out2 <- project(params_npps)
plotBiomass(out2)

## Generate diet partitioning function:
params_npps <- setInitialValues(params_npps,out2)
diets <- apply(getDiet(params_npps,proportion = F),MARGIN = c(1,3),sum)
diets <- subset(diets, select=-c(get("Resource")))
diet_props <- as.vector(diag(1/rowSums(diets)) %*% diets)
plot(sort(diet_props/(1-diet_props),decreasing = T),seq_along(diet_props)/nrow(diets),log="xy",xlim=c(1e-2,1e2),ylim=c(0.05,20))
 
### Now an example with a large number or resource spectra and fish species

nS <- 100  # you can make this arbitrarily large
nR <- 2*nS
nu <- 0.6  # target diet partitioning exponent
interaction_spread_sigma <- sqrt(2*log(0.85*nS))/nu

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
  species_params(params_npps)$gamma * inter_row_max * 3 #10
with(species_params(params_npps),plot(gamma/w_inf^(0.1)~w_inf,log="xy"))

# The species size spectrum
plot(log(rev(seq(nrow(species_params(params_npps)))))~log(sort(species_params(params_npps)$w_inf)))
lm(log(sort(species_params(params_npps)$w_inf))~log(rev(seq(nrow(species_params(params_npps))))))

# First simulation
out3 <- project(params_npps,t_max = 100)
#params_npps_hold <- params_npps
params_npps <- setInitialValues(params_npps, out3)
plotBiomass(out3)

# Remove species that went extinct
species_to_remove <- tail(getBiomass(out3),n = 1) < 1e9
mean(species_to_remove)
params_npps2 <- removeSpecies(params_npps,species_to_remove)
params_npps2@other_params[["n_pps"]]$interaction <-
  params_npps@other_params[["n_pps"]]$interaction[!species_to_remove,]

# Second simulation
out4 <- project(params_npps2,t_max = 1000)
params_npps2 <- setInitialValues(params_npps2, out4)
plotBiomass(out4)

# Remove species that went extinct
species_to_remove2 <- tail(getBiomass(out4),n = 1) < 1e9
mean(species_to_remove2)
params_npps3 <- removeSpecies(params_npps2,species_to_remove2)
params_npps3@other_params[["n_pps"]]$interaction <-
  params_npps2@other_params[["n_pps"]]$interaction[!species_to_remove2,]

# Third simulation
out5 <- project(params_npps3,t_max = 1000)
params_npps3 <- setInitialValues(params_npps3, out5)
plotBiomass(out5)

# Remove species that went extinct
species_to_remove3 <- tail(getBiomass(out5),n = 1) < 1e9
mean(species_to_remove3)
params_npps4 <- removeSpecies(params_npps3,species_to_remove3)
params_npps4@other_params[["n_pps"]]$interaction <-
  params_npps4@other_params[["n_pps"]]$interaction[!species_to_remove3,]

# Fourth simulation (there should be no more extinctions)
out6 <- project(params_npps4,t_max = 1000)
params_npps4 <- setInitialValues(params_npps4, out6)
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
barplot(BB_percentage[,seq(1,ncol(BB_percentage),length.out = 30)],col = coul, border="white", xlab="time", ylab="Percent biomass contribution")
# Compare this e.g. with Fig. 10 here: http://axel.rossberg.net/paper/Shephard2012ICES-JMS_Fishing_drives_LSI.pdf

## Animation of spectra takes a while, but result is interesting!  Note that, contrary to what one might expect, variability at small sizes is larger than at intermediate sizes. 
#animateSpectra(out6,total = T)

# To compute diet partitioning functions (http://axel.rossberg.net/paper/Rossberg2011PRSB_Diet_Partitioning.pdf), we need to make getDiet work by registering params@other_encounter for all background components.  This requires some function closer trickery....


#params_npps4@initial_effort <- params@initial_effort   # !!! Hack to make it work
params_npps4@gear_params$knife_edge_size[] <- 0
for(target_sp in 1:S){
  paramsF <- params_npps4 #make a copy for iterative change.
  ff <- round(c(0,10^seq(log10(0.01),log10(20),length.out = 60)),3)
  rec <- rep(0,length(ff))
  ssb <- rep(0,length(ff))
  no_sp <- nrow(paramsF@species_params)
  species_params(paramsF)$species[target_sp]
  my_gear <- as.character(unique(paramsF@gear_params$gear))[1]
  for(i in seq_along(ff)){
    fishing <- ff[i]
    names(fishing) <- c(my_gear)
    cb <- matrix(as.numeric(seq_len(no_sp)==target_sp),nrow=1,ncol=no_sp,
                 dimnames = list(gear = c(my_gear), sp=row.names(species_params(paramsF))))
    paramsF <-
      setFishing(params = paramsF,
                 catchability = cb)
    outF <- project(paramsF, effort = fishing, t_max = 20)
    ssb[i] <- tail(getBiomass(outF)[,target_sp],1)  # getSSB does not work
    paramsF <- setInitialValues(paramsF,outF)
    rec[i] <- getRecruitment(paramsF)[target_sp]
    if(i > 1 && rec[i] < 0.001*max(rec)) break;
  }
  plot(rec ~ ssb,xlim = c(0,max(ssb)), ylim=c(0,max(rec)), xlab="biomass")
  feasible_fs <- which(rec > 0.01*max(rec))
  tag_fs <- 
    unique(mapply(function(f) which.min(abs(ff-f)),seq(0,max(ff[feasible_fs]),length.out = 7)))
  points(x = ssb[tag_fs], y=rec[tag_fs], col = "blue", pch = 4)
  text(x = ssb[tag_fs], y=rec[tag_fs],labels = paste0("F = ",ff[tag_fs]),pos = ifelse(ssb[tag_fs]>0.6*max(ssb),2,4), col = "blue")
  title(species_params(paramsF)$species[target_sp])
}


params <- newMultispeciesParams(NS_species_params,inter)
paramsEQ <- setInitialValues(params,project(params))
paramsEQ@gear_params$knife_edge_size[] <- 0
for(target_sp in 1:12){
  paramsF <- paramsEQ #make a copy for iterative change.
  ff <- round(c(0,10^seq(log10(0.01),log10(20),length.out = 60)),3)
  rec <- rep(0,length(ff))
  ssb <- rep(0,length(ff))
  no_sp <- nrow(paramsF@species_params)
  species_params(paramsF)$species[target_sp]
  my_gear <- as.character(unique(paramsF@gear_params$gear))[1]
  for(i in seq_along(ff)){
    fishing <- ff[i]
    names(fishing) <- c(my_gear)
    cb <- matrix(as.numeric(seq_len(no_sp)==target_sp),nrow=1,ncol=no_sp,
                 dimnames = list(gear = c(my_gear), sp=row.names(species_params(paramsF))))
    paramsF <-
      setFishing(params = paramsF,
                 catchability = cb)
    outF <- project(paramsF, effort = fishing, t_max = 10)
    ssb[i] <- tail(getBiomass(outF)[,target_sp],1)  # getSSB does not work
    paramsF <- setInitialValues(paramsF,outF)
    rec[i] <- getRecruitment(paramsF)[target_sp]
    if(i > 1 && rec[i] < 0.001*max(rec)) break;
  }
  plot(rec ~ ssb,xlim = c(0,max(ssb)), ylim=c(0,max(rec)), xlab="biomass")
  feasible_fs <- which(rec > 0.01*max(rec))
  tag_fs <- 
    unique(mapply(function(f) which.min(abs(ff-f)),seq(0,max(ff[feasible_fs]),length.out = 7)))
  points(x = ssb[tag_fs], y=rec[tag_fs], col = "blue", pch = 4)
  text(x = ssb[tag_fs], y=rec[tag_fs],labels = paste0("F = ",ff[tag_fs]),pos = ifelse(ssb[tag_fs]>0.6*max(ssb),2,4), col = "blue")
  title(species_params(paramsF)$species[target_sp])
}
