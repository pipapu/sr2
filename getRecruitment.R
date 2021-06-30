#' Calculate the recruitment rate for all species
#' 
#' Calculates the recruitment rate for the species in
#' the `MizerParams` class. The recruitment rate is defined as the rate
#' at which biomass enters maturity.
#' 
#' @references Rossberg, A. G., Houle, J. E. & Hyder, K. 
#'   Stock–recruitment relations controlled by feeding interactions alone. 
#'   Can. J. Fish. Aquat. Sci. 70, 1447–1455 (2013)
#' 
#' @param sim An object of class `MizerParams`.
#' @param n	 matrix of species abundances (species x size).
#' @param n_pp A vector of the resource abundance by size
#'   ecosystem
#' @param t	The time for which to do the calculation (Not used by standard 
#'   mizer rate functions but useful for extensions with time-dependent 
#'   parameters.)
#' @param ... Unused
#'   
#' @return A named vector containing the recruitment rate for each species in
#'   grams per year.
#' @export
#' @examples
#' \dontrun{
#' getRecruitment(NS_params)
#' }
getRecruitment <- function(params, 
                           n = initialN(params),
                           n_pp = initialNResource(params),
                           n_other = initialNOther(params),
                           t = 0,
                           ...) {
  #assert_that(is(params, "MizerParams"))
  no_sp <- nrow(params@species_params)
  # Here we approximate the derivative of the maturity ogive numerically
  diff_maturity <- t(rbind(diff(t(params@maturity)),
                           rep(0, no_sp)))
  growth <- getEGrowth(params, n = n, n_pp = n_pp, n_other = n_other, 
                       t = t, ...)
  # Equation 4 from Axel's paper
  recruitment <- (growth * n * diff_maturity) %*% params@w
  return(recruitment[, 1])
}
