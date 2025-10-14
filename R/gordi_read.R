#' Read ordination output into a standardized structure
#' 
#' @description
#' `gordi_read` is the first function in the family of `gordi_`functions. 
#' It takes the output of `rda`, `cca`, `capscale`, `decorana` and `metaMDS` and 
#' translates in into a format usable by other `gordi_` functions. `gordi_read` 
#' was not tested for partial ordination (with `Condition()`) yet, and axis labels
#' are the default ones from vegan.
#' 
#' @details
#' The function uses several arguments passed directly to \code{\link[vegan]{scores}}
#' to retrieve species, site, and predictor scores. The automatic classification
#' of the ordination type (`type` element in the returned list) is based on the
#' class and arguments used in the original \code{vegan} call.
#'
#' This function was not tested for partial ordination (with \code{Condition()}) yet,
#' and axis labels are the default ones from \code{vegan}.
#'
#'
#' @param m An ordination result: object of class \code{rda}, \code{cca}, \code{capscale},
#'   \code{decorana}, or \code{metaMDS}.
#' @param spe A data frame of species abundances or community data. This is stored
#'   in the output list but not used internally by `gordi_read`. Defaults to \code{NULL}.
#' @param env A data frame with environmental variables of samples; must have
#'   the same number of samples (rows) as the species table. Defaults to \code{NULL}.
#' @param traits A data frame of species traits; one species can have only one
#'   trait value in its dedicated column. Defaults to \code{NULL}.
#' @param choices Numeric vector of selected ordination axes; defaults to the
#'   first two axes.
#' @param scaling Scaling of species and site scores; see \code{\link[vegan]{scores.cca}}
#'   for options. Defaults to \code{'symm'}.
#' @param correlation Logical; default = \code{FALSE}. If \code{TRUE}, predictor
#'   scores are standardized to correlation scores (only relevant for linear methods).
#'   See \code{\link[vegan]{scores.cca}} or \code{\link[vegan]{scores.rda}} for details.
#' @param hill Logical; default = \code{FALSE}. If \code{TRUE}, sites and species are scaled
#'   to show separation for Hill's distance (only relevant for CCA). See \code{\link[vegan]{scores.cca}}
#'   for details.
#' @param const Numeric vector of length 1 or 2. Scaling constant(s) for sites and species
#'   when the method is not CA/CCA/RDA. Passed to \code{\link[vegan]{scores}}. Defaults to \code{c(2, 2)}. 
#' 
#' 
#' @return A list with elements:
#'   \describe{
#'     \item{m}{the original ordination object}
#'     \item{explained_variation}{eigenvalues divided by total inertia, or NA for DCA/NMDS}
#'     \item{site_scores}{tibble of site (sample) scores}
#'     \item{species_scores}{tibble of species scores}
#'     \item{predictor_scores}{tibble of biplot scores for environmental predictors}
#'     \item{env}{environmental data frame (if supplied)}
#'     \item{traits}{trait data frame (if supplied)}
#'     \item{choices}{selected axes}
#'     \item{type}{character string with ordination type}
#'     \item{axis_names}{axis labels}
#'     \item{species_names}{species names}
#'     \item{predictor_names}{predictor names, or NA for DCA/NMDS}
#'   }
#' 
#' @seealso \code{\link[gordi]{gordi_sites}}, \code{\link[gordi]{gordi_predict}}, \code{\link[gordi]{gordi_species}}
#' 
#' @importFrom vegan capscale rda cca decorana metaMDS scores eigenvals
#' @importFrom tibble as_tibble
#' @importFrom dplyr case_when filter
#' 
#' @examples
#' \dontrun{
#' library(vegan)
#' library(tidyverse)
#' data(dune)
#' data(dune.env)
#'
#' # Example 1: db-RDA
#' m_dbRDA <- capscale(dune ~ Management + A1, data = dune.env, distance = "bray")
#' o_dbRDA <- gordi_read(m_dbRDA, env = dune.env)
#' print(o_dbRDA)
#'
#' # Example 2: Simple PCoA (capscale with no predictors)
#' m_pcoa <- capscale(dune ~ 1, distance = "bray")
#' o_pcoa <- gordi_read(m_pcoa)
#' print(o_pcoa$type) # Should be "PCoA"
#' }
#' @export
gordi_read <- function(m,
                       spe = NULL,
                       env = NULL,
                       traits = NULL,
                       choices = 1:2,
                       scaling = 'symm',
                       correlation = F,
                       hill = F,
                       const = c(2,2)) {

  # type of ordination 
  type <- dplyr::case_when(
    inherits(m, 'capscale') & !is.null(m$call$distance) & is.null(m$CCA) ~ 'PCoA', #via capscale
    inherits(m, 'capscale') & !is.null(m$call$distance)                  ~ 'db-RDA',
    inherits(m, 'rda') & is.null(m$call$distance) & is.null(m$CCA)       ~ 'PCA',
    inherits(m, 'rda')                                                   ~ 'RDA',
    inherits(m, 'cca') & is.null(m$call$distance) & is.null(m$CCA)       ~ 'CA',
    inherits(m, 'cca')                                                   ~ 'CCA',
    inherits(m, 'decorana')                                              ~ 'DCA',
    inherits(m, c('metaMDS', 'monoMDS'))                                 ~ 'NMDS',
    TRUE ~ paste(class(m), collapse = '/') # writes just one output
  )

  # create pass object 
  pass <- list(
    m = m,
    explained_variation = if (type %in% c('DCA', 'NMDS')) {NA} else {vegan::eigenvals(m) / m$tot.chi},
    site_scores = tibble::as_tibble(as.data.frame(vegan::scores(m, display = 'sites', scaling = scaling, choices = choices, correlation = correlation, hill = hill, const = const))),
    species_scores = tibble::as_tibble(as.data.frame(vegan::scores(m, display = 'species', scaling = scaling, choices = choices, correlation = correlation, hill = hill, const = const))),
    predictor_scores = if (type %in% c('DCA', 'NMDS')) {NA} else {vegan::scores(m, choices = choices, display = 'all', scaling = scaling, tidy = T) |> tibble::as_tibble() |> dplyr::filter(score %in% c('biplot', 'centroids'))},
    env = env,
    traits = traits,
    choices = choices,
    type = type,
    axis_names = colnames(vegan::scores(m, display = 'sites', choices = choices)),
    species_names = tibble::as_tibble(as.data.frame(vegan::scores(m, scaling = scaling, choices = choices, correlation = correlation, hill = hill)$species), rownames = 'species_names')[1],
    predictor_names = if (type %in% c('DCA', 'NMDS')) {NA} else {tibble::as_tibble(as.data.frame(vegan::scores(m, scaling = scaling, choices = choices, correlation = correlation, hill = hill)$biplot), rownames = 'predictor_names')[1]},
    spe = spe
    )
  
  # Return pass object
  return(pass)
  
}

