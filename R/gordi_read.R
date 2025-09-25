#' Read ordination output into a standardized structure
#' 
#' @description
#' `gordi_read` is the first function in the family of `gordi_`functions. 
#' It takes the output of `rda`, `cca`, `capscale`, `decorana` and `metaMDS` and 
#' translates in into a format usable by other `gordi_` functions. `gordi_read` 
#' was not tested for partial ordination (with `Condition()`) yet, and axis labels
#' are the default ones from vegan.
#'
#'
#' @param m An ordination result: object of class `rda`, `cca`, `capscale`,
#'   `decorana`, or `metaMDS`
#' @param env A data frame with environmental variables of samples; must have
#'   the same number of samples (rows) as the species table
#' @param traits A data frame of species traits; one species can have only one
#'   trait value in its dedicated column
#' @param choices Numeric vector of selected ordination axes; defaults to the
#'   first two axes
#' @param scaling Scaling of species and site scores; see [vegan::scores.cca()]
#' @param correlation Logical; default = FALSE, see [vegan::scores.cca()] or 
#'   [vegan::scores.rda()] for details
#' @param hill Logical; default = FALSE, see [vegan::scores.cca()] for details
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
#' @seealso [gordi_sites()], [gordi_read()], [gordi_predict()], [gordi_species()]
#' 
#' @importFrom vegan capscale rda cca decorana metaMDS scores eigenvals
#' @importFrom tibble as_tibble
#' @importFrom dplyr case_when
#' 
#' @examples
#' library(vegan)
#' library(tidyverse)
#' data(dune)
#' m <- capscale(dune ~ 1)
#' o <- gordi_read(m)
#' o
#' @export
gordi_read <- function(m,
                       env = NULL,
                       traits = NULL,
                       choices = 1:2,
                       scaling = 'symm',
                       correlation = F,
                       hill = F) {

  # type of ordination 
  type <- case_when(
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
    explained_variation = if (type %in% c('DCA', 'NMDS')) {NA} else {eigenvals(m) / m$tot.chi},
    site_scores = as_tibble(as.data.frame(scores(m, display = 'sites', scaling = scaling, choices = choices, correlation = correlation, hill = hill))),
    species_scores = as_tibble(as.data.frame(scores(m, display = 'species', scaling = scaling, choices = choices, correlation = correlation, hill = hill))),
    predictor_scores = as_tibble(as.data.frame(scores(m, scaling = scaling, choices = choices, correlation = correlation, hill = hill)$biplot)),
    env = env,
    traits = traits,
    choices = choices,
    type = type,
    axis_names = colnames(scores(m, display = 'sites')),
    species_names = as_tibble(as.data.frame(scores(m, scaling = scaling, choices = choices, correlation = correlation, hill = hill)$species), rownames = 'species_names')[1],
    predictor_names = if (type %in% c('DCA', 'NMDS')) {NA} else {as_tibble(as.data.frame(scores(m, scaling = scaling, choices = choices, correlation = correlation, hill = hill)$biplot), rownames = 'predictor_names')[1]}
    )
  
  # Return pass object
  return(pass)
  
}

