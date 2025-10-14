#'Filtering species scores based on goodness of fit/frequency.
#' @description
#' The function filters species or site scores, based on either goodness of fit or absolute/relative frequency. To filter scores based on goodness of fit, it uses either [goodness()] or [envfit()] function, chosen according to the ordination type. After the filtering, th efunction rewrites the scores object in the list, created by [gordi_read()].
#' 
#' @details
#' To specify what scores to filter, use argument `what` (what = 'species' or what = 'sites'). If you decide to filter based on goodness of fit you always need to specify `slice_max`, which selects rows with the largest values of a variable (in this case goodness of fit (R²)). The [gordi_fit()] itself chooses based on ordination whether to use [goodness()] or [envfit()] to calculate goodness of fit. If the ordination was not calculated by [capscale()] and falls into one of: CCA, RDA, CA, DCA and PCA, the [goodness()] function is used to calculate goodness of fit. In this case you just need to additionally specify `choices` and `summarize` (otherwise defaults will be used). If the ordination was calculated by [capscale()], R² is automatically calculated by [envfit()]. If this is the case, you need to specify `display`, and optionally `permutations` and `choices`. There is an option to just filter scores based on frequency of occurence. This applies only to species scores and you can use either `abs_frequency` (absolute frequency) or `relat_frequency` (relative frequency). If you want to filter species scores based on frequency as well as goodness of fit, you need to specify both `slice_max` and `abs_frequency`\`relat_frequency` (first the function filters by frequency and then by goodness of fit).
#' 
#' @param pass A list object produced by [gordi_read()].
#' @param what Character; either 'species' or 'sites'.
#' @param display Character; either 'sites' for ordinary site scores and 'lc' for linear combination scores.
#' @param choices Numeric vector of selected ordination axes, defaults to the first two axis.
#' @param summarize Logical; show only the total value for each site/species, does not work though!
#' @param slice_max Numeric; select rows with the largest values of goodness of fit/R².
#' @param abs_frequency Numeric; calculates absolute frequency and selects species with higher frequency than input value.
#' @param relat_frequency Numeric; calculates relative frequency and selects species with higher frequency than input value.
#' @param permutations Number of permutations required.
#' 
#' @return The updated `pass` object with filtered species/site scores.
#' @seealso [gordi_read()], vegan :: [goodness()], vegan :: [envfit()]
#' 
#' @examples
#' # Calculate species goodness of fit when m is CCA and filter 30 best fitting species. 
#' gordi_read(m, env)|> gordi_fit(spe = spe, what = 'species', slice_max = 30, summarize = T)
#' #Calculate species R² when m is db-RDA and filter 30 best fitting species with absolute occurence more than 5.
#' 
#' @importFrom magrittr %>%
#' @importFrom dplyr mutate if_else summarize across everything mutate_all
#' @importFrom dplyr bind_cols left_join filter select slice_max case_when 
#' @importFrom tidyr pivot_longer
#' @importFrom tibble as_tibble
#' @importFrom vegan goodness envfit
#' 
#' @export

gordi_fit <- function(pass,
                      slice_max = NULL,
                      abs_frequency = NULL,
                      relat_frequency = NULL,
                      display = c('sites', 'lc'),
                      choices = 1:2,
                      summarize = T,
                      permutations = 0
){
  
  m <- pass$m

  display <- match.arg(display)
  
  if (is.null(slice_max) && is.null(abs_frequency) && is.null(relat_frequency)){
    stop("If you would like to filter based on goodness of fit/R², please specify `slice_max`. If you would like to filter based on frequency, either specify abs_frequency (absolute frequency) or relat_frequency (relative_frequency).")
  }
  
  pass$spe %>%
    mutate_all(~ if_else(. > 0, 1, 0))|> 
    summarize(across(everything(), sum))|>
    pivot_longer(everything(), names_to = 'species_names', values_to = 'abs_freq')|>
    mutate(relat_freq = abs_freq/nrow(pass$spe))-> frequency
  
  bind_cols(pass$species_names, pass$species_scores)|>
    left_join(frequency, by = 'species_names') -> spe_freq
  
  if (!is.null(abs_frequency)){  
    #select species scores based on species frequency
    spe_freq|>
      filter(abs_freq > abs_frequency) -> species_filtered
  }
  
  if (!is.null(relat_frequency)){
    #select species scores based on species frequency
    spe_freq|>
      filter(relat_freq > relat_frequency) -> species_filtered
  }
  
  type <- case_when(
    inherits(m, 'capscale') & !is.null(m$call$distance) & is.null(m$CCA)                            ~ 'PCoA', #via capscale
    inherits(m, 'capscale') & !is.null(m$call$distance)                                             ~ 'db-RDA',
    inherits(m, 'rda') & is.null(m$call$distance) & is.null(m$CCA) & !inherits(m, 'capscale')       ~ 'PCA',
    inherits(m, 'rda')                                                                              ~ 'RDA',
    inherits(m, 'cca') & is.null(m$call$distance) & is.null(m$CCA)                                  ~ 'CA',
    inherits(m, 'cca')                                                                              ~ 'CCA',
    inherits(m, 'decorana')                                                                         ~ 'DCA',
    inherits(m, c('metaMDS', 'monoMDS'))                                                            ~ 'NMDS',
    TRUE ~ paste(class(m), collapse = '/') # writes just one output
  )
  
  model <- case_when(
    type %in% c('CCA', 'RDA') ~ 'CCA',
    type %in% c('CA', 'PCA') ~ 'CA',
    TRUE ~ NA
  )
  
  spe_fitted <- NULL
  sites_fitted <- NULL
  
  if (!is.na(model)){
    warning("Using function `goodness()` to calculate goodness of fit.")
    if (!is.null(slice_max)){
    goodness_fit <- goodness(m, display = 'species', model = model, choices = choices, summarize = summarize)|>
      as_tibble()
    spe_fitted <- bind_cols(pass$species_names, goodness_fit, pass$species_scores)
    } 
  }
  
  if (inherits(m, 'capscale') || type == 'DCA'){
    warning("Using function `envfit()` to calculate goodness of fit (R²).")
    if (!is.null(slice_max)){
      goodness_fit <- envfit(m, env = pass$spe, display = display, choices = choices, permutations = permutations)
      r <- goodness_fit$vectors$r |> as_tibble()
      spe_fitted <- bind_cols(pass$species_names, r, pass$species_scores)
    }
  } 
  
  if (!is.null(slice_max) && is.null(abs_frequency) && is.null(relat_frequency)){
    warning("Filtering species based on goodness of fit (R²). If you would like to additionally filter based on frequency of occurence, set either `abs_frequency` (absolute frequency) or `relat_frequency` (relative frequency).")
    if (!is.null(spe_fitted)){
    spe_fitted <- spe_fitted|>
      slice_max(order_by = spe_fitted[2], n = slice_max)
    
    spe_fitted|>
      select(-c(1,2)) -> pass$species_scores
    spe_fitted|>
      select(1) -> pass$species_names
    } 
    if (!is.null(sites_fitted)){
      sites_fitted <- sites_fitted|>
        slice_max(order_by = sites_fitted[2], n = slice_max)
      
      sites_fitted|>
        select(-c(1,2)) -> pass$site_scores
    }
  }
  
  if (!is.null(abs_frequency) || !is.null(relat_frequency)){
    if (!is.null(species_filtered)){
      species_filtered|>
        select(c(2,3)) -> pass$species_scores
      species_filtered|>
        select(1) -> pass$species_names
    }
  }
  
  if (!is.null(slice_max) && !is.null(abs_frequency)){
    if (!is.null(spe_fitted)){
      spe_fitted|>
        left_join(frequency, by = 'species_names')|>
        filter(abs_freq > abs_frequency)|>
        slice_max(order_by = value, n = slice_max) -> spe_fitted
      
      spe_fitted|>
        select(c(3,4)) -> pass$species_scores
      spe_fitted|>
        select(1) -> pass$species_names
    }
  }
  if (!is.null(slice_max) && !is.null(relat_frequency)){
    if (!is.null(spe_fitted)){
      spe_fitted|>
        left_join(frequency, by = 'species_names')|>
        filter(relat_freq > relat_frequency)|>
        slice_max(order_by = value, n = slice_max) -> spe_fitted
      
      spe_fitted|>
        select(c(3,4)) -> pass$species_scores
      spe_fitted|>
        select(1) -> pass$species_names
    }
    }
  
  return(pass)
}