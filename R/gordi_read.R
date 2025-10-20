#' Read ordination output into a standardized structure
#' 
#' @description
#' \code{gordi_read()} is the first function in the family of \code{gordi_} functions. 
#' It takes the output of \code{rda}, \code{cca}, \code{capscale}, \code{decorana} and \code{metaMDS} and 
#' translates in into a format usable by other \code{gordi_} functions. It can also read ordination objects
#' which were calculated on dissimilarity matrices.
#' 
#' @details
#' The function uses several arguments passed directly to \code{\link[vegan]{scores}}
#' to retrieve species, sites and predictors scores. The automatic classification
#' of the ordination type (`type` element in the returned list) is based on the
#' class and arguments used in the original \code{vegan} call.
#' 
#' \strong{Handling Distance-Based Ordinations (PCoA, db-RDA, NMDS):}
#' For ordinations calculated directly on a distance matrix (like PCoA or NMDS), 
#' species scores cannot be calculated internally by \code{vegan}. But you supply
#' the \strong{original species data} to the \code{spe} argument, \code{gordi_read}
#' will calculate species scores using an appropriate method (\code{envfit} for
#' LC scores, or \code{wascores} for WA scores) based on the \code{species_score_type} argument.
#'
#' This function was not tested for partial ordination (with \code{Condition()}) yet,
#' and axis labels are the default ones from \code{vegan}.
#'
#'
#' @param m An ordination result: object of class \code{rda}, \code{cca}, \code{capscale},
#'   \code{decorana}, or \code{metaMDS}.
#' @param spe A data frame of species abundances or community data. This is crucial
#'   if the ordination (\code{capscale} or \code{metaMDS}) was calculated on a 
#'   distance matrix and you require species scores. Defaults to \code{NULL}.
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
#' @param species_score_type Character string specifying the type of species scores 
#'   to calculate for distance-based ordinations (db-RDA/PCoA/NMDS) when species 
#'   data is supplied via \code{spe}. Must be one of \code{'lc'} (linear combination, default) or 
#'   \code{'wa'} (weighted average). LC scores are calculated via \code{envfit} and 
#'   WA scores via \code{wascores}. Note: WA scores are not supported for PCoA yet.
#' 
#' 
#' @return A list with elements:
#'  \describe{
#'    \item{m}{The original ordination object.}
#'    \item{type}{Character string with the determined ordination type (e.g., 'RDA', 'db-RDA', 'NMDS').}
#'    \item{explained_variation}{A numeric vector of eigenvalues divided by total inertia. Returns \code{NA} for DCA and NMDS.}
#'    \item{choices}{The numeric vector of selected ordination axes (e.g., \code{c(1, 2)}).}
#'    \item{scaling}{Character string specifying the scaling used for scores (e.g., 'symm').}
#'    \item{correlation}{Logical; value of the \code{correlation} argument used for scoring.}
#'    \item{hill}{Logical; value of the \code{hill} argument used for scoring.}
#'    \item{const}{Numeric vector of scaling constant(s) for scores.}
#'    \item{axis_names}{Character vector of axis labels (e.g., \code{c('RDA1', 'RDA2')}).}
#'    \item{site_scores}{Tibble of site (sample) scores (columns: Axis1, Axis2, ...).}
#'    \item{species_scores}{Tibble of species scores (columns: Axis1, Axis2, ...). Can be \code{NULL} if scores could not be calculated.}
#'    \item{predictor_scores}{Tibble of biplot/centroid scores for environmental predictors. Can be \code{NULL}.}
#'    \item{species_names}{Tibble containing a single column, \code{species_names}, with the names of all species, matching the scores.}
#'    \item{predictor_names}{Tibble containing a single column, \code{predictor_names}, with the names of all predictors/factors. Can be \code{NULL}.}
#'    \item{env}{The environmental data frame (\code{env}) supplied to the function. Can be \code{NULL}.}
#'    \item{spe}{The species data frame (\code{spe}) supplied to the function. Can be \code{NULL}.}
#'    \item{traits}{The trait data frame (\code{traits}) supplied to the function. Can be \code{NULL}.}
#'  }
#' 
#' @seealso \code{\link[gordi]{gordi_sites}}, \code{\link[gordi]{gordi_predict}}, \code{\link[gordi]{gordi_species}}
#' 
#' @importFrom vegan capscale rda cca decorana metaMDS scores eigenvals wascores envfit
#' @importFrom tibble as_tibble
#' @importFrom dplyr case_when filter select all_of
#' @importFrom rlang .data
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
                       env = NULL,
                       spe = NULL,
                       traits = NULL,
                       choices = 1:2,
                       scaling = 'symm',
                       correlation = F,
                       hill = F,
                       const = c(2,2),
                       species_score_type = c('lc', 'wa')) {
  
  # type of ordination 
  type <- dplyr::case_when(
    inherits(m, 'capscale') && is.null(m$CCA)                            ~ 'PCoA', #via capscale
    inherits(m, 'capscale') && !is.null(m$CCA)                           ~ 'db-RDA',
    inherits(m, 'rda') && !inherits(m, 'capscale') && is.null(m$CCA)     ~ 'PCA',
    inherits(m, 'rda') && !inherits(m, 'capscale')                       ~ 'RDA',
    inherits(m, 'cca') && is.null(m$call$distance) && is.null(m$CCA)     ~ 'CA',
    inherits(m, 'cca')                                                   ~ 'CCA',
    inherits(m, 'decorana')                                              ~ 'DCA',
    inherits(m, c('metaMDS', 'monoMDS'))                                 ~ 'NMDS',
    TRUE ~ paste(class(m), collapse = '/') # writes just one output
  )
  
  
  # force 'lc' as default when NA is provided
  species_score_type <- match.arg(species_score_type)
  
  # --- create empty pass object ---
  pass <- NULL
  
  # --- pre-extract predictor scores ---
  # Only constrained ordinations (RDA, CCA, db-RDA) have predictor scores.
  
  
  # Warning when interactions are present, but env table is not provided
  if (type %in% c('RDA', 'CCA', 'db-RDA') && 
      scores(m,
             display = 'all', 
             choices = choices, 
             scaling = scaling, 
             correlation = correlation, 
             hill = hill, 
             const = const, 
             tidy = T) |> 
      as_tibble() |> 
      filter(str_detect(label, ':|*')) |> 
      nrow() > 0 && is.null(env)) {
    stop('The model contains interaction terms, but for the correct calculation of their scores, the `env` data frame must be supplied to the `gordi_read()` function.')
  }
  
  
  pred_scores <- if (type %in% c('RDA', 'db-RDA')) {
    scores(m, 
           display = 'all', 
           choices = choices, 
           scaling = scaling, 
           correlation = correlation, 
           hill = hill, 
           const = const, 
           tidy = T) |> 
      as_tibble() |> 
      dplyr::filter(score %in% c('biplot', 'centroids'))
  } else if (type == 'CCA') {
    scores(m, 
           display = 'all', 
           choices = choices, 
           scaling = scaling, 
           correlation = correlation, 
           hill = hill, 
           const = const, 
           tidy = T) |> 
      as_tibble() |> 
      dplyr::filter(score %in% c('biplot', 'centroids')) |> 
      dplyr::select(-weight)
  } else {
    NULL
  }
  
  
  # --- if else block starts ---
  # 1. for PCA, RDA, CA and CCA - create pass object without hesitation
  if (type %in% c('PCA', 'RDA', 'CA', 'CCA', 'DCA')) {
    
    # DCA specific check
    if (type == 'DCA' && isTRUE(all.equal(rownames(scores(m, display = "sites")), rownames(scores(m, display = "species"))))) {
      stop('It seems like the DCA was calculated on a distance matrix. Even though decorana() allows you to do that, it is nonsense. DCA belongs to a group of ordination methods calculating with raw species data. Please calculate the DCA on species data table.')
    }
    
    # species and sites score extraction
    si_scores <- scores(m, scaling = scaling, choices = choices, correlation = correlation, hill = hill, const = const, tidy = T) |> as_tibble() |> dplyr::filter(score == 'sites') |> select(-c(score))
    spe_scores <- scores(m, scaling = scaling, choices = choices, correlation = correlation, hill = hill, const = const, tidy = T) |> as_tibble() |> dplyr::filter(score == 'species') |> select(-c(score))
    
    pass <- list(m = m,
                 type = type,
                 explained_variation = if (type %in% c('DCA', 'NMDS')) {NA} else {vegan::eigenvals(m) / m$tot.chi},
                 choices = choices,
                 scaling = scaling, 
                 correlation = correlation, 
                 hill = hill, 
                 const = const,
                 axis_names = colnames(scores(m, display = 'sites', choices = choices)),
                 site_scores = if(type %in% c('CA', 'CCA', 'DCA')) {si_scores |> dplyr::select(-c(weight, label))} else {si_scores |> dplyr::select(-label)},
                 species_scores = if(type %in% c('CA', 'CCA', 'DCA')) {spe_scores |> dplyr::select(-c(weight, label))} else {spe_scores |> dplyr::select(-label)},
                 predictor_scores = pred_scores,
                 species_names = spe_scores |> dplyr::select(species_names = .data$label),
                 predictor_names = if (!is.null(pred_scores)) {pred_scores |> dplyr::select(predictor_names = label)} else {NULL},
                 env = env,
                 spe = spe,
                 traits = traits)
    
    # 2. PCoA, db-RDA                 
  } else if (type %in% c('PCoA', 'db-RDA')) {
    
    # Core Site Scores (Required in all branches)
    si_scores <- scores(m, scaling = scaling, choices = choices, correlation = correlation, hill = hill, const = const, tidy = T) |> as_tibble() |> dplyr::filter(score == 'sites') |> select(-c(score))
    
    # 3.1 If PCoA or db-RDA were calculated on species table, continue and create pass object.
    if (isFALSE(all(is.na(vegan::scores(m, display = "species"))))) {
      
      # Species scores extraction
      spe_scores <- scores(m, scaling = scaling, choices = choices, correlation = correlation, hill = hill, const = const, tidy = T) |> as_tibble() |> dplyr::filter(score == 'species') |> select(-c(score))
      
      pass <- list(
        m = m,
        type = type,
        explained_variation = if (type %in% c('DCA', 'NMDS')) {NA} else {vegan::eigenvals(m) / m$tot.chi},
        choices = choices,
        scaling = scaling, 
        correlation = correlation, 
        hill = hill, 
        const = const,
        axis_names = colnames(scores(m, display = 'sites', choices = choices)),
        site_scores = si_scores |> dplyr::select(-label),
        species_scores = spe_scores |> select(-label),
        predictor_scores = pred_scores,
        species_names = spe_scores |> select(species_names = .data$label),
        predictor_names = if (!is.null(pred_scores)) {pred_scores |> dplyr::select(predictor_names = label)} else {NULL},
        env = env,
        spe = spe,
        traits = traits
      )
      
      # 3.2 If PCoA or db-RDA were calculated on distance matrix, check is spe table was supplied.       
    } else if (isTRUE(all(is.na(scores(m, display = "species"))))) {
      
      # 3.2.1 If spe table was not supplied, warn the user that pass will include just site scores, but species scores will be empty.
      if (is.null(spe)) {
        
        warning(paste0('The ', type, ' was calculated on distance matrix, but you didn\'t supply the species data into the argument `spe`. Species scores were not calculated, but only site scores are available.'))
        
        pass <- list(
          m = m,
          type = type,
          explained_variation = if (type %in% c('DCA', 'NMDS')) {NA} else {vegan::eigenvals(m) / m$tot.chi},
          choices = choices,
          scaling = scaling, 
          correlation = correlation, 
          hill = hill, 
          const = const,
          axis_names = colnames(scores(m, display = 'sites', choices = choices)),
          site_scores = si_scores |> dplyr::select(-label),
          species_scores = NULL,
          predictor_scores = pred_scores,
          species_names = NULL,
          predictor_names = if (!is.null(pred_scores)) {pred_scores |> dplyr::select(predictor_names = label)} else {NULL},
          env = env,
          spe = spe,
          traits = traits
        )
        
        # 3.2.2 If spe table was supplied, check, if it is db-RDA or PCoA    
      } else {
        
        # 3.2.2.1 For db-RDA, there are two options - lc and wa species scores. LC are the default.      
        if(type == 'db-RDA') {
          
          # 3.2.2.1.1 If species_score_type was defined as lc (or was not defined) calculate lc scores and create pass object.                    
          if(species_score_type == 'lc') {
            
            dbrda_spe_scores_lc <- envfit(m, env = spe, display = 'lc', choices = choices) |> 
              scores(display = 'vectors') |> 
              as_tibble(rownames = 'species_names') 
            
            pass <- list(
              m = m,
              type = type,
              explained_variation = if (type %in% c('DCA', 'NMDS')) {NA} else {vegan::eigenvals(m) / m$tot.chi},
              choices = choices,
              scaling = scaling, 
              correlation = correlation, 
              hill = hill, 
              const = const,
              axis_names = colnames(scores(m, display = 'sites', choices = choices)),
              site_scores = si_scores |> dplyr::select(-label),
              species_scores = dbrda_spe_scores_lc |> select(-species_names),
              predictor_scores = pred_scores,
              species_names =  dbrda_spe_scores_lc |> select(species_names),
              predictor_names = if (!is.null(pred_scores)) {pred_scores |> dplyr::select(predictor_names = label)} else {NULL},
              env = env,
              spe = spe,
              traits = traits
            )
            
            warning(paste0('For this ',  type, ' `lc` scores of species were calculated.'))
            
            # 3.2.2.1.2 If species_score_type was defined as wa calculate wa scores and create pass object.                
          } else if (species_score_type == 'wa') {
            
            dbrda_spe_scores_wa <- wascores(scores(m, display = 'lc', choices = choices), spe) |> 
              as_tibble(rownames = 'species_names') 
            
            pass <- list(
              m = m,
              type = type,
              explained_variation = if (type %in% c('DCA', 'NMDS')) {NA} else {vegan::eigenvals(m) / m$tot.chi},
              choices = choices,
              scaling = scaling, 
              correlation = correlation, 
              hill = hill, 
              const = const,
              axis_names = colnames(scores(m, display = 'sites', choices = choices)),
              site_scores = si_scores |> dplyr::select(-label),
              species_scores = dbrda_spe_scores_wa |> select(-species_names),
              predictor_scores = pred_scores,
              species_names = dbrda_spe_scores_wa |> select(species_names),
              predictor_names = if (!is.null(pred_scores)) {pred_scores |> dplyr::select(predictor_names = label)} else {NULL},
              env = env,
              spe = spe,
              traits = traits
            )
            
            warning(paste0('For this ',  type, ' `wa` scores of species were calculated.'))
            
          }
          
          # 3.2.2.2 For PCoA, there are theoretically also two options - lc and wa, but currently, only lc scores will be calculated.          
        } else if(type == 'PCoA') {
          
          # 3.2.2.2.1 If species_score_type was defined as lc (or was not defined) calculate lc scores and create pass object.                  
          if(species_score_type == 'lc') {
            
            pcoa_spe_scores_lc <- envfit(m, env = spe, display = 'sites', choices = choices) |> 
              scores(display = 'vectors') |> 
              as_tibble(rownames = 'species_names') 
            
            pass <- list(
              m = m,
              type = type,
              explained_variation = if (type %in% c('DCA', 'NMDS')) {NA} else {vegan::eigenvals(m) / m$tot.chi},
              choices = choices,
              scaling = scaling, 
              correlation = correlation, 
              hill = hill, 
              const = const,
              axis_names = colnames(scores(m, display = 'sites', choices = choices)),
              site_scores = si_scores |> dplyr::select(-label),
              species_scores = pcoa_spe_scores_lc |> select(-species_names),
              predictor_scores = pred_scores,
              species_names = pcoa_spe_scores_lc |> select(species_names),
              predictor_names = if (!is.null(pred_scores)) {pred_scores |> dplyr::select(predictor_names = label)} else {NULL},
              env = env,
              spe = spe,
              traits = traits
            )
            
            warning(paste0('For this ',  type, ' `lc` scores of species were calculated.'))
            
            # 3.2.2.2.1 If species_score_type was defined as wa, don't calculate species scores and create pass object without them.            
          } else if (species_score_type == 'wa') {
            
            warning('WA (weighted average) scores for species are not available for PCoA calculated on distance matrix yet. The discussion about their application is still going on, so they may be available in the future. At this moment, species scores were not calculated.')  
            
            pass <- list(
              m = m,
              type = type,
              explained_variation = if (type %in% c('DCA', 'NMDS')) {NA} else {vegan::eigenvals(m) / m$tot.chi},
              choices = choices,
              scaling = scaling, 
              correlation = correlation, 
              hill = hill, 
              const = const,
              axis_names = colnames(scores(m, display = 'sites', choices = choices)),
              site_scores = si_scores |> dplyr::select(-label),
              species_scores = NULL,
              predictor_scores = pred_scores,
              species_names = NULL,
              predictor_names = if (!is.null(pred_scores)) {pred_scores |> dplyr::select(predictor_names = label)} else {NULL},
              env = env,
              spe = spe,
              traits = traits
            )
            
          }
        }
      }
    }
    
    # 4. NMDS  
  } else if (type == 'NMDS') {
    
    # core site scores (required in all branches)
    si_scores <- scores(m, display = 'sites', scaling = scaling, choices = choices, correlation = correlation, hill = hill, const = const, tidy = T) |> as_tibble() |> dplyr::filter(score == 'sites') |> dplyr::select(-score)
    
    # 4.1 If species scores exist in m, create pass object
    if(!is.null(scores(m, display = 'species', tidy = T))) {
      
      spe_scores <- scores(m, scaling = scaling, choices = choices, correlation = correlation, hill = hill, const = const, tidy = T) |> as_tibble() |> dplyr::filter(score == 'species') |> dplyr::select(-score)
      
      pass <- list(
        m = m,
        type = type,
        explained_variation = if (type %in% c('DCA', 'NMDS')) {NA} else {vegan::eigenvals(m) / m$tot.chi},
        choices = choices,
        scaling = scaling, 
        correlation = correlation, 
        hill = hill, 
        const = const,
        axis_names = colnames(scores(m, display = 'sites', choices = choices)),
        site_scores = si_scores |> dplyr::select(-label),
        species_scores = spe_scores |> select(-label),
        predictor_scores = pred_scores,
        species_names = spe_scores |> select(species_names = .data$label),
        predictor_names = if (!is.null(pred_scores)) {pred_scores |> dplyr::select(predictor_names = label)} else {NULL},
        env = env,
        spe = spe,
        traits = traits
      ) 
      
      # 4.2 If species scores does not exist, calculate wa scores and create pass object.      
    } else if (is.null(scores(m, display = 'species', tidy = T))) {
      
      # 4.2.1 If spe table was not supplied, warn the user that pass will include just site scores and species scores will be empty. 
      if(is.null(spe)) {
        
        warning(paste0('The ', type, ' was calculated on distance matrix, but you didn\'t supply the species data into the argument `spe`. Species scores were not calculated, but site scores are available.'))
        
        pass <- list(
          m = m,
          type = type,
          explained_variation = if (type %in% c('DCA', 'NMDS')) {NA} else {vegan::eigenvals(m) / m$tot.chi},
          choices = choices,
          scaling = scaling, 
          correlation = correlation, 
          hill = hill, 
          const = const,
          axis_names = colnames(vegan::scores(m, display = 'sites', choices = choices)),
          site_scores = si_scores |> dplyr::select(-label),
          species_scores = NULL,
          predictor_scores = pred_scores,
          species_names = NULL,
          predictor_names = if (!is.null(pred_scores)) {pred_scores |> dplyr::select(predictor_names = label)} else {NULL},
          env = env,
          spe = spe,
          traits = traits
        ) 
      } else
        
        # 4.3. If spe table was supplied, calculate wa species scores      
        if (!is.null(spe)) {
          
          nmds_spe_scores_wa <- wascores(scores(m, display = 'si', choices = choices), spe, expand = T) |> 
            as.data.frame() |> 
            as_tibble(rownames = 'species_names')
          
          pass <- list(
            m = m,
            type = type,
            explained_variation = if (type %in% c('DCA', 'NMDS')) {NA} else {vegan::eigenvals(m) / m$tot.chi},
            choices = choices,
            scaling = scaling, 
            correlation = correlation, 
            hill = hill, 
            const = const,
            axis_names = colnames(vegan::scores(m, display = 'sites', choices = choices)),
            site_scores = si_scores |> dplyr::select(-label),
            species_scores = nmds_spe_scores_wa |> select(-species_names),
            predictor_scores = pred_scores,
            species_names = nmds_spe_scores_wa |> select(species_names),
            predictor_names = if (!is.null(pred_scores)) {pred_scores |> dplyr::select(predictor_names = label)} else {NULL},
            env = env,
            spe = spe,
            traits = traits
          )  
          
          warning(paste0('For this ',  type, ' `wa` scores of species were calculated. `lc` are not an option for this type of ordination.'))
          
        }
    }
  }
  
  # Return pass object
  return(pass)
  
}