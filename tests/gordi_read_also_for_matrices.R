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
    si_scores <- scores(m, scaling = scaling, choices = choices, correlation = correlation, hill = hill, const = const, tidy = T) |> as_tibble() |> dplyr::filter(score == 'sites') |> dplyr::select(-score)
    
  # 4.1 If species scores exist in m, create pass object
    if(!is.null(scores(m, display = 'species', tidy = T))) {
    
      spe_scores <- scores(m, scaling = scaling, choices = choices, correlation = correlation, hill = hill, const = const, tidy = T) |> as_tibble() |> dplyr::filter(score == 'species') |> dplyr::select(-score)

      pass <- list(
        m = m,
        type = type,
        explained_variation = if (type %in% c('DCA', 'NMDS')) {NA} else {vegan::eigenvals(m) / m$tot.chi},
        choices = choices,
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




# Testing -----------------------------------------------------------------


data(dune)
data(dune.env)

dune.dist <- vegdist(dune, method = 'eucli')

m <- rda(dune ~ 1)
m <- rda(dune ~ A1, dune.env)
m <- cca(dune ~ 1, dune.env)
m <- cca(dune ~ A1, dune.env)
m <- decorana(dune)
m <- decorana(dune.dist)
m <- capscale(dune ~ 1)
m <- capscale(dune.dist ~ 1)
m <- capscale(dune ~ A1, dune.env)
m <- capscale(dune.dist ~ A1, dune.env)
m <- metaMDS(dune, k = 3)
m <- metaMDS(dune.dist, k = 3)


gordi_read(m)
gordi_read2(m, spe = dune, env = as_tibble(dune.env, rownames = 'ID'), species_score_type = 'lc') |> 
  gordi_sites() |> 
  gordi_species() |> 
  gordi_corr(variables = c('A1', 'Use', 'Management'), permutations = 999, p_val_adjust = T, p_val_adjust_method = 'holm', colour = 'name') |> 
  gordi_label(what = 'sites', repel_label = T)

scores(m, scaling = 'species', choices = 1:4, correlation = T, hill = F, const = c(1,1), tidy = T) |> as_tibble()



scores(m, tidy = T)

as_tibble(as.data.frame(scores(m, display = 'species', scaling = 'species', choices = 1:2,
                               correlation = T, hill = T, const = c(1,2))))

envfit(m, env = dune, display = 'sites', choices = 1:2, correlation = T, hill = T) |> 
  scores(display = 'vectors')



tibble::as_tibble(scores(m, scaling = 'sym', display = 'species',choices = 1:2, correlation = T, hill = F), rownames = 'species_names') |> select(species_names)

inherits(m, 'capscale') & !is.null(m$call$distance) & is.null(m$CCA)

class(m)

as_tibble(NULL)


m <- capscale(dune.dist ~ A1 + Use, dune.env)
m <- capscale(dune.dist ~ 1)
m <- rda(dune.dist ~ 1) # nope
m <- cca(dune.dist ~ 1) # nope
m <- decorana(dune.dist) # je to blbost, nema se to tak vubec pocitat
m <- decorana(dune)
m1 <- metaMDS(dune.dist, k = 3)
m2 <- metaMDS(dune, k = 3)

# db-RDA, weighted averages (min preferovana moznost)
wascores(scores(m, display = 'lc'), dune)

wascores(scores(m, display = 'lc'), dune) |> 
  as_tibble(rownames = 'species_names') |> 
  select(-species_names)

wascores(scores(m, display = 'lc'), dune) |> 
  as_tibble(rownames = 'species_names') |> 
  select(species_names)

# db-RDA, 'lc' scores = korelace druhu s LC scores samplu
envfit(m, env = dune, display = 'lc') 

envfit(m, env = dune, display = 'lc') |> 
  scores(display = 'vectors') |> 
  as_tibble(rownames = 'species_names') |> 
  select(-species_names)

envfit(m, env = dune, display = 'lc') |> 
  scores(display = 'vectors') |> 
  as_tibble(rownames = 'species_names') |> 
  select(species_names)



# PCoA wa zatim nedelame, ale mozna to neni blbost
#...
# PCoA, 'lc' = linearni korelace, sipky ('lc')
envfit(m, env = dune, display = 'sites')

envfit(m, env = dune, display = 'sites') |> 
  scores(display = 'vectors') |> 
  as_tibble(rownames = 'species_names') |> 
  select(-species_names)

envfit(m, env = dune, display = 'sites') |> 
  scores(display = 'vectors') |> 
  as_tibble(rownames = 'species_names') |> 
  select(species_names)

# NMDS - wa scores = druhova optima ('lc' nemaji smysl)
wa_m1 <- wascores(scores(m1, display = 'si'), dune, expand = T)

wascores(scores(m1, display = 'si', choices = 1:2), dune, expand = T) |> 
  as.data.frame() |> 
  as_tibble(rownames = 'species_names')

is.null(scores(m2, display = 'species', tidy = T))



# PCA, RDA, CA, CCA
scores(rda(dune ~ 1), display = 'species') |> 
  as.data.frame() |> 
  summarise(sum1 = sum(PC1),
            sum2 = sum(PC2))

scores(m, tidy = T) |> 
  as_tibble() 

all.equal(
  rownames(vegan::scores(m, display = "sites")),
  rownames(vegan::scores(m, display = "species")))

m <- capscale(dune.dist ~ 1) 

all(is.na(vegan::scores(m, display = "species")))
