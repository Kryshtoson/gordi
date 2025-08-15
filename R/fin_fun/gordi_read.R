# gordi_read()

# arguments
gordi_read <- function(m,
                       env,
                       traits = NULL,
                       choices = 1:2,
                       scaling = 'symm',
                       correlation = F,
                       hill = F) {

  # type of ordination 
  type <- case_when(
    inherits(m, 'capscale') & !is.null(m$call$distance) & is.null(m$CCA) ~ 'PCoA (capscale)', #via capscale
    inherits(m, 'capscale') & !is.null(m$call$distance) ~ 'db-RDA (capscale)',
    inherits(m, 'rda') & !is.null(m$call$distance) & is.null(m$CCA) ~ 'PCoA (rda)', #if rda and distance but no constrainned variable PCoA
    inherits(m, 'rda') & !is.null(m$call$distance) ~ 'db-RDA (by RDA argument)',
    inherits(m, 'rda') & is.null(m$call$distance) & is.null(m$CCA) ~ 'PCA', inherits(m, 'rda') ~ 'RDA constrained',
    inherits(m, 'cca') & is.null(m$call$distance) & is.null(m$CCA) ~ 'CA', inherits(m, 'cca') ~ 'CCA',
    TRUE ~ paste(class(m), collapse = '/') # writes just one output
  )

  # create pass object 
  pass <- list(
    m = m,
    explained_variation = m$CA$eig/m$tot.chi,
    site_scores = as_tibble(as.data.frame(scores(m, scaling = scaling, choices = choices, correlation = correlation, hill = hill)$sites)),
    species_scores = as_tibble(as.data.frame(scores(m, scaling = scaling, choices = choices, correlation = correlation, hill = hill)$species)),
    env = env,
    traits = traits,
    choices = choices,
    type = type,
    species_names = as_tibble(as.data.frame(scores(m, scaling = scaling, choices = choices, correlation = correlation, hill = hill)$species), rownames = 'species_names')[1]
  )
  
  # Return pass object
  return(pass)
  
}