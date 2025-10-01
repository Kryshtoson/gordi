#goodness of fit

gordi_fit <- function(pass,
                      spe = NULL,
                      what = c('species', 'sites'),
                      display = c('sites', 'lc'),
                      choices = 1:2,
                      summarize = T,
                      slice_max = NULL,
                      permutations = 0
){
  
  m <- pass$m
  pass$spe <- spe #add spe to pass list
  #env is already in pass list
  
  what <- match.arg(what)
  display <- match.arg(display)
  
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
    type %in% c('CCA') ~ 'CCA',
    type %in% c('CA', 'DCA', 'PCA') ~ 'CA',
    TRUE ~ NA
  )
  
  spe_fitted <- NULL
  sites_fitted <- NULL
  
  if (!is.na(model)){
    if (what == 'species'){
    goodness_fit <- goodness(m, display = 'species', model = model, choices = choices, summarize = summarize)|>
      as_tibble()
    spe_fitted <- bind_cols(pass$species_names, goodness_fit, pass$species_scores)
    } else if (what == 'sites'){
      goodness_fit <- goodness(m, display = 'sites', model = model, choices = choices, summarize = summarize)|>
        as_tibble()
      sites_fitted <- bind_cols(pass$env[1], goodness_fit, pass$site_scores)
    
    }
  }
  
  if (inherits(m, 'capscale')){
    if (what == 'species'){
      goodness_fit <- envfit(m, env = spe, display = display, choices = choices, permutations = permutations)
      r <- goodness_fit$vectors$r
      spe_fitted <- bind_cols(pass$species_names, r, pass$species_scores)
    } else if (what == 'sites'){
      goodness_fit <- envfit(m, env = pass$env, display = display, choices = choices, permutations = permutations)
      r <- goodness_fit$vectors$r
      sites_fitted <- bind_cols(pass$env[1], r, pass$site_scores)
    }
  } 
  
  if (!is.null(slice_max)){
    if (!is.null(spe_fitted)){
    spe_fitted <- spe_fitted|>
      slice_max(order_by = spe_fitted[2], n = slice_max)
    } 
    if (!is.null(sites_fitted)){
      sites_fitted <- sites_fitted|>
        slice_max(order_by = sites_fitted[2], n = slice_max)
    }
  }
  
  if (!is.null(spe_fitted)){
    spe_fitted|>
      select(-(1:2)) -> pass$species_scores
    spe_fitted[1] -> pass$species_names
  }
  if (!is.null(sites_fitted)){
    pass$site_scores <- sites_fitted
  }
  
  return(pass)
}

gordi_read(m, env) -> o
envfit(m, env = spe, display = 'si', choices = 1:2, permutations = 0) -> ef
bind_cols(o$species_names, ef$vectors$r, o$species_scores,) -> result
result|>
  slice_max(order_by = result[[2]], n = 30)|>
  select(-(1:2))

gordi_read(m, env)|>
  gordi_fit(spe = spe, what = 'species', slice_max = 30, display = 'sites')

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

model <- case_when(
  type %in% c('CCA', 'RDA') ~ 'CCA',
  type %in% c('CA', 'DCA', 'PCA') ~ 'CA',
  TRUE ~ NA
)

envfit(m, env = env, display = 'sites', choice = 1:2)

goodness_fit <- goodness(m, display = 'species', model = model, choices = 1:2, summarize = T)|>
  as_tibble()
