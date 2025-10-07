#goodness of fit

gordi_fit <- function(pass,
                      spe = NULL,
                      what = c('species', 'sites'),
                      display = c('sites', 'lc'),
                      choices = 1:2,
                      summarize = T,
                      slice_max = NULL,
                      abs_frequency = NULL,
                      relat_frequency = NULL,
                      permutations = 0
){
  
  m <- pass$m
  pass$spe <- spe #add spe to pass list
  #env is already in pass list
  
  what <- match.arg(what)
  display <- match.arg(display)
  
  if (what == 'sites' && !is.null(abs_frequency)||!is.null(relat_frequency)){
    stop("Cannot calculate frequency of occurence for sites. For sites, only filtering based on goodness of fit is available. Please specify `slice_max` if you would like to filter based on goodness of fit.")
  }
  
  if (!is.null(what)){
    warning("If you want to filter based on goodness of fit, this function chooses between ")
  }
  
  spe %>%
    mutate_all(~ if_else(. > 0, 1, 0))|> 
    summarize(across(everything(), sum))|>
    pivot_longer(everything(), names_to = 'species_names', values_to = 'abs_freq')|>
    mutate(relat_freq = abs_freq/nrow(spe))-> frequency
  
  bind_cols(pass$species_names, pass$species_scores)|>
    left_join(frequency, by = 'species_names') -> spe_freq
  
  if (!is.null(abs_frequency) && is.null(slice_max)){  
    #select species scores based on species frequency
    spe_freq|>
      filter(abs_freq > abs_frequency) -> species_filtered
  }
  
  if (!is.null(relat_frequency) && is.null(slice_max)){
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
    type %in% c('CA', 'DCA', 'PCA') ~ 'CA',
    TRUE ~ NA
  )
  
  spe_fitted <- NULL
  sites_fitted <- NULL
  
  if (!is.na(model)){
    warning("Using function `goodness()` to calculate goodness of fit.")
    if (what == 'species' && !is.null(slice_max)){
    goodness_fit <- goodness(m, display = 'species', model = model, choices = choices, summarize = summarize)|>
      as_tibble()
    spe_fitted <- bind_cols(pass$species_names, goodness_fit, pass$species_scores)
    } else if (what == 'sites' && !is.null(slice_max)){
      goodness_fit <- goodness(m, display = 'sites', model = model, choices = choices, summarize = summarize)|>
        as_tibble()
      sites_fitted <- bind_cols(pass$env[1], goodness_fit, pass$site_scores)
    
    }
  }
  
  if (inherits(m, 'capscale')){
    warning("Using function `envfit()` to calculate goodness of fit (R²).")
    if (what == 'species'){
      goodness_fit <- envfit(m, env = spe, display = display, choices = choices, permutations = permutations)
      r <- goodness_fit$vectors$r |> as_tibble()
      spe_fitted <- bind_cols(pass$species_names, r, pass$species_scores)
    } else if (what == 'sites'){
      goodness_fit <- envfit(m, env = pass$env, display = display, choices = choices, permutations = permutations)
      r <- goodness_fit$vectors$r|> as_tibble()
      sites_fitted <- bind_cols(pass$env[1], r, pass$site_scores)
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
