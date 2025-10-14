library(tidyverse)
library(vegan)

devtools::document()

data(dune)
data("dune.env")

dune.env |> as_tibble()

m <- capscale(sqrt(dune) ~ 1, env = dune.env, sqrt.dist = T)

envfit(m ~ Use + Management + A1, data = dune.env, permutations = 0, choices = 1:2) |> 
  p.adjust.envfit(method = 'holm') |> 
  scores(display = 'cn')

o <- envfit(m, env = dune.env[, c('A1', 'Management', 'Use', 'Manure')], data = dune.env, perm = 0, choices = 1:2) 

o

scores(o, display = 'bp') |> as_tibble(rownames = 'variable')
scores(o, display = 'cn') |> as_tibble(rownames = 'variable')

o$vectors$r
o$vectors$pvals

o$factors$r
o$factors$pvals

oo <- p.adjust.envfit(o)

scores(oo, display = 'bp') |>
  as_tibble(rownames = 'covariate') |>
  mutate(type = 'vector',
         name = names(oo$vectors$r))

scores(oo, display = 'cn') |>
  as_tibble(rownames = 'covariate') |>
  mutate(
    type = 'factor',
    name = {
      factor_names <- names(oo$factors$r)
      # matrix of matches: rows = covariates, columns = factor names
      matches <- outer(covariate, factor_names, Vectorize(function(cv, f) grepl(f, cv, fixed = TRUE)))
      # pick first matching factor name per covariate
      factor_names[apply(matches, 1, function(x) which(x)[1])]
    }
  )

         
oo$vectors$pvals

oo$vectors$r
oo$vectors$pvals

oo$factors$r
oo$factors$pvals

??base::logic



o$vectors |> is.null()
oo$vectors$pvals
oo$factors$pvals

o$vectors$arrows |> as_tibble()

o$vectors$r
o$factors$r

scores(o, display = 'cn')

m |> 
  gordi_read()


# gordi_corr --------------------------------------------------------------

gordi_corr <- function(
    pass,
    variables = '',
    permutations = 0,
    p_val_adjust = TRUE,
    p_val_adjust_method = 'bonferroni',
    strata = NULL,
    colour = '', 
    fill = '', 
    alpha = '',
    linetype = '',
    linewidth = '',
    shape = '',
    size = '',
    arrow_size = 0.3) {
  
  if(is.null(variables) || length(variables) == 0 || (length(variables) == 1 && identical(variables, ''))) {
    stop('You must provide at least one variable you want correlate with the ordination axes.')}
  
  if(permutations == 0) {warning('You did not specify the number of permutations. No p-values were computed.')}
  
  # ---- Helper function by D. Zeleny: p.adjust.envfit ----
  # Function p.adjust.envfit
  # Calculates adjusted P values for results stored in envfit object,
  # created using envfit function from vegan, which fits supplementary variables
  # onto axes of unconstrained ordination. 
  # Arguments: 
  # x - envfit object
  # method - method for correction of multiple testing issue (default = 'bonferroni',
  #          see ''?p.adjust' for more options)
  # n - optional, number of tests for which to correct; if not given, the number is
  #          taken as the number of tests conducted by envfit function (for both vectors and factors).
  # Author: David Zeleny
  p.adjust.envfit <- function (x, method = 'bonferroni', n)
  {
    x.new <- x
    if (!is.null (x$vectors)) pval.vectors <- x$vectors$pvals else pval.vectors <- NULL
    if (!is.null (x$factors)) pval.factors <- x$factors$pvals else pval.factors <- NULL
    if (missing (n)) n <- length (pval.vectors) + length (pval.factors)
    if (!is.null (x$vectors)) x.new$vectors$pvals <- p.adjust (x$vectors$pvals, method = method, n = n)
    if (!is.null (x$factors)) x.new$factors$pvals <- p.adjust (x$factors$pvals, method = method, n = n)
    message('Adjustment of significance by ', method, ' method')
    return (x.new)
  }
  
  # --- Continuation ---
  
  # fit environmental variables onto ordination 
  ef <- envfit(ord = pass$m,
               env = pass$env[ , variables, drop = FALSE],
               permutations = permutations,
               strata = strata,
               choices = pass$choices)
  if (p_val_adjust) {
    ef <- p.adjust.envfit(ef, method = p_val_adjust_method)
  } else {
    warning ('Significance values were not adjusted')
  }
  
  
  # --- Collect plotting coordinates ---
  
  # prepare empty list to store variable coordinates
  coords <- list()
  
  # extract vector scores
  if (!is.null(ef$vectors)) {
    vec_tbl <- dplyr::as_tibble(vegan::scores(ef, display = 'bp'), rownames = 'covariate') |> 
      dplyr::mutate(type = 'vector',
                    name = names(ef$vectors$r))
    coords <- append(coords, list(vec_tbl))
  }
  
  # extract factor scores
  if (!is.null(ef$factors)) {
    fct_tbl <- dplyr::as_tibble(vegan::scores(ef, display = 'cn'), rownames = 'covariate') |> 
      dplyr::mutate(type = 'factor',
                    name = if(length(covariate) > 0 && length(names(ef$factors$r)) > 0) {
                      factor_names <- names(ef$factors$r)
                      # matrix of matches: rows = covariates, columns = factor names
                      matches <- outer(covariate, factor_names, Vectorize(function(cv, f) grepl(f, cv, fixed = TRUE)))
                      # pick first matching factor name per covariate
                      factor_names[apply(matches, 1, function(x) which(x)[1])]
                    } else {NA_character_})
    coords <- append(coords, list(fct_tbl))
  }
  
  pass$corr_coords <- dplyr::bind_rows(coords)
  
  
  
  # --- Collect stats ---
  
  # create empty list for stats
  stats <- list()
  
  if (!is.null(ef$vectors)) {
    # 1. Determine the length of the vector results
    N_vct <- length(ef$vectors$r)
    
    # 2. Define the NA vector as the default
    pvals_vct <- rep(NA_real_, N_vct)
    
    # 3. OVERWRITE only if p-values exist
    if(length(ef$vectors$pvals) > 0) {pvals_vct <- ef$vectors$pvals} 
    
    vct_stats <- tibble::tibble(
        variable = names(ef$vectors$r),
        type = "vector",
        r2 = ef$vectors$r,
        p_value_adj = pvals_vct)

    stats <- append(stats, list(vct_stats))
  }
  
  
  
  if (!is.null(ef$factors)) {
    # 1. Determine the length of the factor results
    N_fct <- length(ef$factors$r)
    
    # 2. Define the NA vector as the default
    pvals_fct <- rep(NA_real_, N_fct)
    
    # 3. OVERWRITE only if p-values exist
    if(length(ef$factors$pvals) > 0) {pvals_fct <- ef$factors$pvals} 
    
    fct_stats <- tibble::tibble(
        variable = names(ef$factors$r),
        type = "factor",
        r2 = ef$factors$r,
        p_value_adj = pvals_fct)
    
    stats <- append(stats, list(fct_stats))
  }
  
  pass$corr_stats <- dplyr::bind_rows(stats)
  
  
  
  # --- Ordination axis labels ---
  
  if(pass$type %in% c('DCA', 'NMDS')) {actual_labs <- paste0(pass$axis_names)} else 
  {actual_labs <- paste0(pass$axis_names, " (", round(pass$explained_variation[pass$choices]*100, 2), "%)")}
  
  
  
  # --- Add correlated variables to plot ---
  
  ### axis names used in spe_df
  names(pass$corr_coords)[2:3] <- paste0("Axis_corr", 1:2)
  
  ### create corr_df which is then called in the ggplot
  corr_df <- pass$corr_coords
  
  ### plot
  # Creates blank plot if this function is used as the first one after gordi_read()
  # or passes already existing plot
  
  if (is.null(pass$plot)) { # checks whether p exists in pass, if not it draws plot
    p <- ggplot() +
      theme_bw() +
      labs(x = actual_labs[1], y = actual_labs[2]) +
      theme(
        text = element_text(size = 15),
        panel.grid = element_blank(),
        legend.justification = c(1, 1))
  } else {p <- pass$plot}
  
  
  ### Detect mapped vs constant aesthetics
  # colour
  # if colour != "" AND ALSO colour represents a colname present in spe_df, then map_colour is TRUE, otherwise is FALSE
  map_colour <- !identical(colour, '') && has_name(corr_df, colour) 
  # if map_colour is FALSE AND ALSO the thing inputed in arguments is a HEX code or is included in colours() or in palette() (word or number), then use it as const_colour
  const_colour <- !map_colour && (grepl("^#(?:[A-Fa-f0-9]{6}[A-Fa-f0-9]{3})$", colour) || colour %in% grDevices::colours()) || (is.character(colour) && colour %in% palette()) || (is.numeric(colour) && colour %in% seq_along(palette()))
  # fill
  map_fill <- !identical(fill, '') && has_name(corr_df, fill)
  const_fill <- !map_fill && (grepl("^#(?:[A-Fa-f0-9]{6}[A-Fa-f0-9]{3})$", fill) || fill %in% grDevices::colours()) || (is.character(fill) && fill %in% palette()) || (is.numeric(colour) && colour %in% seq_along(palette()))
  # alpha
  map_alpha <- !identical(alpha, '') && has_name(corr_df, alpha)
  const_alpha <- !map_alpha && is.numeric(alpha)
  # linetype
  map_linetype <- !identical(linetype, '') && has_name(corr_df, linetype)
  const_linetype <- !map_linetype && (is.numeric(linetype) || !identical(linetype,''))
  # linewidth
  map_linewidth <- !identical(linewidth, '') && has_name(corr_df, linewidth)
  const_linewidth <- !map_linewidth && is.numeric(linewidth)
  # shape
  map_shape <- !identical(shape, '') && has_name(corr_df, shape)
  const_shape <- !map_shape && is.numeric(shape)
  # size
  map_size <- !identical(size, '') && has_name(corr_df, size)
  const_size <- !map_size && is.numeric(size)
  


  ### Prepare aes arguments for geom_segment()
  # Start with fixed x/y for the base (0,0) and end at the species scores
  aes_args_segment <- list(
    x = 0, y = 0,
    xend = sym("Axis_corr1"),
    yend = sym("Axis_corr2"))
  
  if(map_colour) aes_args_segment$colour <- sym(colour)
  if(map_alpha) aes_args_segment$alpha <- sym(alpha)
  if(map_linewidth) aes_args_segment$linewidth <- sym(linewidth)
  if(map_linetype) aes_args_segment$linetype <- sym(linetype)


  ### Prepare constant arguments for geom_point() (mapped first, then defaults if nothing)
  const_args_segment <- list()
  
  # Add constant arguments for geom_segment() if not mapped
  # colour
  if(!map_colour){
    if(!identical(colour, '')) {const_args_segment$colour <- colour} else {const_args_segment$colour <- 'gray10'}}
  # alpha
  if(!map_alpha){
    if(!identical(alpha, '')) {const_args_segment$alpha <- alpha} else {const_args_segment$alpha <- 1}}
  # linetype
  if(!map_linetype){
    if(!identical(linetype, '')) {const_args_segment$linetype <- linetype} else {const_args_segment$linetype <- 1}}
  # linewidth
  if(!map_linewidth){
    if(!identical(linewidth, '')) {const_args_segment$linewidth <- linewidth} else {const_args_segment$linewidth <- 1}}
  # arrow_size (does not make sense to have map_arrow)
  const_args_segment$arrow <- arrow(
    length = unit(
      if (!identical(arrow_size, '')) as.numeric(arrow_size) else 0.3, "cm"))
  
 
  
  
  ### add to plot  
  
  if (!is.null(corr_df |> filter(type == 'vector'))) {
    p <- p + do.call(geom_segment,
                     c(list(data = corr_df |> filter(type == 'vector'),
                            mapping = do.call(aes, aes_args_segment)),
                       const_args_segment)) 
    
    p <- p + ggplot2::labs(colour = 'Covariate',
                           fill = 'Covariate',
                           alpha = 'Covariate',
                           linewidth = 'Covariate',
                           linetype = 'Covariate')
  }
  
  ### Categorical variables
  
  ### Prepare aes arguments for geom_point()
  aes_args_point <- list(
    x = sym("Axis_corr1"),
    y = sym("Axis_corr2")
  )
  
  if(map_colour) aes_args_point$colour <- sym(colour)
  if(map_fill) aes_args_point$fill <- sym(fill)
  if(map_alpha) aes_args_point$alpha <- sym(alpha)
  if(map_shape) aes_args_point$shape <- sym(shape)
  if(map_size) aes_args_point$size <- sym(size)
  
  
  
  ### Prepare constant arguments for geom_point()
  const_args_point <- list()
  
  if(!map_colour){
    if(!identical(colour, '')) {const_args_point$colour <- colour} else {const_args_point$colour <- 'gray30'}}
  if(!map_fill){
    if(!identical(fill, '')) {const_args_point$fill <- fill} else {const_args_point$fill <- 2}}
  if(!map_alpha){
    if(!identical(alpha, '')) {const_args_point$alpha <- alpha} else {const_args_point$alpha <- 1}}
  if(!map_shape){
    if(!identical(shape, '')) {const_args_point$shape <- shape} else {const_args_point$shape <- 16}}
  if(!map_size){
    if(!identical(size, '')) {const_args_point$size <- size} else {const_args_point$size <- 2}}
  
  # add to plot
  if (!is.null(corr_df |> filter(type == 'factor'))) {
    p <- p + do.call(geom_point,
                     c(list(data = corr_df |> filter(type == 'factor'),
                            mapping = do.call(aes, aes_args_point)),
                       const_args_point)) 
    
    p <- p + ggplot2::labs(
      colour = 'Covariate',
      fill = 'Covariate',
      alpha = 'Covariate',
      shape = 'Covariate',
      size = 'Covariate')
  }

  pass$plot <- p
  
  
  # --- Return object ---
  
  return(pass)
}

m |> 
  gordi_read(env = dune.env, scaling = 'species', correlation = T, hill = T, choices = 1:2) |> 
  gordi_species() |> 
  gordi_corr(variables = c('Manure', 'Use', 'A1'), p_val_adjust = T,
             p_val_adjust_method = 'holm', perm = 999,
             shape = 'covariate', size = 5) |> 
  gordi_shape(scale = 'discrete', values = c(0,1,2,3,4,15,16,17))


vars <- c('Management', 'Use')

is.null(vars) || length(vars) == 0 || (length(vars) == 1 && identical(vars, ''))

obj <- m |> 
  gordi_read(env = dune.env, scaling = 'species', correlation = T, hill = T, choices = 1:2) |> 
  gordi_species() |> 
  gordi_colour(scale = 'discrete', family = 'manual',
               values = c('red', 'blue', 'blue', 'blue')) |> 
  gordi_plot() +
  labs(colour = 'Covariates')

obj



# a taky labels

# pridat correlations (covariates) to gordi_label()