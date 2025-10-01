library(tidyverse)
library(vegan)

devtools::document()

data(dune)
data("dune.env")

dune.env |> as_tibble()

m <- capscale(sqrt(dune) ~ 1, env = dune.env, sqrt.dist = T)

envfit(m ~ Use + Management + A1, data = dune.env, perm = 999, choices = 2:3) |> 
  p.adjust.envfit(method = 'holm') |> 
  scores(display = 'cn')

o <- envfit(m ~ Use + Management + A1, data = dune.env, perm = 999, choices = 2:3) 

scores(oo, display = 'bp') |> as_tibble(rownames = 'variable')
scores(oo, display = 'cn') |> as_tibble(rownames = 'variable')

table(dune.env$Use)
table(dune.env$Management)

o$factors$r
o$factors$pvals

oo <- p.adjust.envfit(o)

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
    p_val_adjust_method = 'bonferroni',
    strata = NULL,
    colour = '', 
    fill = '', # jeste neni
    alpha = '',
    linetype = '',
    linewidth = '',
    arrow_size = 0.3) {
  
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
    cat ('Adjustment of significance by', method, 'method')
    return (x.new)
  }
  
  # ---- Continuation ----
  
  # fit environmental variables onto ordination 
  ef <- envfit(ord = pass$m,
               env = pass$env[ , variables, drop = FALSE],
               permutations = permutations,
               strata = strata,
               choices = pass$choices) |> 
    p.adjust.envfit(method = p_val_adjust_method)
  
  # --- Collect plotting coordinates ---
  
  # prepare empty list to store variable coordinates
  coords <- list()
  
  # extract vector scores
  if (!is.null(ef$vectors)) {
    vec_tbl <- as_tibble(scores(ef, display = 'bp'), rownames = 'variable') |> 
      dplyr::mutate(type = 'vector')
    coords <- append(coords, list(vec_tbl))
    }
  
  # extract factor scores
  if (!is.null(ef$factors)) {
    fct_tbl <- as_tibble(scores(ef, display = 'cn'), rownames = 'variable') |> 
      dplyr::mutate(type = 'factor')
    coords <- append(coords, list(fct_tbl))
  }
  
  pass$corr_coords <- dplyr::bind_rows(coords)
  
  
  
  # --- Collect stats ---
  
  # create empty list for stats
  stats <- list()
  
  if (!is.null(ef$vectors)) {
    stats <- append(stats, list(
      tibble::tibble(variable = names(ef$vectors$r),
                     type = "vector",
                     r2 = ef$vectors$r,
                     p_value_adj = if (!is.null(ef$vectors$pvals)) ef$vectors$pvals else NA_real_)
    ))
  }
  
  if (!is.null(ef$factors)) {
    stats <- append(stats, list(
      tibble::tibble(variable = names(ef$factors$r),
                     type = "factor",
                     r2 = ef$factors$r,
                     p_value_adj = if (!is.null(ef$factors$pvals)) ef$factors$pvals else NA_real_)
    ))
  }
  
  pass$corr_stats <- dplyr::bind_rows(stats)
  
  ### ordination axis labels
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
  # alpha
  map_alpha <- !identical(alpha, '') && has_name(corr_df, alpha)
  const_alpha <- !map_alpha && is.numeric(alpha)
  # linetype
  map_linetype <- !identical(linetype, '') && has_name(corr_df, linetype)
  const_linetype <- !map_linetype && (is.numeric(linetype) || !identical(linetype,''))
  # linewidth
  map_linewidth <- !identical(linewidth, '') && has_name(corr_df, linewidth)
  const_linewidth <- !map_linewidth && is.numeric(linewidth)
  
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
    if(!identical(colour, '')) {const_args_segment$colour <- colour} else {const_args_segment$colour <- 2}}
  # alpha
  if(!map_alpha){
    if(!identical(alpha, '')) {const_args_segment$alpha <- alpha} else {const_args_segment$alpha <- 1}}
  # linetype
  if(!map_linetype){
    if(!identical(linetype, '')) {const_args_segment$linetype <- linetype} else {const_args_segment$linetype <- 1}}
  # linewidth
  if(!map_linewidth){
    if(!identical(linewidth, '')) {const_args_segment$linewidth <- linewidth} else {const_args_segment$linewidth <- 0.7}}
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
  }
  
  ### Categorical variables
  
  ### Prepare aes arguments for geom_point()
  aes_args_point <- list(
    x = sym("Axis_corr1"),
    y = sym("Axis_corr2")
  )
  
  if(map_colour) aes_args_point$colour <- sym(colour)
  if(map_alpha) aes_args_point$alpha <- sym(alpha)
  
  ### Prepare constant arguments for geom_point()
  const_args_point <- list()
  
  if(!map_colour){
    if(!identical(colour, '')) {const_args_point$colour <- colour} else {const_args_point$colour <- 2}}
  if(!map_alpha){
    if(!identical(alpha, '')) {const_args_point$alpha <- alpha} else {const_args_point$alpha <- 1}}
 
  # add to plot
  if (!is.null(corr_df |> filter(type == 'factor'))) {
    p <- p + do.call(geom_point,
                     c(list(data = corr_df |> filter(type == 'factor'),
                            mapping = do.call(aes, aes_args_point)),
                       const_args_point)) 
  }

  pass$plot <- p
  
  
  # --- Return object ---
  
  return(pass)
  
}


obj <- m |> 
  gordi_read(env = dune.env, choices = 1:2) |> 
  gordi_species() |> 
  gordi_corr(variables = c('Use', 'Management', 'A1'), permutations = 999,
             colour = 'green', arrow_size = 0.5, linewidth = 2)


obj


# jeste map_shape, map_size and map_fill 
# aby fungovaly pro geom_point

# mozna by bylo dobre, aby colour a fill existovaly zvlast pro point a segment

# a taky labels