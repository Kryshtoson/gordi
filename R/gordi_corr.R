#' Correlate Environmental Variables with Ordination Axes
#'
#' @description
#' Passively fits environmental variables (vectors for numeric, centroids for factors) onto
#' the current ordination plot stored in \code{pass$plot} using \code{\link[vegan]{envfit}}. 
#' If the number of permutations is provided, it can also calculate R2, p-values (or adjusted p-values)
#' for the correlation. These values are then stored in \code{pass$corr_stats}.
#'
#' @details
#' The function calculates the correlation scores (vectors and factor centroids)
#' and optionally adjusts the p-values for multiple testing.
#'
#' \strong{Plotting:}
#' \itemize{
#'   \item \strong{Vectors} (for numeric variables) are plotted as arrows using \code{\link[ggplot2]{geom_segment}}.
#'   \item \strong{Centroids} (for factor levels) are plotted as points using \code{\link[ggplot2]{geom_point}}.
#' }
#'
#' Aesthetic arguments (\code{colour}, \code{fill}, \code{alpha}, etc.) can be set in two ways:
#' \enumerate{
#'   \item \strong{Mapping:} If the argument matches a column name in the internal \code{pass$corr_coords} data frame (e.g., 'covariate', 'variable_level', 'variable'), the aesthetic is mapped to that variable.
#'   \item \strong{Constant:} Otherwise, the value is used as a constant aesthetic (e.g., \code{colour = 'red'}).
#' }
#'
#' Use \code{\link{gordi_colour}()} or \code{\link{gordi_shape}()} or other similar functions immediately after \code{gordi_corr()}
#' to define custom legends and scales for mapped aesthetics. If you want to change `fill`, use \code{gordi_colour(fill = T)}.
#'
#' @param pass A list object produced by \code{\link{gordi_read}()} and containing the ordination results (\code{pass$m}) and environmental data (\code{pass$env}).
#' @param variables Character vector; names of the environmental variables (columns in \code{pass$env}) to fit onto the ordination. Required.
#' @param permutations Numeric; The number of permutations to use for significance testing with \code{\link[vegan]{envfit}}. Set to 0 to skip testing. Default is 0 (no test).
#' @param p_val_adjust Logical; Should p-values be adjusted for multiple comparisons? Uses an internal wrapper by David Zelený around \code{\link[stats]{p.adjust}}, called \code{p.adjust.envfit}. Default is \code{TRUE}.
#' @param p_val_adjust_method Character; The method for p-value adjustment. See \code{\link[stats]{p.adjust}}. Default is 'bonferroni'.
#' @param strata Grouping variable for permutations. See details in \code{\link[vegan]{envfit}}.
#' @param label Character; The column from \code{pass$corr_coords} to use for the text labels. Must be one of 'covariate', 'variable_level', or 'variable'. Default is 'covariate'.
#' @param show_label Logical; If \code{TRUE}, plot text labels for the vectors and factor centroids. Default is \code{FALSE}.
#' @param repel_label Logical; If \code{TRUE} (and \code{show_label} is \code{TRUE}), use \code{\link[ggrepel]{geom_text_repel}} to prevent label overlap. Default is \code{FALSE}.
#' @param colour Aesthetic; A column name to map to the colour of vectors (line) and factor centroids (point or point outline for shapes 21-25), or a constant colour.
#' @param fill Aesthetic; A column name to map to the fill colour of factor centroids (only works with shapes 21-25), or a constant colour.
#' @param alpha Aesthetic; A column name to map to transparency, or a constant numeric value (0 to 1).
#' @param linetype Aesthetic; A column name to map to linetype (vectors only), or a constant linetype (e.g., 1, 2, 'dotted', 'dashed').
#' @param linewidth Aesthetic; A column name to map to line width (vectors only), or a constant numeric value.
#' @param shape Aesthetic; A column name to map to point shape (factors only), or a constant numeric value (point shape).
#' @param size Aesthetic; A column name to map to point size (factors only), or a constant numeric value.
#' @param arrow_size Numeric; Size of the arrow heads for vectors, in 'cm'. Default is 0.3.
#'
#' @return The updated \code{pass} object with the environmental correlation layers added to \code{pass$plot} and the correlation statistics added to \code{pass$corr_stats}.
#'
#' @examples
#' \dontrun{
#' data(dune)
#' data(dune.env)
#' 
#' m <- capscale(sqrt(dune) ~ 1)
#' 
#' # 1. Plot correlated variables
#' gordi_read(m, env = dune.env) |> 
#'   gordi_corr(variables = c('A1', 'Management', 'Use'))
#' 
#' # 2. Check stats
#' o <- gordi_read(m, env = dune.env) |> 
#'   gordi_corr(variables = c('A1', 'Management', 'Use'))    
#'   
#' o$corr_stats
#'
#' # 3. Map colour to 'covariate' (Discrete)
#' gordi_read(m, env = dune.env) |> 
#'   gordi_corr(variables = c('A1', 'Management', 'Use'), colour = 'covariate') |> 
#'   gordi_colour(scale = 'discrete', family = 'viridis')
#' }
#' @importFrom vegan envfit scores
#' @importFrom dplyr as_tibble bind_rows mutate select starts_with relocate
#' @importFrom ggplot2 ggplot theme_bw labs theme element_text element_blank arrow unit aes geom_segment geom_point geom_text
#' @importFrom ggrepel geom_text_repel
#' @importFrom tibble tibble
#' @importFrom rlang sym expr
#' @importFrom stringr str_replace
#' @export
gordi_corr <- function(
    pass,
    variables = '',
    permutations = 0,
    p_val_adjust = TRUE,
    p_val_adjust_method = 'bonferroni',
    strata = NULL,
    scaling_coefficient = 0.9,
    label = c('covariate', 'variable_level', 'variable'),
    show_label = F,
    repel_label = F,
    colour = '',
    fill = '',
    alpha = '',
    linetype = '',
    linewidth = '',
    shape = '',
    size = '',
    arrow_size = 0.3) {
  
  if(is.null(variables) || length(variables) == 0 || (length(variables) == 1 && identical(variables, ''))) {
    stop('You must provide at least one variable you want to fit to the ordination.')}
  
  if(permutations == 0) {warning('You did not specify the number of permutations. No p-values were computed.')}
  
  label <- match.arg(label)
  
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
               env = pass$env[, variables, drop = FALSE],
               permutations = permutations,
               strata = strata,
               choices = pass$choices)
  if (p_val_adjust) {
    ef <- p.adjust.envfit(ef, method = p_val_adjust_method)
  } else {
    warning ('Significance values were not adjusted.')
  }
  
  
  # --- Collect plotting coordinates ---
  
  # prepare empty list to store variable coordinates
  coords <- list()
  
  # extract vector scores
  if (!is.null(ef$vectors)) {
    vec_tbl <- dplyr::as_tibble(vegan::scores(ef, display = 'bp'), rownames = 'covariate') |> 
      dplyr::mutate(score = 'vector',
                    variable = names(ef$vectors$r),
                    variable_level = names(ef$vectors$r)) |> 
      relocate(starts_with('Axis'), .before = 1) 
      
    coords <- append(coords, list(vec_tbl))
  } else {vec_tbl <- NULL}
  
  # extract factor scores
  if (!is.null(ef$factors)) {
    fct_tbl <- dplyr::as_tibble(vegan::scores(ef, display = 'cn'), rownames = 'covariate') |> 
      dplyr::mutate(score = 'factor',
                    variable = if(length(covariate) > 0 && length(names(ef$factors$r)) > 0) {
                      factor_names <- names(ef$factors$r)
                      # matrix of matches: rows = covariates, columns = factor names
                      matches <- outer(covariate, factor_names, Vectorize(function(cv, f) grepl(f, cv, fixed = TRUE)))
                      # pick first matching factor name per covariate
                      factor_names[apply(matches, 1, function(x) which(x)[1])]
                    } else {NA_character_},
                    variable_level = str_replace(covariate, variable, "")) |>  
      relocate(starts_with('Axis'), .before = 1)
    
    coords <- append(coords, list(fct_tbl))
  } else {fct_tbl <- NULL}
  
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
    if (length(ef$vectors$pvals) > 0) {pvals_vct <- ef$vectors$pvals} 
    
    vct_stats <- tibble::tibble(
      variable = names(ef$vectors$r),
      type = "vector",
      r2 = ef$vectors$r,
      p_value_adj = pvals_vct)
    
    if (p_val_adjust == F) {
      vct_stats <- vct_stats |> 
        rename('p_value' = 'p_value_adj')
    }
    
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
    
    if (p_val_adjust == F) {
      fct_stats <- fct_stats |> 
        rename('p_value' = 'p_value_adj')
    }
    
    stats <- append(stats, list(fct_stats))
  }
  
  pass$corr_stats <- dplyr::bind_rows(stats)
  
  
  
  # --- Ordination axis labels ---
  
  if(pass$type %in% c('DCA', 'NMDS')) {actual_labs <- paste0(pass$axis_names)} else 
  {actual_labs <- paste0(pass$axis_names, " (", round(pass$explained_variation[pass$choices]*100, 1), "%)")}
  
  
  
  # --- Add correlated variables to plot ---
  
  ### axis names used in spe_df
  names(pass$corr_coords)[2:3] <- paste0("Axis_corr", 1:2)
  
  ### create corr_df which is then called in the ggplot
  corr_df <- pass$corr_coords |> 
    select(Axis_corr1, Axis_corr2, score, covariate, variable_level, variable)
  
  pass$corr_coords <- pass$corr_coords |> 
    select(Axis_corr1, Axis_corr2, score, covariate, variable_level, variable)
  
  
  ### plot
  # Creates blank plot if this function is used as the first one after gordi_read()
  # or passes already existing plot
  
  if (is.null(pass$plot)) { # checks whether p exists in pass, if not it draws plot
    p <- ggplot() +
      theme_bw() +
      geom_vline(aes(xintercept = 0), linetype = 3, linewidth = 0.2, colour = 'gray15', alpha = 0.6) +
      geom_hline(aes(yintercept = 0), linetype = 3, linewidth = 0.2, colour = 'gray15', alpha = 0.6) +
      labs(x = actual_labs[1], y = actual_labs[2]) +
      theme(
        text = element_text(size = 15),
        panel.grid = element_blank(),
        legend.justification = c(1, 1))
  } else {p <- pass$plot}
  
  
  ### Detect mapped vs constant aesthetics
  # colour
  # if colour != "" AND ALSO colour represents a colname present in corr_df, then map_colour is TRUE, otherwise is FALSE
  map_colour <- !identical(colour, '') && (colour %in% names(corr_df))
  # if map_colour is FALSE AND ALSO the thing inputed in arguments is a HEX code or is included in colours() or in palette() (word or number), then use it as const_colour
  const_colour <- !map_colour && (grepl("^#(?:[A-Fa-f0-9]{6}[A-Fa-f0-9]{3})$", colour) || colour %in% grDevices::colours()) || (is.character(colour) && colour %in% palette()) || (is.numeric(colour) && colour %in% seq_along(palette()))
  # fill
  map_fill <- !identical(fill, '') && (fill %in% names(corr_df))
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
  
  ### Set scaling coefficient
  # extract plot frame size (x and y axis lengths)
  if (!is.null(vec_tbl)) {
    p_build <- ggplot_build(p)
    
    plot_range <- c(xmin_plot = p_build$layout$panel_params[[1]]$x.range[1],
                    xmax_plot = p_build$layout$panel_params[[1]]$x.range[2],
                    ymin_plot = p_build$layout$panel_params[[1]]$y.range[1],
                    ymax_plot = p_build$layout$panel_params[[1]]$y.range[2])
    
    corr_range <- c(xmin_pred = min(corr_df[corr_df$score == 'vector',1]),
                    xmax_pred = max(corr_df[corr_df$score == 'vector',1]),
                    ymin_pred = min(corr_df[corr_df$score == 'vector',2]),
                    ymax_pred = max(corr_df[corr_df$score == 'vector',2]))
    
    coef <- (max(abs(plot_range)) / max(abs(corr_range))) * scaling_coefficient
  } else {coef <- 1}
  
  
  ### Prepare aes arguments for geom_segment()
  # Start with fixed x/y for the base (0,0) and end at the species scores
  aes_args_segment <- list(
    x = 0, y = 0,
    xend = expr(Axis_corr1 * !!coef),
    yend = expr(Axis_corr2 * !!coef)
  )
  
  if(map_colour) aes_args_segment$colour <- sym(colour)
  if(map_alpha) aes_args_segment$alpha <- sym(alpha)
  if(map_linewidth) aes_args_segment$linewidth <- sym(linewidth)
  if(map_linetype) aes_args_segment$linetype <- sym(linetype)
  
  
  ### Prepare constant arguments for geom_point() (mapped first, then defaults if nothing)
  const_args_segment <- list()
  
  # Add constant arguments for geom_segment() if not mapped
  # colour
  if(!map_colour){
    if(!identical(colour, '')) {const_args_segment$colour <- colour} else {const_args_segment$colour <- 'gray20'}}
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
  
  if (!is.null(corr_df |> filter(score == 'vector'))) {
    p <- p + do.call(geom_segment,
                     c(list(data = corr_df |> filter(score == 'vector'),
                            mapping = do.call(aes, aes_args_segment)),
                       const_args_segment)) 
    
    if (isTRUE(show_label)){
      
      aes_args_text <- list(
        x = expr(Axis_corr1 * !!coef),
        y = expr(Axis_corr2 * !!coef),
        label = sym(label)
      )
      
      if(map_colour) aes_args_text$colour <- sym(colour)
      if(map_alpha) aes_args_text$alpha <- sym(alpha)
      
      ### Prepare constant arguments for geom_text() (mapped first, then defaults if nothing)
      const_args_text <- list(
        show.legend = F
      )

      if(!map_colour){ if(!identical(colour, '')) {const_args_text$colour <- colour} else {const_args_text$colour <- 'gray20'}}
      if(!map_alpha){ if(!identical(alpha, '')) {const_args_text$alpha <- alpha} else {const_args_text$alpha <- 1}}
      
      
        if (isTRUE(repel_label)){
        p <- p + do.call(ggrepel::geom_text_repel,
                           c(list(data = corr_df |> filter(score == 'vector'),
                                  mapping = do.call(aes, aes_args_text)),
                             const_args_text))
          
        } else {
        p <- p + do.call(geom_text,
                         c(list(data = corr_df |> filter(score == 'vector'),
                                mapping = do.call(aes, aes_args_text)),
                           const_args_text))
      }
    }
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
    if(!identical(colour, '')) {const_args_point$colour <- colour} else {const_args_point$colour <- 'gray62'}}
  if(!map_fill){
    if(!identical(fill, '')) {const_args_point$fill <- fill} else {const_args_point$fill <- 'gray62'}}
  if(!map_alpha){
    if(!identical(alpha, '')) {const_args_point$alpha <- alpha} else {const_args_point$alpha <- 1}}
  if(!map_shape){
    if(!identical(shape, '')) {const_args_point$shape <- shape} else {const_args_point$shape <- 16}}
  if(!map_size){
    if(!identical(size, '')) {const_args_point$size <- size} else {const_args_point$size <- 2}}
  
  # add to plot
  if (!is.null(corr_df |> filter(score == 'factor'))) {
    
    p <- p + do.call(geom_point,
                     c(list(data = corr_df |> filter(score == 'factor'),
                            mapping = do.call(aes, aes_args_point)),
                       const_args_point)) 
    
    
    if (isTRUE(show_label)){
      
      aes_args_text <- list(
        x = sym("Axis_corr1"),
        y = sym("Axis_corr2"),
        label = sym(label)
      )
      
      if(map_colour) aes_args_text$colour <- sym(colour)
      if(map_alpha) aes_args_text$alpha <- sym(alpha)
      
      ### Prepare constant arguments for geom_text() (mapped first, then defaults if nothing)
      const_args_text <- list(
        show.legend = F
      )
      
      if(!map_colour){ if(!identical(colour, '')) {const_args_text$colour <- colour} else {const_args_text$colour <- 'gray20'}}
      if(!map_alpha){ if(!identical(alpha, '')) {const_args_text$alpha <- alpha} else {const_args_text$alpha <- 1}}
      
      
      if (isTRUE(repel_label)){
        p <- p + do.call(ggrepel::geom_text_repel,
                         c(list(data = corr_df |> filter(score == 'factor'),
                                mapping = do.call(aes, aes_args_text)),
                           const_args_text))
      } else {
        p <- p + do.call(ggrepel::geom_text_repel,
                         c(list(data = corr_df |> filter(score == 'factor'),
                                mapping = do.call(aes, aes_args_text)),
                           const_args_text))
      }
    }
  }
  
  
  pass$plot <- p
  
  
  # --- Return object ---
  
  return(pass)
}