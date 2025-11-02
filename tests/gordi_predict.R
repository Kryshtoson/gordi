#' Extracts predictor scores and relevant stuff and creates plot - so far working for capscale
#' 
#' @description
#' [gordi_predict()] takes the result of [gordi_read()] and creates plot with
#' predictor arrows (for continuous variables) and points (for categorical
#' variables, i.e., centroids). So far, it can work with continuous explanatory 
#' variables. It also can't display explanatory variables passively (i.e., not 
#' involved in the ordination model).
#' 
#' There will be a special function for this once, maybe something as
#' `gordi_corr()`. Similarly to [gordi_species()] and [gordi_sites()],
#' you can also set a wide range of graphing parameters, such as colour, fill, size,
#' shape, alpha, stroke, and more traditional ggplot arguments, which can read
#' both, static and dynamic variable (e.g., 'red' or 'elevation').
#' 
#' @param pass Object from [gordi_read()] function.
#' @param label Logical; default = F, whether to display label by each point/arrow 
#'   or not. Labels use the `predictor` name for arrows and the combined 
#'   `predictor: level` for centroids. If you want to customize the labels, you can 
#'   with [gordi_label()] function, which overrides this setting.
#' @param colour Colour can be defined statically as word from the [colours()]
#'   list (e.g. 'red'), HEX code (e.g. #5d782e), or number from [palette()] (e.g. 3).
#'   It can also defined dynamically (according to some variable, e.g. elevation
#'   or vegetation type). This variable has to be present in the `env` table and its
#'   name has to be written in "quotation marks". Default colour of arrows is 2, 
#'   and 2 for points.
#' @param alpha Transparency of symbols/arrows. Numeric value. Default is 1 (for arrows)
#'   and 1 (for points). Can be set statically or dynamically (use "").       
#' @param arrow_size Numeric; length of arrow in cm. Default is 0.3 cm. 
#' @param linewidth Self-explanatory. Can be numerical value. Used for arrows. See
#'   [aes_linetype_size_shape()] for more details. Default is 0.7.
#' @param linetype Changes type of line used in arrows. Can be specified with either
#'   an integer (0-6), a name (0 = blank, 1 = solid, 2 = dashed, 3 = dotted,
#'   4 = dotdash, 5 = longdash, 6 = twodash). See [aes_linetype_size_shape()] for 
#'   more details. Default is 1.
#' @param fill Fill colour for centroid points. Can be defined statically (e.g., 'red', '#ff0000', or 3) 
#'   or dynamically (e.g., "elevation"). Default is 3.
#' @param shape Shape of centroid points. Can be specified with either 
#'   an integer (0-25) or a name. See [aes_linetype_size_shape()] for more details. 
#'   Default is 16.
#' @param size Size of centroid points. Numeric value. Can be set statically 
#'   or dynamically (use ""). Default is 4.
#' @param scaling_coefficient Numeric; edits the lengts of predictor arrow. Default value
#'   is `0.9` which means that the predictor arrow will be as long as 0.9 of
#'   the longest axis displayed in the plotframe.
#' @param repel_label Logical; repels labels of predictors for better readability.
#'   Default is F. If you want to customize the labels, you can do it with
#'   [gordi_label()] function, which overrides this setting.    
#' 
#' 
#' 
#' @return The input object `pass` with an updated `plot` element that includes
#'   predictor arrows and labels.
#' 
#' 
#' @importFrom rlang is_empty sym syms expr
#' @importFrom stringr str_detect regex str_extract str_remove
#' @importFrom dplyr filter mutate bind_cols rename
#' @importFrom tidyr unite
#' @importFrom purrr discard keep
#' @importFrom grDevices colours palette
#' @importFrom ggplot2 ggplot theme_bw labs theme element_text element_blank element_rect aes geom_segment arrow unit geom_point ggplot_build
#' @importFrom ggrepel geom_text_repel
#' @importFrom utils hasName
#' @importFrom tidyselect last_col
#' 
#'
#' @examples 
#' library(vegan)
#' library(tidyverse)
#' library(ggrepel)
#' 
#' data(dune)
#' data(dune.env)
#' 
#' m <- capscale(dune ~ A1 + Management, data = dune.env)
#' gordi_read(m, env = dune.env, scaling = 'species', correlation = T) |>
#'   gordi_species(label = F) |>
#'   gordi_predict(scaling_coefficient = 1)    
#' @export
gordi_predict <- function(
    pass,
    label = F,
    colour = '',
    alpha = '',
    arrow_size = '',
    linewidth = '',
    linetype = '',
    fill = '',
    shape = '',
    size = '',
    scaling_coefficient = 0.9,
    repel_label = F) {
  
  warning("So far, `gordi_predict()` can't calculate and plot scores for interactions of two continuous predictors and for interactions of continuous and categorical predictors.")
  
  if (!pass$type %in% c('RDA', 'CCA', 'db-RDA') || rlang::is_empty(pass$predictor_scores)) {
    stop("You provided an unconstrained ordination. Predictors can't be displayed. If you want to passively plot environmental variables, try `gordi_corr()`.")
  }
  
  ### ordination axis labels
  if(pass$type %in% c('DCA', 'NMDS')) {actual_labs <- paste0(pass$axis_names)} else 
  {actual_labs <- paste0(pass$axis_names, " (", round(pass$explained_variation[pass$choices]*100, 1), "%)")}
  
  
  
  
  ### --- MAIN TERMS ---
  # create a tibble with centroid scores
  if (length(factor_predictors) > 0) {
    
    factor_scores <- pass$predictor_scores |> 
      filter(score == 'centroids') |> 
      filter(str_detect(label, regex(paste(factor_predictors, collapse = '|')))) |> 
      mutate(predictor = str_extract(label, paste(factor_predictors, collapse = '|')),
             level = str_remove(label, paste(factor_predictors, collapse = '|')),
             ordered = as.character(grepl(paste(ordered_factors, collapse = '|'), label))) |> 
      rename(Axis_pred1 = 1,
             Axis_pred2 = 2) } else {factor_scores <- NULL}
  
  # create a tibble with arrow ends
  # Inside gordi_predict, after 'vector_predictors' is defined (around line 100)
  
  if (length(vector_predictors) > 0) {
    # If there ARE vector predictors, run your current logic:
    vector_scores <- pass$predictor_scores |>
      filter(score == "biplot") |>
      filter(str_detect(label, ":", negate = T)) |> 
      filter(str_detect(label, regex(paste(vector_predictors, collapse = '|')))) |> 
      bind_cols(vector_predictors) |> 
      rename('predictor' = last_col()) |> 
      mutate(level = NA_character_,
             ordered = NA_character_) |> # <--- Ensure NA_character_ here
      rename(Axis_pred1 = 1,
             Axis_pred2 = 2) 
  } else {vector_scores <- NULL}
  
  
  
  ### --- INTERACTION TERMS ---
  if (length(factor_predictors) > 1 && grepl(':|\\*', as.character(pass$m$call$formula)[3])) {
    
    # categorical variables in interaction
    # pre extract arguments
    m <- pass$m
    env <- pass$env
    choices <- pass$choices
    scaling <- pass$scaling
    correlation <- pass$correlation
    hill <- pass$hill
    const <- pass$const
    
    
    # extract centroids from LC site scores by defined interaction groups
    interaction_scores <- scores(m,
                                 display = 'all', 
                                 choices = choices, 
                                 scaling = scaling, 
                                 correlation = correlation, 
                                 hill = hill,
                                 const = const, 
                                 tidy = T) |> 
      as_tibble() |> 
      filter(score == 'constraints') |>
      bind_cols(env) |> 
      rename('Axis_pred1' = 1,
             'Axis_pred2' = 2) |> 
      select(starts_with('Axis'), matches(factor_predictors)) |> 
      group_by(across(where(is.factor) | where(is.character))) |>
      summarise(Axis_pred1 = mean(Axis_pred1),
                Axis_pred2 = mean(Axis_pred2)) |>
      ungroup() |> 
      relocate(where(is.numeric), .before = 1) |> 
      mutate(score = 'interaction_factors_centroids',
             label = paste(!!!syms(factor_predictors), sep = ":"),
             predictor_level = paste0(names(m$terminfo$xlev)[1],
                                      "_",
                                      !!sym(names(m$terminfo$xlev)[1]),
                                      ":", names(m$terminfo$xlev)[2],
                                      "_",
                                      !!sym(names(m$terminfo$xlev)[2])),
             predictor = paste(factor_predictors, collapse = ":"),
             level = paste(!!!syms(factor_predictors), sep = ":"),
             ordered = NA_character_) |> 
      select(-all_of(factor_predictors)) } else {interaction_scores <- NULL}
  
  
  pred_df <- bind_rows(vector_scores, factor_scores) |> 
    unite('predictor_level', predictor, level, sep = '_', remove = F, na.rm = T) |> 
    bind_rows(interaction_scores)
  
  pass$pred_df <- pred_df
  
  pass$predictor_scores <- pred_df
  
  # interaction terms???
  
  
  
  
  ### --- PLOT ---
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
  
  ### create factor_scores and vector_scores
  # containing only scores for vector predictors - to be plotted as arrows
  # and only scores for factor predictors - to be plotted as points
  
  # get names of categorical predictors
  factor_predictors <- names(pass$m$terminfo$xlev)
  
  # get names of continuous predictors
  vector_predictors <- names(pass$m$terminfo$ordered) |> 
    discard(~ .x %in% names(pass$m$terminfo$xlev))
  
  # get names of ordered categorical predictors (it might not be used, idk)
  ordered_factors <- pass$m$terminfo$ordered |> 
    keep(~ .x) |> 
    names()
  
  
  ### Detect mapped vs constant aesthetics
  
  # --- Shared Aesthetics (colour, alpha) ---
  # For 'vector_scores' (geom_segment)
  map_colour <- !identical(colour, '') && has_name(pred_df, colour)
  map_alpha <- !identical(alpha, '') && has_name(pred_df, alpha)
  
  # For 'factor_scores' (geom_point)
  map_colour <- !identical(colour, '') && has_name(pred_df, colour)
  map_alpha <- !identical(alpha, '') && has_name(pred_df, alpha)
  
  # Constant detection (only needs to be done once, as it checks the input value)
  is_valid_colour <- (grepl("^#(?:[A-Fa-f0-9]{6}|[A-Fa-f0-9]{3})$", colour) || colour %in% grDevices::colours()) || (is.character(colour) && colour %in% palette()) || (is.numeric(colour) && colour %in% seq_along(palette()))
  const_colour <- !identical(colour, '') && !map_colour && !map_colour && is_valid_colour
  const_alpha <- !identical(alpha, '') && !map_alpha && !map_alpha && is.numeric(alpha)
  
  
  # --- Segment-Only Aesthetics (linetype, linewidth) ---
  # Check mapping against 'vector_scores'
  map_linetype <- !identical(linetype, '') && has_name(pred_df, linetype)
  map_linewidth <- !identical(linewidth, '') && has_name(pred_df, linewidth)
  
  # Constant detection
  const_linetype <- !identical(linetype, '') && !map_linetype && (is.numeric(linetype) || is.character(linetype))
  const_linewidth <- !identical(linewidth, '') && !map_linewidth && is.numeric(linewidth)
  
  
  # --- Point-Only Aesthetics (fill, shape, size) ---
  # Check mapping against 'factor_scores'
  map_fill <- !identical(fill, '') && has_name(pred_df, fill)
  map_shape <- !identical(shape, '') && has_name(pred_df, shape)
  map_size <- !identical(size, '') && has_name(pred_df, size)
  
  # Constant detection
  is_valid_fill <- (grepl("^#(?:[A-Fa-f0-9]{6}|[A-Fa-f0-9]{3})$", fill) || fill %in% grDevices::colours()) || (is.character(fill) && fill %in% palette()) || (is.numeric(fill) && fill %in% seq_along(palette()))
  const_fill <- !identical(fill, '') && !map_fill && is_valid_fill
  const_shape <- !identical(shape, '') && !map_shape && (is.numeric(shape) || is.character(shape))
  const_size <- !identical(size, '') && !map_size && is.numeric(size)
  
  
  
  ### Set scaling coefficient
  # extract plot frame size (x and y axis lengths)
  if (!is.null(vector_scores)) {
    p_build <- ggplot_build(p)
    
    plot_range <- c(xmin_plot = p_build$layout$panel_params[[1]]$x.range[1],
                    xmax_plot = p_build$layout$panel_params[[1]]$x.range[2],
                    ymin_plot = p_build$layout$panel_params[[1]]$y.range[1],
                    ymax_plot = p_build$layout$panel_params[[1]]$y.range[2])
    
    predictor_range <- c(xmin_pred = min(pred_df[pred_df$score == 'biplot',1]),
                         xmax_pred = max(pred_df[pred_df$score == 'biplot',1]),
                         ymin_pred = min(pred_df[pred_df$score == 'biplot',2]),
                         ymax_pred = max(pred_df[pred_df$score == 'biplot',2]))
    
    coef <- (max(abs(plot_range)) / max(abs(predictor_range))) * scaling_coefficient
  } else {coef <- 1}
  
  
  ### Prepare aes arguments for geom_segment()
  aes_args_segment <- list(
    x = 0, y = 0,
    xend = expr(Axis_pred1 * !!coef),
    yend = expr(Axis_pred2 * !!coef)
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
  p <- p + do.call(geom_segment,
                   c(list(data = pred_df |> filter(score == 'biplot'),
                          mapping = do.call(aes, aes_args_segment)),
                     const_args_segment))
  
  
  if (isTRUE(label)){
    if (isTRUE(repel_label)){
      p <- p + ggrepel::geom_text_repel(data = pred_df |> filter(score == 'biplot'), aes(Axis_pred1 * coef, Axis_pred2 * coef, label = predictor_level), colour = colour)
    } else {
      p <- p + geom_text(data = pred_df |> filter(score == 'biplot'), aes(Axis_pred1 * coef, Axis_pred2 * coef, label = predictor_level), colour = colour)
    }
  }
  
  
  
  ### Prepare aes arguments for geom_point()
  # Start with fixed x/y for the base (0,0) and end at the species scores
  aes_args_point <- list(
    x = sym("Axis_pred1"),
    y = sym("Axis_pred2")
  )
  
  if(map_colour) aes_args_point$colour <- sym(colour)
  if(map_alpha) aes_args_point$alpha <- sym(alpha)
  if(map_fill) aes_args_point$fill <- sym(fill)
  if(map_shape) aes_args_point$shape <- sym(shape)
  if(map_size) aes_args_point$size <- sym(size)
  
  
  ### Prepare constant arguments for geom_point() (mapped first, then defaults if nothing)
  const_args_point <- list()
  
  # Add constant arguments for geom_point() if not mapped
  # colour
  if(!map_colour){
    if(!identical(colour, '')) {const_args_point$colour <- colour} else {const_args_point$colour <- 2}}
  # alpha
  if(!map_alpha){
    if(!identical(alpha, '')) {const_args_point$alpha <- alpha} else {const_args_point$alpha <- 1}}
  # linetype
  if(!map_fill){
    if(!identical(fill, '')) {const_args_point$fill <- fill} else {const_args_point$fill <- 3}}
  # linewidth
  if(!map_shape){
    if(!identical(shape, '')) {const_args_point$shape <- shape} else {const_args_point$shape <- 16}}
  # linewidth
  if(!map_size){
    if(!identical(size, '')) {const_args_point$size <- size} else {const_args_point$size <- 4}}
  
  
  ### add to plot
  p <- p + do.call(geom_point,
                   c(list(data = pred_df |> filter(score %in% c('centroids', 'interaction_factors_centroids')),
                          mapping = do.call(aes, aes_args_point)),
                     const_args_point))
  
  if (isTRUE(label)){
    if (isTRUE(repel_label)){
      p <- p + ggrepel::geom_text_repel(data = pred_df |> filter(score %in% c('centroids', 'interaction_factors_centroids')), aes(Axis_pred1, Axis_pred2, label = predictor_level), colour = 2)
    } else {
      p <- p + geom_text(data = pred_df |> filter(score == 'centroids'), aes(Axis_pred1, Axis_pred2, label = predictor_level), colour = 2)
    }
  }
  
  
  
  ### save plot
  pass$plot <- p
  
  
  return(pass)
}

