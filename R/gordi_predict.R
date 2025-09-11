#' Extracts predictor scores and relevant stuff and creates plot
#' 
#' @description
#' [gordi_predict()] takes the result of [gordi_read()] and creates plot with
#' predictor arrows. So far, it can work with continuous explanatory variables,
#' as for categorical we need to come up with a way how to calculate scores
#' for the last category. It also can't display explanatory variables passively.
#' 
#' There will be a special function for this once, maybe something as
#' `gordi_passive_agressive()`. Similarly to [gordi_species()] and [gordi_sites()],
#' you can also set a wide range of graphing parameters, such as colour, fill, size,
#' shape, alpha, stroke, and more traditional ggplot arguments, which can read
#' both, static and dynamic variable (e.g., 'red' or 'elevation').
#' 
#' @param pass Object from [gordi_read()] function.
#' @param label Logical; default = T, whether to display label by each point/arrow 
#'   or not. In [gordi_predict()], label can be only `predictor_names`, which are 
#'   displayed as a full name. If you want to customize the labels, you can it
#'   with [gordi_label()] function, which overrides this setting.
#' @param colour Colour can be defined statically as word from the [colours()]
#'   list (e.g. 'red'), HEX code (e.g. #5d782e), or number from [palette()] (e.g. 3).
#'   It can also defined dynamically (according to some variable, e.g. elevation
#'   or vegetation type). This variable has to be present in the `env` table and its
#'   name has to be written in "quotation marks". Default colour of arrows is 4.
#' @param alpha Transparency of symbols. Numeric value. Default is 0.6. Can be set
#'   statically or dynamically (use "").      
#' @param arrow_size Numeric; length of arrow in cm. Default is 0.3 cm. 
#' @param linewidth Self-explanatory. Can be numerical value. See
#'   [aes_linetype_size_shape()] for more details.
#' @param linetype Changes type of line used in arrows. Can be specified with either
#'   an integer (0-6), a name (0 = blank, 1 = solid, 2 = dashed, 3 = dotted,
#'   4 = dotdash, 5 = longdash, 6 = twodash). See [aes_linetype_size_shape()] for 
#'   more details.
#' @param scaling_coefficient Numeric; edits the lengts of predictor arrow. Default value
#'   is `0.9` which means that the predictor arrow will be as long as 0.9 of
#'   the longest axis displayed in the plotframe.
#' @param repel_label Logical; repels labels of species for better readability.
#'   Default is F. If you want to customize the labels, you can do it with
#'   [gordi_label()] function, which overrides this setting.   
#' 
#' 
#' 
#' @return The input object `pass` with an updated `plot` element that includes
#'   predictor arrows and labels.
#' 
#' @examples 
#' library(vegan)
#' library(tidyverse)
#' library(ggrepel)
#' 
#' data(dune)
#' data(dune.env)
#' 
#' m <- capscale(dune ~ A1, data = dune.env)
#' gordi_read(m, env = dune.env, scaling = 'species', correlation = T) |>
#'   gordi_species(label = F) |>
#'   gordi_predict(scaling_coefficient = 1)   
#' @export
gordi_predict <- function(
    pass,
    label = T,
    colour = '',
    alpha = '',
    arrow_size = '',
    linewidth = '',
    linetype = '',
    scaling_coefficient = 0.9,
    repel_label = T) {
  
  # Check if there are any predictors
  if (rlang::is_empty(pass$predictor_scores)) {
    stop("You did not use any predictor. Add some before calling gordi_predict().")
  }
  
  ### ordination axis labels
  if(pass$type %in% c('DCA', 'NMDS')) {actual_labs <- paste0(pass$axis_names)} else 
  {actual_labs <- paste0(pass$axis_names, " (", round(pass$explained_variation[pass$choices]*100, 2), "%)")}
  
  
  ### axis names used in spe_df
  names(pass$predictor_scores) <- paste0("Axis_pred", 1:2)
  
  
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
  
  ### create spe_df which is then called in the ggplot (spe_df exists only in this function and does not pass to the next)
  pred_df <- bind_cols(pass$predictor_scores, pass$predictor_names)
  
  
  ### Detect mapped vs constant aesthetics
  # colour
  # if colour != "" AND ALSO colour represents a colname present in spe_df, then map_colour is TRUE, otherwise is FALSE
  map_colour <- !identical(colour, '') && has_name(pred_df, colour) 
  # if map_colour is FALSE AND ALSO the thing inputed in arguments is a HEX code or is included in colours() or in palette() (word or number), then use it as const_colour
  const_colour <- !map_colour && (grepl("^#(?:[A-Fa-f0-9]{6}[A-Fa-f0-9]{3})$", colour) || colour %in% grDevices::colours()) || (is.character(colour) && colour %in% palette()) || (is.numeric(colour) && colour %in% seq_along(palette()))
  # alpha
  map_alpha <- !identical(alpha, '') && has_name(pred_df, alpha)
  const_alpha <- !map_alpha && is.numeric(alpha)
  # linetype
  map_linetype <- !identical(linetype, '') && has_name(pred_df, linetype)
  const_linetype <- !map_linetype && (is.numeric(linetype) || !identical(linetype,''))
  # linewidth
  map_linewidth <- !identical(linewidth, '') && has_name(pred_df, linewidth)
  const_linewidth <- !map_linewidth && is.numeric(linewidth)
  
  
  ### Set scaling coefficient
  # extract plot frame size (x and y axis lengths)
  p_build <- ggplot_build(p)
  
  plot_range <- c(xmin_plot = p_build$layout$panel_params[[1]]$x.range[1],
                  xmax_plot = p_build$layout$panel_params[[1]]$x.range[2],
                  ymin_plot = p_build$layout$panel_params[[1]]$y.range[1],
                  ymax_plot = p_build$layout$panel_params[[1]]$y.range[2])
  
  predictor_range <- c(xmin_pred = min(pred_df[,1]),
                       xmax_pred = max(pred_df[,1]),
                       ymin_pred = min(pred_df[,2]),
                       ymax_pred = max(pred_df[,2]))
  
  coef <- (max(abs(plot_range)) / max(abs(predictor_range))) * scaling_coefficient
  
  
  ### Prepare aes arguments for geom_segment()
  # Start with fixed x/y for the base (0,0) and end at the species scores
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
                   c(list(data = pred_df,
                          mapping = do.call(aes, aes_args_segment)),
                     const_args_segment)) 
  
  if (isTRUE(label)){
    labcol <- names(pred_df)[3]
    if (isTRUE(repel_label)){
      p <- p + geom_text_repel(data = pred_df, aes(Axis_pred1*coef, Axis_pred2*coef, label = !!sym(labcol)), colour = 'black') 
    } else {
      p <- p + geom_text(data = pred_df, aes(Axis_pred1*coef, Axis_pred2*coef, label = !!sym(labcol)), colour = 'black')
    }
  }
  
  ### save plot
  pass$plot <- p
  
  
  return(pass)
}

