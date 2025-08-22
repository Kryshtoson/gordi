# gordi_predict()
library(vegan)
library(tidyverse)
library(ggrepel)


gordi_predict <- function(
    pass,
    label = '',
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
  #' if colour != "" AND ALSO colour represents a colname present in spe_df, then map_colour is TRUE, otherwise is FALSE
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
                     const_args_segment)) +
    geom_text_repel(data = pred_df,
                    aes(x = Axis_pred1 * coef,
                        y = Axis_pred2 * coef,
                        label = predictor_names), colour = 2)
  
  ### save plot
  pass$plot <- p
  
  
  return(pass)
}