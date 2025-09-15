#' Draw ordination sites. 
#' 
#' @description
#' The function [gordi_sites()] takes the result of [gordi_read()] and adds sites (as points) 
#' and labels (optional) to a ggplot object created by the function itself. 
#' In addition to drawing sites you can change a wide range of parameters.
#' In the process it renames ordination axes to Axis_site1 and Axis_site2, while the original vegan names 
#' are passed to 'actual_labs' used for plotting.
#' 
#' @details
#' The function builds a dataframe 'site_df' by binding 'pass$env' and 'pass$site_scores'.
#' Aesthetics can be mapped by passing a column name present in the site_df dataframe or by using literal values (e.g. hex colour, numeric size).
#' Fresh colour/fill scales are started for the site layer by ggnewscale::new_scale_colour() and/or ggnewscale::new_scale_fill(), so that later calls (e.g. gordi_colour, layers addded after sites) can define their own independent colour/fill scales.
#' Additional fresh scales (for size, alpha, stroke, shape) are created at the end of the function, so subsequent layers do not inherit site mappings.
#' If label = TRUE, site labels are drawn using the first column in site_df dataframe.
#' The default setting for repel_label is TRUE, so the labels are drawn using `geom_text_repel()`. 
#' For changing colour and label aesthetics please refer to [gordi_colour()] and [gordi_label()].
#' 
#' @param pass A list object produced by [gordi_read()]
#' @param label Logical; Default is `TRUE`. Draws site labels (`TRUE`) or not (`FALSE`). Site labels can use only the first column of 'site_df', if you want custom labels, please use function [gordi_label()], which overrides label layer settings. 
#' @param fill Changes site fill colour. Can be either: a variable = a specific column name of 'site_df', a constant value: R colour name, numeric value or hex (e.g. `'green'`, `1`, `'#6aa84f'`). When left as `''`, the internal default colours are used.
#' @param colour Changes site colours. Can be either: a variable = a specific column name of 'site_df', a constant value: R colour name, numeric value or hex (e.g. `'green'`, `1`, `'#6aa84f'`). When left as `''`, the internal default colours are used.
#' @param alpha Integer; Changes transparency of points. When left as `''`, the internal default values are used (alpha = 1).
#' @param stroke Integer; Change shape or arrow width. Defined as numeric values. When left as `''`, the internal default values are used (stroke = 0.5).
#' @param shape Integer; Change the shape of points. Either by a specific mapping column from site_df or by a constant (numeric value). When left as `''`, the internal default values are used (shape = 16).
#' @param size Integer; Change the size of points. Either by a specific mapping column from site_df or by a constant (numeric value). When left as `''`, the internal default values are used (size = 3).
#' @param repel_label Logical; The labels are drawn by [geom_text_repel()]. The default is set to repel_label = TRUE.
#' 
#' @return The updated 'pass' object with site layers appended to `pass$plot`.
#' 
#' @examples
#' o <- gordi_read(m, dune.env)
#' 
#' # basic sites plot with default styling and labels
#' o |> gordi_sites()
#' 
#' # mapped colour and size from env dataset, constant alpha, shape, fill
#' o |> gordi_sites(colour = 'elevation', size = 'elevation', alpha = 0.7, shape = 21, fill = 'white')
#' 
#' # no labels, constant green coloured points
#' o |> gordi_sites(label = FALSE, colour = 'green')
#' 
#' @seealso [gordi_read()], [gordi_species()], [gordi_label()], [gordi_colour()], [gordi_predict()], [ggplot2::ggplot()], [ggrepel::geom_text_repel()]
#' @export
gordi_sites <- function(pass,
                        label = TRUE,
                        fill = '',
                        alpha = '',
                        stroke = '',
                        shape = '',
                        size = '',
                        colour = '',
                        repel_label = T) {
  
  # axis names
  names(pass$site_scores) <- paste0("Axis_site", 1:2)
  
  # actual labs
  if(pass$type %in% c('DCA', 'NMDS')) {actual_labs <- paste0(pass$axis_names)} else 
  {actual_labs <- paste0(pass$axis_names, " (", round(pass$explained_variation[pass$choices]*100, 2), "%)")}
  
  # plot set up  
  if (is.null(pass$plot)) { 
    p <- ggplot() +
      theme_bw() +
      labs(x = actual_labs[1], y = actual_labs[2]) +
      theme(
        text = element_text(size = 15),
        panel.grid = element_blank(),
        legend.justification = c(1, 1)
      )
  } else {
    p <- pass$plot
  }
  
  p <- p + ggnewscale::new_scale_colour()
  p <- p + ggnewscale::new_scale_fill()
  
  site_df <- bind_cols(pass$env, pass$site_scores)
  
  # Detect mapped vs constant aesthetics
  map_colour <- !identical(colour, '') && has_name(site_df, colour)
  const_colour <- const_colour <- !map_colour && (grepl("^#(?:[A-Fa-f0-9]{6}[A-Fa-f0-9]{3})$", colour) || colour %in% grDevices::colours()) || (is.character(colour) && colour %in% palette()) || (is.numeric(colour) && colour %in% seq_along(palette()))
  
  if(!identical(colour, '')){
    message("To customize colours (similarly to `ggplot2::scale_colour_()` functions), please use gordi_colour() right after `gordi_sites()`.")
    if(!map_colour && !const_colour){
      warning("`colour` must be either a column in the `env` dataframe, a valid R colour name/hex code, or a numeric code! Ignoring input, default is being used.")
      colour <- ''
    }
  }
  
  map_size <- !identical(size, '') && has_name(site_df, size)
  const_size <- !map_size && is.numeric(size)
  
  if(!identical(size, '')){
    if(!map_size && !const_size){
      warning("`size` must be either a numeric constant or a numeric column in the `env` dataframe! Ignoring input, default is being used.")
      size <- ''
    }
  }
  
  map_shape <- !identical(shape, '') && has_name(site_df, shape)
  const_shape <- !map_shape && is.numeric(shape)
  
  if(!identical(shape, '')){
    message("To customize shapes (similarly to `ggplot2::scale_shape_()` functions), please use gordi_shape() right after `gordi_sites()`.")
    if(!map_shape && !const_shape){
      warning("`shape` must be either a numeric constant (0-25) or a numeric column in the `env` dataframe! Ignoring input, default is being used.")
      shape <- ''
    }
  }
  
  map_fill <- !identical(fill, '') && has_name(site_df, fill)
  const_fill <- !map_fill && (grepl("^#(?:[A-Fa-f0-9]{6}[A-Fa-f0-9]{3})$", fill) || fill %in% grDevices::colours())|| (is.character(fill) && fill %in% palette()) || (is.numeric(fill) && fill %in% seq_along(palette()))
  
  if(!identical(fill, '')){
    message("To customize fill colours (similarly to `ggplot2::scale_fill_()` functions), please use gordi_colour() right after `gordi_sites()`.")
    if(!map_fill && !const_fill){
      warning("`fill` must be either a column in the `env` dataframe, a valid R colour name/hex code, or a numeric code! Ignoring input, default is being used.")
      fill <- ''
    }
  }
  
  map_alpha <- !identical(alpha, '') && has_name(site_df, alpha)
  const_alpha <- !map_alpha && is.numeric(alpha)
  
  if(!identical(alpha, '')){
    if(!map_alpha && !const_alpha){
      warning("`alpha` must be either a numeric constant (0-1) or a numeric column in the `env` dataframe! Ignoring input, default is being used.")
      alpha <- ''
    }
  }
  
  map_stroke <- !identical(stroke, '') && has_name(site_df, stroke)
  const_stroke <- !map_stroke && is.numeric(stroke)
  
  if(!identical(stroke, '')){
    if(!map_stroke && !const_stroke){
      warning("`stroke` must be either a numeric constant or a numeric column in the `env` dataframe! Ignoring input, default is being used.")
      stroke <- ''
    }
  }
  
  # Prepare aes arguments (only mapped)
  aes_args_point <- list(x = quote(Axis_site1), y = quote(Axis_site2))
  if(map_colour) aes_args_point$colour <- sym(colour)
  if(map_size) aes_args_point$size <- sym(size)
  if(map_shape) aes_args_point$shape <- sym(shape)
  if(map_fill) aes_args_point$fill <- sym(fill)
  if(map_alpha) aes_args_point$alpha <- sym(alpha)
  if(map_stroke) aes_args_point$stroke <- sym(stroke)
  
  # Prepare constant arguments (mapped first, then defaults if nothing)
  const_args_point<- list()
  if(!map_colour) {
    if(!identical(colour, '')) { 
      const_args_point$colour <- colour} else {const_args_point$colour <- 1}}
  
  if(!map_size) {
    if(!identical(size, '')) {
      const_args_point$size <- size} else {const_args_point$size <- 3}}
  
  if(!map_shape) {
    if(!identical(shape, '')) { 
      const_args_point$shape <- shape} else {const_args_point$shape <- 16}}
  
  if(!map_fill) {
    if(!identical(fill, '')) { 
      const_args_point$fill <- fill} else {const_args_point$fill <- 'white'}}
  
  if(!map_alpha) {
    if(!identical(alpha, '')) { 
      const_args_point$alpha <- alpha} else {const_args_point$alpha <- 1}}
  
  if(!map_stroke){
    if(!identical(stroke, '')){
      const_args_point$stroke <- stroke} else {const_args_point$stroke <- 0.5}}  
  
  # plot  
  p <- p + do.call(geom_point, c(list(mapping = do.call(aes, aes_args_point), data = site_df), const_args_point))
 
  if (isTRUE(label)){
    message("Labels have been drawn. To customize labels, please use `gordi_label()` right after `gordi_sites()`.")
    labcol <- names(site_df)[1]
    if (isTRUE(repel_label)){
      p <- p + geom_text_repel(data = site_df, aes(Axis_site1, Axis_site2, label = !!sym(labcol)), colour = 'black') 
    } else {
      p <- p + geom_text(data = site_df, aes(Axis_site1, Axis_site2, label = !!sym(labcol)), colour = 'black')
    }
  } else if(is.character(label)){
      warning("`label` must be logical (TRUE/FALSE). To customize labels please use `gordi_label()` right after `gordi_sites()`. Ignoring `label` input, setting `label = FALSE`.")
    label <- FALSE
    }
   
  p <- p + ggnewscale::new_scale("size") 
  p <- p + ggnewscale::new_scale("shape")  
  p <- p + ggnewscale::new_scale("alpha")
  p <- p + ggnewscale::new_scale("stroke") 
  
  pass$plot <- p
  
  return(pass)
}

