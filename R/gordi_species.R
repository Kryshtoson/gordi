#' Extracts species scores and relevant stuff and creates plot
#' 
#' @description
#' [gordi_species()] takes the result of [gordi_read()] and creates plot with
#' species scores. PCA, RDA, PCoA, and db-RDA use arrows by default, and CA, CCA,
#' DCA and NMDS use points. This can be changed manually. You can also set a wide
#' range of graphing parameters, such as colour, fill, size, shape, alpha, stroke,
#' and more traditional ggplot arguments, which can read both, static and dynamic
#' variable (e.g., 'red' or 'elevation').
#' 
#' In the process it renames ordination
#' axes to Axis_spe1 and Axis_spe2, and the original vegan names passes
#' to "actual_labs" which are used just for plotting.
#' 
#' @param pass Object from [gordi_read()] function.
#' @param label Logical; default = TRUE, whether to display label by each point/arrow or not. In [gordi_species()], label can be only `species_names`, which are displayed as a full name. If you want to customize the labels, use [gordi_label()] which overrides this setting.
#' @param symbol `c('default', 'point', 'arrow')`; How should species scores be 
#'   displayed. 
#' @param colour Colour can be defined statically as word from the [colours()]
#'   list (e.g. 'red'), HEX code (e.g. #5d782e), or number from [palette()] (e.g. 3).
#'   It can also defined dynamically (according to some variable, e.g. elevation
#'   or vegetation type). This variable has to be present in the `env` table and its
#'   name has to be written in "quotation marks". Default colour of arrows is 4.
#' @param size Changes size of points. Can be defined as a number (statically) or 
#'   by a variable (dynamically) - the name of the variable has to be in "quotation
#'   marks". Default size is 3.
#' @param shape Defines shape of points. Can be statical (numeric value from 0 to 25)
#'   or dynamical (by categorical variable written in "quotation marks"). Default
#'   shape is 16. Symbols 0-20 have only colour, symbols 21-25 have both, colour 
#'   and fill, which can be defined separately. 
#' @param fill Fill can be defined statically as word from the [colours()]
#'   list (e.g. 'red'), HEX code (e.g. #5d782e), or number from [palette()] (e.g. 3).
#'   It can also defined dynamically (according to some variable, e.g. elevation
#'   or vegetation type). This variable has to be present in the `env` table and
#'   its name has to be written in "quotation marks". Default fill "white".
#' @param alpha Transparency of symbols. Numeric value. Default is 0.6. Can be set
#'   statically or dynamically (use "").
#' @param stroke Shape or arrow outline width. Numeric value. Can be set statically
#'   or dynamically (use "").
#' @param linetype Changes type of line used in arrows. Can be specified with either
#'   an integer (0-6), a name (0 = blank, 1 = solid, 2 = dashed, 3 = dotted,
#'   4 = dotdash, 5 = longdash, 6 = twodash). See [aes_linetype_size_shape()] for 
#'   more details.
#' @param linewidth Self-explanatory. Can be numerical value. See
#'   [aes_linetype_size_shape()] for more details.
#' @param arrow_size Changes the size of arrow. Numeric value. Default is 0.3 cm.
#'   Can't be set dynamically. Just can't. Why would you do that.
#' @param repel_label Logical; repels labels of species for better readability.
#'   Default is F. If you want to customize the labels, you can it with
#'   [gordi_label()] function, which overrides this setting.
#' 
#' 
#' @return The input `pass` object with an updated element:
#'   \describe{
#'     \item{plot}{A `ggplot` object containing the ordination plot with
#'     species scores added.}
#'   }
#'
#' @seealso [gordi_sites()], [ggplot2::ggplot()], [ggrepel::geom_text_repel()], 
#'   [gordi_read()], [gordi_predict()] 

#' @examples 
#' library(vegan)
#' library(tidyverse)
#' library(ggrepel)
#' 
#' data(dune)
#' m <- capscale(dune ~ 1)
#' o <- gordi_read(m) |> gordi_species()
#' o
#' @export
gordi_species <- function(pass,
                          label = FALSE,
                          symbol = c('default', 'point', 'arrow'),
                          colour = '',
                          size = '',
                          shape = '',
                          fill = '',
                          alpha = '',
                          stroke = '',
                          linetype = '',
                          linewidth = '',
                          arrow_size = '',
                          repel_label = FALSE) {
  
  ### axis names used in spe_df 
  names(pass$species_scores) <- paste0("Axis_spe", 1:2)
  
  
  ### ordination types -> later used in axis labels 
  if(pass$type %in% c('DCA', 'NMDS')) {actual_labs <- paste0(pass$axis_names)} else 
  {actual_labs <- paste0(pass$axis_names, " (", round(pass$explained_variation[pass$choices]*100, 2), "%)")}
  
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
  spe_df <- bind_cols(pass$species_names, pass$species_scores)
  
  # joins traits - only one trait value for one plant in one column is permitted
  if (!is.null(pass$traits)) {
    spe_df <- spe_df |>
      left_join(pass$traits,
                by = join_by(!!sym(names(spe_df)[1]) == !!sym(names(pass$traits)[1])))
  }
  
  p <- p + ggnewscale::new_scale_colour() 
  p <- p + ggnewscale::new_scale_fill()
  
  ### Detect mapped vs constant aesthetics
  
  ### arguments working in both, points and arrows (colour, alpha)
  # colour
  # if colour != "" AND ALSO colour represents a colname present in spe_df, then map_colour is TRUE, otherwise is FALSE
  map_colour <- !identical(colour, '') && has_name(spe_df, colour) 
  # if map_colour is FALSE AND ALSO the thing inputed in arguments is a HEX code or is included in colours() or in palette() (word or number), then use it as const_colour
  const_colour <- !map_colour && (grepl("^#(?:[A-Fa-f0-9]{6}[A-Fa-f0-9]{3})$", colour) || colour %in% grDevices::colours()) || (is.character(colour) && colour %in% palette()) || (is.numeric(colour) && colour %in% seq_along(palette()))
 
  if(!identical(colour, '')){
    message("To customize colours (similarly to `ggplot2::scale_colour_()` functions), please use gordi_colour() right after `gordi_species()`.")
    if(!map_colour && !const_colour){
      warning("`colour` must be either a column in the `env` dataframe, a valid R colour name/hex code, or a numeric code! Ignoring input, default is being used.")
      colour <- ''
    }
  }
  
   # alpha
  map_alpha <- !identical(alpha, '') && has_name(spe_df, alpha)
  const_alpha <- !map_alpha && is.numeric(alpha)
  
  if(!identical(alpha, '')){
    if(!map_alpha && !const_alpha){
      warning("`alpha` must be either a numeric constant (0-1) or a numeric column in the `env` dataframe! Ignoring input, default is being used.")
      alpha <- ''
    }
  }
  
  ### arguments working only in geom_point
  # size
  map_size <- !identical(size, '') && has_name(spe_df, size)
  const_size <- !map_size && is.numeric(size)
  
  if(!identical(size, '')){
    if(!map_size && !const_size){
      warning("`size` must be either a numeric constant or a numeric column in the `env` dataframe! Ignoring input, default is being used.")
      size <- ''
    }
  }
  
  # shape
  map_shape <- !identical(shape, '') && has_name(spe_df, shape)
  const_shape <- !map_shape && is.numeric(shape)
  
  if(!identical(shape, '')){
    message("To customize shapes (similarly to `ggplot2::scale_shape_()` functions), please use gordi_shape() right after `gordi_species()`.")
    if(!map_shape && !const_shape){
      warning("`shape` must be either a numeric constant (0-25) or a numeric column in the `env` dataframe! Ignoring input, default is being used.")
      shape <- ''
    }
  }
  
  # fill
  map_fill <- !identical(fill, '') && has_name(spe_df, fill)
  const_fill <- !map_fill && (grepl("^#(?:[A-Fa-f0-9]{6}[A-Fa-f0-9]{3})$", fill) || fill %in% grDevices::colours())
 
  if(!identical(fill, '')){
    message("To customize fill colours (similarly to `ggplot2::scale_fill_()` functions), please use gordi_colour() right after `gordi_species()`.")
    if(!map_fill && !const_fill){
      warning("`fill` must be either a column in the `env` dataframe, a valid R colour name/hex code, or a numeric code! Ignoring input, default is being used.")
      fill <- ''
    }
  }
  
   # stroke
  map_stroke <- !identical(stroke, '') && has_name(spe_df, stroke)
  const_stroke <- !map_stroke && is.numeric(stroke)
  
  if(!identical(stroke, '')){
    if(!map_stroke && !const_stroke){
      warning("`stroke` must be either a numeric constant or a numeric column in the `env` dataframe! Ignoring input, default is being used.")
      stroke <- ''
    }
  }
  
  ### arguments working only in geom_segment
  # linetype
  line_types <- c('solid', 'dashed', 'dotted', 'dotdash', 'longdash', 'twodash')
  map_linetype <- !identical(linetype,'') && has_name(spe_df, linetype)
  const_linetype <- !map_linetype && (is.numeric(linetype) || !identical(linetype,'')) && linetype %in% line_types
  
  if(!identical(linetype, '')){
    if(!map_linetype && !const_linetype){
      warning("`linetype` must be either a numeric constant, a column name or a valid linetype (e.g. 'solid', 'dashed', 'dotted', 'dotdash', 'longdash', 'twodash'. Ignoring input, default is being used")
      linetype <- ''
    }
  }
  
  # linewidth
  map_linewidth <- !identical(linewidth, '') && has_name(spe_df, linewidth)
  const_linewidth <- !map_linewidth && is.numeric(linewidth)
  
  if(!identical(linewidth, '')){
    if(!map_linewidth && !const_linewidth){
      warning("`linewidth` must ve either a numeric constant or a numeric column. Ignoring input, default is being used.")
      linewidth <- ''
      }
  }
  
  
  ### Prepare aes arguments for geom_point()
  aes_args_point <- list(x = quote(Axis_spe1), y = quote(Axis_spe2))
  
  if(map_colour) aes_args_point$colour <- sym(colour)
  if(map_size) aes_args_point$size <- sym(size)
  if(map_shape) aes_args_point$shape <- sym(shape)
  if(map_fill) aes_args_point$fill <- sym(fill)
  if(map_alpha) aes_args_point$alpha <- sym(alpha)
  if(map_stroke) aes_args_point$stroke <- sym(stroke)
  
  
  ### Prepare constant arguments for geom_point() (mapped first, then defaults if nothing)
  const_args_point <- list()
  
  # colour 
  if(!map_colour){
    if(!identical(colour, '')) {const_args_point$colour <- colour} else {const_args_point$colour <- 4}}
  # alpha 
  if(!map_alpha){
    if(!identical(alpha, '')) {const_args_point$alpha <- alpha} else {const_args_point$alpha <- 1}}
  # size 
  if(!map_size){
    if(!identical(size, '')) {const_args_point$size <- size} else {const_args_point$size <- 3}}
  # shape 
  if(!map_shape){
    if(!identical(shape, '')) {const_args_point$shape <- shape} else {const_args_point$shape <- 16}}
  # fill 
  if(!map_fill){
    if(!identical(fill, '')) {const_args_point$fill <- fill} else {const_args_point$fill <- "white"}}
  # stroke
  if(!map_stroke){
    if(!identical(stroke, '')){const_args_point$stroke <- stroke} else {const_args_point$stroke <- 0.5}}
  
  
  ### Prepare aes arguments for geom_segment()
  # Start with fixed x/y for the base (0,0) and end at the species scores
  aes_args_segment <- list(
    x = 0, y = 0,
    xend = quote(Axis_spe1),
    yend = quote(Axis_spe2)
  )
  
  if(map_colour) aes_args_segment$colour <- sym(colour)
  if(map_linewidth) aes_args_segment$linewidth <- sym(linewidth)
  if(map_linetype) aes_args_segment$linetype <- sym(linetype)
  if(map_alpha) aes_args_segment$alpha <- sym(alpha)
  
  
  ### Prepare constant arguments for geom_point() (mapped first, then defaults if nothing)
  const_args_segment <- list()
  
  # Add constant arguments for geom_segment() if not mapped
  # colour
  if(!map_colour){
    if(!identical(colour, '')) {const_args_segment$colour <- colour} else {const_args_segment$colour <- 4}}
  # alpha
  if(!map_alpha){
    if(!identical(alpha, '')) {const_args_segment$alpha <- alpha} else {const_args_segment$alpha <- 0.6}}
  # linetype
  if(!map_linetype){
    if(!identical(linetype, '')) {const_args_segment$linetype <- linetype} else {const_args_segment$linetype <- 1}}
  # linewidth
  if(!map_linewidth){
    if(!identical(linewidth, '')) {const_args_segment$linewidth <- linewidth} else {const_args_segment$linewidth <- 0.5}}
  # arrow_size (does not make sense to have map_arrow)
  const_args_segment$arrow <- arrow(
    length = unit(
      if (!identical(arrow_size, '')) as.numeric(arrow_size) else 0.3, "cm"))
  
  
  
  # Add the layer
  # If linear ordination is used (PCA, RDA, PCoA, db-RDA), arrows are used
  # if unimodal (CA, CCA, DCA, NMDS), points are used
  
  if (is.null(symbol) || any(symbol == "default")) {
    if (pass$type %in% c('CA', 'CCA', 'DCA', 'NMDS')) {
      p <- p + do.call(geom_point, c(list(mapping = do.call(aes, aes_args_point), data = spe_df), const_args_point))
    } else {
      p <- p + do.call(geom_segment, c(list(data = spe_df, mapping = do.call(aes, aes_args_segment)), const_args_segment))
    }
  } else if (symbol == 'point') {
    p <- p + do.call(geom_point, c(list(mapping = do.call(aes, aes_args_point), data = spe_df), const_args_point))
  } else if (symbol == 'arrow') {
    p <- p + do.call(geom_segment, c(list(data = spe_df, mapping = do.call(aes, aes_args_segment)), const_args_segment))
  }

  # 
  if (isTRUE(label)){
    message("Labels have been drawn. To customize labels, please use `gordi_label()` right after `gordi_sites()`.")
    if (isTRUE(repel_label)){
      p <- p + geom_text_repel(data = spe_df, aes(Axis_spe1, Axis_spe2, label = species_names), colour = 'black') 
    } else {
      p <- p + geom_text(data = spe_df, aes(Axis_spe1, Axis_spe2, label = species_names), colour = 'black')
    }
  } else if(is.character(label)){
    warning("`label` must be logical (TRUE/FALSE). To customize labels please use `gordi_label()` right after `gordi_sites()`. Ignoring `label` input, setting `label = FALSE`.")
    label <- FALSE
  }
  
  
  # More scales possibility (e.g. one colour in sites and other in species)
  p <- p + ggnewscale::new_scale('size')
  #p <- p + ggnewscale::new_scale('shape')
  p <- p + ggnewscale::new_scale('alpha')
  p <- p + ggnewscale::new_scale('stroke')
  
  
  # Save plot into pass
  pass$plot <- p
  
  
  # Return pass
  return(pass)
}

