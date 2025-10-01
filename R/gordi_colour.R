#'Set colour/fill scales for gordi plots
#'
#'@description
#'Adds a colour or fill scale to the current ggplot plot stored in `pass$plot` (created by gordi_*() functions). 
#'After applying the colour scale, it starts a new independent scale using `ggnewscale`, 
#'so later layers can have their own scale without any interference.
#'
#'@details
#' use fill = FALSE to modify colour scales (default) and fill = TRUE to modify the fill scales
#' Compulsory arguments:
#' - `scale` selects the scale type: `continuous`, `discrete`, `binned` or `auto`. 
#' With `scale = 'auto'`, the function chooses scale `discrete` when `values` are specified (e.g. manual palette), 
#' otherwise `continuous`.
#' - `family` selects the scale family within the scale type:
#' for scale `discrete` -> `'manual'`, `'viridis'`, `'brewer'`, `'default'`
#' for scale `continuous` ->  `'viridis'`, `'gradient'`, `'steps'`, `'brewer'`, `'default'`
#' for scale `binned` -> `'viridis'`, `'brewer'`, `'steps'`, `'default'`
#' For `scale = 'discrete'`, `family = 'manual'` you must provide `values` (vector of colours)
#' Argument `na.values` sets colours for NA values in the dataset.
#' After adding a certain scale, the function set a new colour/fill scale with `ggnewscale`, so subsequent layers can use a fresh scale.
#' 
#' Call [gordi_colour()] right after you have added a layer that maps colour/fill (e.g. `gordi_species(colour = ...)`, `gordi_sites(fill = ...)`).
#' If no layer uses colour/fill aesthetic, this function has no effect.
#' 
#' @param pass A list object produced by [gordi_read()] 
#' @param scale Character; One of the `c('continuous', 'discrete','auto', 'binned')`.
#'    Character.
#' @param family Character; One of the palette families `c('viridis', 'brewer', 'gradient', 'steps', 'manual', 'default')`.
#'    Character.
#' @param fill Logical; If fill = FALSE operates on colour scale, if fill = TRUE operates on fill scale.
#'    Default is FALSE.
#' @param breaks Scale breaks. See `scale_colour_...` or `scale_fill_...`.
#' @param name Scale name. See `scale_colour_...` or `scale_fill_...`.
#' @param labels Labels for scale. See `scale_colour_...` or `scale_fill_...`.
#' @param limits Scale limits. See `scale_colour_...` or `scale_fill_...`.
#' @param guide Guide for the scale. See `scale_colour_...` or `scale_fill_...`.
#' @param trans Transformation for the scale. See `scale_colour_...` or `scale_fill_...`.
#' @param na.value Colour used for NA values.
#' @param values Vector of colours for `manual` discrete scale, and for continuous gradientn and stepsn variants.
#' @param option, direction, begin, end, alpha  Options for `viridis` scales.
#' @param low, mid, midpoint, high Colours for `gradient`/`gradient2` or `steps`/`steps2` scales.
#' @param palette_name Arguments for `brewer` scales.
#' @param bins, n.breaks Arguments for `binned`/`steps` scales.
#' 
#' @return The updated `pass` object with label layers added to `pass$plot`. Followed by a fresh `ggnewscale` for the same aesthetic (colour or fill).
#' 
#' @examples
#' # discrete species colours with a brewer palette
#' gordi_read(m, env, traits) |>
#'   gordi_species(colour = 'form') |>
#'   gordi_colour(scale = 'discrete', family = 'brewer', palette_name = 'Set2')
#' 
#' # continuous viridis fill for sites
#' gordi_read(m, env) |>
#'   gordi_sites(fill = 'elevation') |>
#'   gordi_colour(fill = TRUE, scale = 'continuous', family = 'viridis')
#' @export
gordi_colour <- function(pass,
                         scale = c('continuous', 'discrete','auto', 'binned'), #what colour scale to use
                         family = c('viridis', 'brewer', 'gradient', 'steps', 'manual', 'default'),
                         fill = FALSE,
                         breaks = waiver(), name = waiver(), labels = waiver(), #original scale_ arguments
                         limits = NULL, na.value = 'grey50', guide = waiver(), trans = 'identity', 
                         values = NULL, #used by manual (discrete) and continuous (gradientn, stepsn)
                         option = NULL, direction = 1, begin = 0, end = 1, alpha = 1, #viridis
                         low = NULL, mid = NULL, high = NULL, midpoint = NULL, #gradient, gradient2, steps, steps2
                         palette_name = NULL, type = NULL, # brewer palettes
                         bins = NULL, n.breaks = NULL #binned/steps colours
){
  
  scale <- match.arg(scale)
  family <- match.arg(family)
  
  sc <- function(suffix){
    paste0('scale_', if (isTRUE(fill)) 'fill' else 'colour', '_', suffix)
  }
  
  is_hex   <- function(x) 
    is.character(x) && grepl("^#(?:[A-Fa-f0-9]{3}|[A-Fa-f0-9]{6})$", x) #checks whether x is character as well as hex colour code
  
  is_colname <- function(x) 
    is.character(x) && x %in% grDevices::colours() #checks whether x is a character as well as colour name in r ('red', 'blue'...)
  
  is_numcol <- function(x)
    is.numeric(x)  && x <= length(palette()) #checks if x is a colour number
  
  is_col   <- function(x) 
    is_hex(x) || is_colname(x) || is_numcol(x) #checks if x is hex colour code or colour name in R
    
    all_cols <- function(v) {
      if(is.numeric(v)){
        all(vapply(v, is_numcol,, logical(1)))
        } else if(is.character(v)){
            length(v) > 0 && all(vapply(v, is_col, logical(1)))} else {FALSE}
          } #when values = c('red', '#FF00', 1,..)
    
    if(scale == 'discrete' && family == 'manual'){
      message("To customise discrete manual colours, supply `values` (e.g. values = c('red', 'green',...).")
      if(is.null(values) || !all_cols(values)){
        stop("For `scale = 'discrete' and family = 'manual'`, `values` must be a valid vector of colours (names or hex code) or a numeric code. Ignoring input.")
      }
    }
    
    viridis_options <- c('magma', 'inferno', 'plasma', 'viridis', 'cividis', 'rocket', 'mako', 'turbo', 'A', 'B', 'C', 'D', 'E', 'F', 'G', 'H')
    if (family == 'viridis' && !is.null(option)){
      message("To customize viridis colour scales, please use argument `option`, valid values are magma, inferno, plasma, viridis, cividis, rocket, mako, turbo and letters A-H, e.g. `option = 'viridis'`")
    if (!option %in% viridis_options){
      warning("`option = ", option, "` is not a recognized viridis colour palette. Default palette is being used.")
      option <- NULL
    }
    }
    
    if (family ==  'brewer'){
      message("To customise brewer colour palettes, please set `palette_name` (e.g. palette_name = 'Set1'). See RColorBrewer::display.brewer.all() for more colour options. ")
      if (!is.null(palette_name) && !palette_name %in% rownames(RColorBrewer::brewer.pal.info)){
        warning("`palette_name = ", palette_name, "` not found in RColorBrewer. Using default palette colour ('Set1' for discrete, 'YlGnBu' for continuous).")
        palette_name <- NULL
      }
    }
    
    #direction of colour palette
    if (!is.null(direction) && !direction %in% c(-1, 1)){
      message("`direction` controls palette order. Use 1 (default) or -1 (reversed).")
      warning("`direction = ", direction, "` is invalid. Default `direction = 1` is being used.")
      direction <- 1
    }
    
    #corrected  hue of the colour palette
    if (!is.null(begin) && (begin < 0 || begin > 1)){
      warning("`begin` must be between 0 and 1. `begin = ", begin, "` is out of range. Setting `begin` to default 0.")
      begin <- 0
    }
    
    if (!is.null(end) && (end < 0 || end > 1)){
      warning("`end` must be between 0 and 1. `end = ", end, "` is out of range. Setting `end` to default 0.")
      end <- 0
    }
    
    #checking for alpha
    if (!is.null(alpha) && (alpha < 0 || alpha > 1)){
      warning("`alpha` controls transparency and must be between 0 and 1. `alpha = ", alpha, "` is out of range. Ignoring input, default `alpha = 1` is being used.")
      alpha <- 1
    }
    
    if (!is.null(n.breaks) && (!is.numeric(n.breaks) || n.breaks <= 1)){
      warning("`n.breaks` controls the number of breaks used in steps colour scales. `n.breaks = `", n.breaks, "` is invalid. Ignoring input, ggplot will choose the breaks.")
      n.breaks <- NULL
    }
    
    #type in scale_colour_fermenter and scale_colour_brewer
    valid_types <- c('seq', 'div', 'qual')
    
    if (!is.null(type) && !type %in% valid_types){
      warning("For brewer palettes, `type` should match the palette family: 'seq' (sequential), 'div' (diverging), 'qual' (qualitative). `type`must be one of 'seq', 'div', 'qual'. Ignoring input, letting ggplot infer the type based on palette.")
      type <- NULL
    }
    
    if (isFALSE(fill)){
      message("`fill = FALSE` -> `gordi_colour` works with `scale_colour_...()` functions.")
    } else {
      message("`fill = TRUE` -> `gordi_colour` works with `scale_fill_...()` functions.")
    }

  
  #function that builds ggplot scale_colour object
  build_scale_obj <- function() {
    #if scale = auto, check values and decide whether discrete or continuous based on them
    if (scale == 'auto'){
      scale <- if(!is.null(values) && is.character(values)) 'discrete' else 'continuous'
    }
    
    if (scale == 'discrete'){
      if (family == 'manual'){
        #values must be provided for family = manual
        #stopifnot(!is.null(values))
        if (is.null(values)){
          stop("gordi_colour with scale = 'discrete' and family = 'manual': you must provide `values`, e.g. values = c('red', 'blue').")
        }
        args <- purrr::compact(list(values = values, name = name, breaks = breaks, labels = labels,
                                    limits = limits, guide = guide, na.translate = TRUE, drop = TRUE,
                                    na.value = na.value))
        return(do.call(sc('manual'), args))
      }
      if (family == 'viridis'){
        args <- purrr::compact(list(option = option, direction = direction, begin = begin, end = end, alpha = alpha,
                                    name = name, breaks = breaks, labels = labels, limits = limits, guide = guide,
                                    na.value = na.value))
        return(do.call(sc('viridis_d'), args))
      }
      if (family == 'brewer'){
        #if family brewer and no palette_name provided set it to Set1
        args <- purrr::compact(list(palette = if (is.null(palette_name)) "Set1" else palette_name, type = type, direction = direction,
                                    name = name, breaks = breaks, labels = labels, limits = limits, guide = guide,
                                    na.value = na.value))
        return(do.call(sc('brewer'), args))
      }
      if (family == 'default'){
        #if family = default just use default ggplot for discrete colours scale_colour_discrete
        args <- purrr::compact(list(name = name, breaks = breaks, na.value = na.value, labels = labels, limits = limits, guide = guide))
        return(do.call(sc('discrete'), args))
      }
      stop('For discrete scales please use family = c(manual, viridis, brewer, default)')
    }
    
    if (scale == 'continuous'){
      #normal viridis
      if(family == 'viridis'){
        args <- purrr::compact(list(option = option, direction = direction, begin = begin, end = end, alpha = alpha,
                                    name = name, breaks = breaks, labels = labels, limits = limits, guide = guide,
                                    trans = trans, na.value = na.value))
        return(do.call(sc('viridis_c'), args))
      }
      if (family == 'gradient'){
        #if values defined uses scale_colour_gradientn
        if (!is.null(values)){
          args <- purrr::compact(list(colours = values, name = name, breaks = breaks, labels = labels,
                                      limits = limits, guide = guide, trans = trans, na.value = na.value))
          return(do.call(sc('gradientn'), args))
        }
        if (!is.null(mid) || !is.null(midpoint)){
          #if mid/midpoint define use gradient2 otherwise gradient
          args <- purrr::compact(list(low = if (is.null(low)) "#2166AC" else low, mid = if (is.null(mid)) 'grey90' else mid, high = if (is.null(high)) "#B2182B" else high,
                                      midpoint = if (is.null(midpoint)) 0 else midpoint,
                                      name = name, breaks = breaks, labels = labels, limits = limits,
                                      guide = guide, trans = trans, na.value = na.value))
          return(do.call(sc(gradient2), args))
        }
        args <- purrr::compact(list(low = if (is.null(low)) "#132B43" else low, high = if (is.null(high)) "#56B1F7" else high,
                                    name = name, breaks = breaks, labels = labels, limits = limits,
                                    guide = guide, trans = trans, na.value = na.value))
        return(do.call(sc('gradient'), args))
      }
    }
    
    if (family == 'brewer'){
      args <- purrr::compact(list(palette = if(is.null(palette_name)) "YlGnBu" else palette_name, type = type, direction = direction,
                                  name = name, breaks = breaks, labels = labels, limits = limits,
                                  guide = guide, trans = trans, na.value = na.value))
      return(do.call(sc('distiller'), args))
    }
    
    if (family == 'steps'){
      #steps with custom colour values
      if (!is.null(values)) {
        args <- purrr::compact(list(colours = values, name = name, breaks = breaks, labels = labels, limits = limits,
                                    guide = guide, trans = trans, na.value = na.value, n.breaks = n.breaks))
        return(do.call(sc('stepsn'), args))
      }
      if (!is.null(mid) || !is.null(midpoint)) {
        #multiple colour steps (as in gradient2)
        args <- purrr::compact(list(low = if(is.null(low)) "#2166AC" else low, mid = if (is.null(mid)) "grey90" else mid, high = if (is.null(high)) "#B2182B" else high,
                                    midpoint = if (is.null(midpoint)) 0 else midpoint,
                                    name = name, breaks = breaks, labels = labels, limits = limits,
                                    guide = guide, trans = trans, na.value = na.value, n.breaks = n.breaks))
        return(do.call(sc('steps2'), args))
      }
      #steps as in gradient
      args <- purrr::compact(list(low = low, high = high,
                                  name = name, breaks = breaks, labels = labels, limits = limits,
                                  guide = guide, trans = trans, na.value = na.value, n.breaks = n.breaks))
      return(do.call(sc('steps'), args))
    }
    
    if (family == "default") {
      args <- purrr::compact(list(name = name, breaks = breaks, labels = labels, limits = limits,
                                  guide = guide, trans = trans, na.value = na.value))
      return(do.call(sc('continuous'), args))
      stop("For continuous scales, use family = c('viridis', 'brewer', 'gradient', 'steps', 'default'.")
    }
    
    #binned colour scales
    if (scale == "binned") {
      if (family == "viridis") {
        args <- purrr::compact(list(option = option, direction = direction, begin = begin, end = end, alpha = alpha,
                                    name = name, breaks = breaks, labels = labels, limits = limits, guide = guide,
                                    trans = trans, na.value = na.value))
        return(do.call(sc('viridis_b'), args))
      }
      if (family == "brewer") {
        # brewer binned -> fermenter in ggplot2
        args <- purrr::compact(list(palette = if (is.null(palette_name)) "YlOrRd" else palette_name, type = type, direction = direction,
                                    bins = bins,
                                    name = name, breaks = breaks, labels = labels, limits = limits,
                                    guide = guide, trans = trans, na.value = na.value))
        return(do.call(sc('fermenter'), args))
      }
      if (family == "steps") {
        # steps family already returns binned scales (same mapping as above)
        if (!is.null(values)) {
          args <- purrr::compact(list(colours = values, name = name, breaks = breaks, labels = labels, limits = limits,
                                      guide = guide, trans = trans, na.value = na.value, n.breaks = n.breaks))
          return(do.call(sc('stepsn'), args))
        }
        if (!is.null(mid) || !is.null(midpoint)) {
          args <- purrr::compact(list(low = if (is.null(low)) "#2166AC" else low, mid = if (is.null(mid)) "grey90" else mid, high = if (is.null(high)) "#B2182B" else high,
                                      midpoint = if (is.null(midpoint)) 0 else midpoint,
                                      name = name, breaks = breaks, labels = labels, limits = limits,
                                      guide = guide, trans = trans, na.value = na.value, n.breaks = n.breaks))
          return(do.call(sc('steps2'), args))
        }
        args <- purrr::compact(list(low = low, high = high,
                                    name = name, breaks = breaks, labels = labels, limits = limits,
                                    guide = guide, trans = trans, na.value = na.value, n.breaks = n.breaks))
        return(do.call(sc('steps'), args))
      }
      if (family == "default") {
        args <- purrr::compact(list(name = name, breaks = breaks, labels = labels, limits = limits,
                                    guide = guide, trans = trans, na.value = na.value, n.breaks = n.breaks))
        return(do.call(sc('binned'), args))
      }
      stop("For binned scales, use family = c('viridis', 'brewer', 'steps', 'default').")
    }
    
    stop("Unknown combination of scale and family. Please check `scale` and `family`.")
    
  }
  
  scale_obj <- build_scale_obj()
  pass$plot <- pass$plot + scale_obj
  
  pass$plot <- pass$plot + if (isTRUE(fill)) ggnewscale::new_scale_fill()
                               else ggnewscale::new_scale_colour()
  
  return(pass)
}

