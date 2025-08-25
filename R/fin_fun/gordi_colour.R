gordi_colour <- function(pass,
                         scale = c('continuous', 'discrete','auto', 'binned'), #what colour scale to use
                         family = c('viridis', 'brewer', 'gradient', 'steps', 'manual', 'default'),
                         fill = FALSE,
                         breaks = waiver(), name = waiver(), labels = waiver(), #original scale_ arguments
                         limits = NULL, na.value = NA, guide = waiver(), trans = 'identity', 
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

