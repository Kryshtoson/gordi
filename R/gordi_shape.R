#'@description
#'The function [gordi_shape()] wraps ggplot functions: [ggplot2::scale_shape_manual()], [ggplot2::scale_shape_identity()] and [ggplot2::scale_shape_discrete()] to enable change of point shapes in the ggplot object stored in `pass$plot` (created by `gordi_*()` functions). After applying the function, it creates a new independent shape scale with `ggnewscale`, so later layer can have their own scale without any interference.
#'
#'@details
#'Use `scale = 'discrete'` when your mapped shape aesthetic is a categorical variable. You can either let ggplot pick shapes automatically (`values = NULL` -> uses [scale_shape_discrete()]) or supply a manual vector of shape codes (0-25) -> uses [scale_shape_manual()]. Note that shape discrete palette can deal with a maximum of 6 discrete values. Use `scale = 'identity'` only when the mapped shape column already contains shapes you want to map. This takes the data as shapes and maps them.
#'
#'@param pass A list object produced by [gordi_read()].
#'@param scale Character; 'discrete' to treat mapped values as categorical or 'identity' to treat them as actual shapes.
#'@param values Optional vector of shape codes (0-25) for manual discrete scale. Note that shape discrete palette can deal with maximum of 6 discrete values.
#'@param breaks Passed to ggplot2 scale functions. Controls legend breaks.
#'@param name Passed to ggplot2 scale functions. Controls legend title.
#'@param labels Passed to ggplot2 scale functions. Controls legend labels.
#'@param limits Passed to ggplot2 scale functions. Order or subset of levels to show.
#'@param guide Passed to ggplot2 scale functions. Legend guide specification.
#'@param drop Logical; drop unused levels from the legend (default = `TRUE`).
#'@param na.translate Logical; include an `NA` key in the legend when present (default = `TRUE`).
#'
#'@return The updated `pass$plot` object with the added shape scale and a fresh shape scale started for subsequent layers.
#'
#'@examples
#' # manually change discrete shapes of species points
#' gordi_read(m, env, traits)|> gordi_species(shape = 'group')|> gordi_shape(scale = 'discrete', values = c(12, 16, 18))
#' # mapping site point shapes as actual data values
#' gordi_read(m, env, traits)|> gordi_sites(shape = 'logger_ID')|> gordi_shape(scale = 'identity')
#' @seealso [gordi_read()], [gordi_sites()], [gordi_species()]
#' @export
gordi_shape <- function(pass,
                        scale = c('discrete', 'identity'),
                        values = NULL,
                        breaks = waiver(),
                        name = waiver(),
                        labels = waiver(),
                        limits = NULL,
                        guide = waiver(),
                        drop = TRUE,
                        na.translate = TRUE) {
  
  scale <- match.arg(scale)
  sc <- function(suffix) paste0('scale_shape_', suffix)

  
  build_scale_obj <- function(){
    if (scale == 'identity'){
      args <- purrr::compact(list(name = name, breaks = breaks, labels = labels, limits = limits, guide = guide))
      return(do.call(sc('identity'), args))
    }
    
    
    if(!is.null(values)){
      args <- purrr::compact(list(values = values, name = name, breaks = breaks, labels = labels, limits = limits, guide = guide, drop = drop, na.translate = na.translate))
      return(do.call(sc('manual'), args))
    }
    else {
      args <- purrr::compact(list(name = name, breaks = breaks, labels = labels, limits = limits, guide = guide, drop = drop, na.translate = na.translate))
      return(do.call(sc('discrete'), args))
    }
  }
  
  scale_obj <- build_scale_obj()
  pass$plot <- pass$plot + scale_obj
  
  pass$plot <- pass$plot + ggnewscale::new_scale('shape')
  
  return(pass)
  
}
