#' @description
#' The function [gordi_plot()] returns the ggplot object created by functions:[gordi_sites()] and/or [gordi_species()] after [gordi_read()].
#' 
#' @details
#' You can either save the plot via `filename` (wrapper over [ggplot2::ggsave()]), or just get the ggplot object to customize it further with ggplot syntax (e.g. `+ theme(...)`, etc.).
#' 
#' @param pass A list object produced by [gordi_read()]
#' @param filename Optional path to save your plot (e.g. "plot.png"), similar to [ggplot2::ggsave()]. If NULL, no file is created.
#' @param customize Logical; if TRUE prints a short message about customization.
#' @param width Integer; passed to [ggplot2::ggsave()] when `filename` is provided. Used to customize plot size expressed by `units`  argument when saving, as in [ggplot2::ggsave()].
#' @param height Integer; passed to [ggplot2::ggsave()] when `filename` is provided. Used to customize plot size expressed by `units` argument when saving, as in [ggplot2::ggsave()].
#' @param units String; allowed units: 'in', 'cm', 'mm', 'px', in which the `width` and `height` arguments are expressed. Passed to [ggplot2::ggsave()] when `filename` is provided. 
#' @param dpi Integer; passed to [ggplot2::ggsave()] when `filename` is provided. Used to customize plot resolution when saving, as in [ggplot2::ggsave()].
#' @param device Device to use. Similar to `device` argument in [ggplot2::ggsave()]. If NULL (default), the device is guessed based on the `filename` extension.
#' 
#' @return A ggplot object (even when saving to filename).
#' 
#' @examples
#' # get a sites ggplot object and further customize
#' p <- gordi_read(m, dune.env)|> gordi_sites()|> gordi_plot()
#' p + theme_minimal()
#' # customization of sites ggplot object in one line
#' gordi_read(m, dune.env)|> gordi_sites()|> gordi_plot()+ theme_minimal()
#' # saving species ggplot object as png file in the directory of the current R-project (still returns the ggplot object)
#' gordi_read(m, dune.env)|> gordi_species()|> gordi_plot(filename = 'plot.png')
#' 
#' @seealso [gordi_read()], [gordi_sites()], [gordi_species()], [ggplot2::ggsave()]
#' @export
gordi_plot <- function(pass,
                       filename = NULL,
                       customize = TRUE,
                       width = 8,
                       height = 6,
                       units = 'in',
                       dpi = 300, 
                       device = NULL){
  p <- pass$plot
  
  valid_devices <- c('png', 'pdf', 'jpeg', 'tiff', 'svg', 'bmp')
  if (!is.null(device) && !device %in% valid_devices){
    warning("Specified `device` is not valid (valid devices: png, pdf, jpeg, tiff, svg, bmp). Invalid input, using default device 'png'.")
    device <- 'png'
  }
  
  if(!is.null(filename)){
    ggplot2::ggsave(filename = filename, plot = p, width = width, height = height, units = units, dpi = dpi, device = device)
  if(isTRUE(customize)) message("Plot saved to `", filename, "` Returning ggplot object for further customization.")
  } else if (isTRUE(customize)){
    message("Returning ggplot object for customization. You can simply continue with `ggplot2` syntax -> `+ theme(...)`, etc.")
  }
  
  p
}
