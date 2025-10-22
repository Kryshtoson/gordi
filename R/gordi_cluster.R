#'Connect observations and draw spiders/hulls
#'@description 
#'The function adds connections between observations using ggplot2::geom_path() (via `group`), and/or draws spider/hull segments from cluster centroid to each observation (via `cluster` and `spider`).
#'
#' @details 
#' To connect observations, you need to set argument `group` to a column in `env` dataframe (supplied in [gordi_read()]) that identifies observations to be connected. Lines are drawn with [ggplot2::geom_path()]. For drawing spiders, set argument `cluster` to a column in `env` dataframe (supplied in [gordi_read()]), and `spider = TRUE`. Centroids are computed as the mean of site scores within each cluster, and segments from centroids to each observation are drawn with [ggplot2::geom_segment()]. The aesthetics of connecting lines can be modified by multiple arguments.
#' 
#' @param pass A list object produced by [gordi_read()].
#' @param group Character; a grouping variable = a column in `env` dataframe (supplied in [gordi_read()]), used to connect observations.
#' @param cluster Character; column in `env` dataframe (supplied in [gordi_read()]), defining clusters for spiders.
#' @param spider Logical; draw spider segment for each centroid.
#' @param hull Logical; draw hulls. Does not work yet.
#' @param label Logical; when `TRUE` and `spider = TRUE`, draw one label per cluster at the centroid location.
#' @param linetype Modify the appearance of lines (paths, segments). Can be specified by numeric values (0-6) or by a name ('blank', 'solid', 'dashed', 'dotted', 'dotdash', 'longdash', 'twodash'). Default `'solid'`.
#' @param linewidth Numeric; set line width.
#' @param colour Changes site colours. Can be either: a variable = a specific column name of 'env' dataset, a constant value: R colour name, numeric value or hex (e.g. `'green'`, `1`, `'#6aa84f'`). When left as `''`, the internal default colours are used.
#' @param arrow Logical; draws arrows instead of lines.
#' @param arrow_type Character; type of arrow, either 'open' (default) or 'closed'.
#' @param arrow_length Numeric; set length of the arrowhead in centimeters. Deafault `0.3`.
#' @param arrow_ends Character; which ends to draw arrowheads on. One of 'last' (default), 'first', 'both'.
#' @param arrow_angle Numeric; angle of the arrowhead in degrees. Default `30`.
#' 
#' @return The updated `pass` object with added layers (either lines or spiders) to `pass$plot`.
#' 
#' @examples
#' # display site scores colour sites based on treatment, and connect the same sites
#' gordi_read(pco.bc, env = auch.env)|> gordi_sites(colour = 'Treatment')|> gordi_colour(scale = 'discrete', family = 'manual', values = c('darkorange', 'darkgreen'))|> gordi_cluster(group = 'site', linetype = 'dotted')
#'
#' #display NMDS, draw spiders to hydrology clusters and label centroids
#' gordi_read(nmds.3, spe = chiro, env = chiro.env)|> gordi_cluster(cluster = 'hydr', spider = T, colour = 'hydr', label = T)
#' 
#' @seealso [gordi_read()], [gordi_sites()], [gordi_colour()], [ggplot2::geom_path()], [ggplot2::geom_segment()], [ggplot2::geom_label()]
#' 
#' @importFrom ggplot2 ggplot aes geom_path geom_segment geom_label theme_bw labs theme element_text
#' @importFrom dplyr group_by summarise left_join
#' @importFrom rlang sym .data
#' @importFrom ggnewscale new_scale_colour
#' @importFrom grid arrow unit
#' @export

gordi_cluster <- function(pass,
                          group = '', #column in env table to connect observations (resurvey etc.)
                          cluster = '', #column in env table for clustering
                          spider = FALSE, #T/F if u want spiders or hulls
                          hull = NULL, #T/F if u want spiders or hulls
                          label = FALSE,
                          linetype = 'solid',
                          linewidth = 0.6,
                          colour = '',
                          arrow = FALSE,
                          arrow_type = 'open',
                          arrow_length = 0.3,
                          arrow_ends = 'last',
                          arrow_angle = 30,
                          show.legend = TRUE
){
  
  site_df <- bind_cols(pass$site_scores, pass$env)
  
  x_col <- names(pass$site_scores)[1]
  y_col <- names(pass$site_scores)[2]
  
  # actual labs
  if(pass$type %in% c('DCA', 'NMDS')) {actual_labs <- paste0(pass$axis_names)} else 
  {actual_labs <- paste0(pass$axis_names, " (", round(pass$explained_variation[pass$choices]*100, 2), "%)")}
  
  # plot set up  
  if (is.null(pass$plot)) { 
    p <- ggplot2::ggplot() +
      theme_bw() +
      geom_vline(aes(xintercept = 0), linetype = 3, linewidth = 0.2, colour = 'gray15', alpha = 0.6) +
      geom_hline(aes(yintercept = 0), linetype = 3, linewidth = 0.2, colour = 'gray15', alpha = 0.6) +
      labs(x = actual_labs[1], y = actual_labs[2]) +
      ggplot2::theme(
        text = element_text(size = 15),
        panel.grid = element_blank(),
        legend.justification = c(1, 1)
      )
  } else {
    p <- pass$plot
  }
  
  p <- p + ggnewscale::new_scale_colour()
  
  if (is.null(pass$env)){
    stop("This function works with `env` data. If you want to use this function, please provide `env` data in gordi_read().")
  }
  
  if (!group %in% names(site_df) && !identical(group, '')){
    stop("Column `", group, "` not found in the `env` data. Please provide valid input.")
  }
  
  if (!is.null(hull)){
    stop("Sorry, this argument is under construction...")
  }
  
  map_colour <- !identical(colour, '') && has_name(site_df, colour)
  # const_colour <- !identical(colour, '') && !map_colour && (grepl("^#(?:[A-Fa-f0-9]{6}[A-Fa-f0-9]{3})$", colour) || colour %in% grDevices::colours()) || (is.character(colour) && colour %in% palette()) || (is.numeric(colour) && colour %in% seq_along(palette()))
  
  const_colour <- !identical(colour, '') && !map_colour && (
    grepl("^#(?:[A-Fa-f0-9]{6}|[A-Fa-f0-9]{3})$", colour) ||
      (is.character(colour) && colour %in% grDevices::colours()) ||
      (is.numeric(colour)   && colour %in% seq_along(palette()))
  )
  
  if(!identical(colour, '')){
    message("To customize colours (similarly to `ggplot2::scale_colour_()` functions), please use gordi_colour() right after `gordi_sites()`.")
    if(!map_colour && !const_colour){
      warning("`colour` must be either a column in the `env` dataframe, a valid R colour name/hex code, or a numeric code! Ignoring input, default is being used.")
      colour <- ''
    }
  }
  
  arrow_spec <- if (isTRUE(arrow)){
    grid::arrow(length = grid::unit(arrow_length, 'cm'),
                type = arrow_type,
                ends = arrow_ends,
                angle = arrow_angle)
  } else NULL
  
  if (!identical(group, '')){
    aes_args <- list(
      x = rlang::sym(x_col),
      y = rlang::sym(y_col),
      group = rlang::sym(group)
    )
    if (map_colour) aes_args$colour <- sym(colour)
    
    const_args<- list()
    if(!map_colour) {
      if(!identical(colour, '')) {
        const_args$colour <- colour} else {const_args$colour <- 1}}
    
    # mapping <-  if (map_colour) {
    #   ggplot2::aes(
    #     x = .data[[x_col]],
    #     y = .data[[y_col]],
    #     group = .data[[group]],
    #     colour = !!rlang::sym(colour))}
    # else {
    #   ggplot2::aes(
    #     x = .data[[x_col]],
    #     y = .data[[y_col]],
    #     group = .data[[group]]
    #   )
    # }
    # p <- p + do.call(geom_segment, c(list(mapping = do.call(aes, aes_args), data = site_centroids), const_args, linetype = linetype, linewidth = linewidth))
    p <- p + do.call(geom_path, c(list(mapping = do.call(aes, aes_args), data = site_df, linetype = linetype, linewidth = linewidth, arrow = arrow_spec, inherit.aes = FALSE, show.legend = show.legend), const_args))
    #  p <- p + geom_path(data = site_df, mapping = mapping, colour = if (const_colour) colour else NULL, linetype = linetype, linewidth = linewidth, arrow = arrow_spec, inherit.aes = FALSE)
  }
  
  if (isTRUE(spider) && !identical(cluster, '')){
    centroids <- site_df|>
      group_by(.data[[cluster]])|>
      summarise(.xc = mean(.data[[x_col]]), .yc = mean(.data[[y_col]]), .groups = 'drop')
    
    
    site_centroids <- site_df|>
      #select(all_of(s.cols))|>
      left_join(centroids, by = cluster)
    
    
    aes_args <- list(
      x = quote(.xc),
      y = quote(.yc),
      xend = rlang::sym(x_col),
      yend = rlang::sym(y_col),
      group = rlang::sym(cluster)
    )
    if (map_colour) aes_args$colour <- sym(colour)
    
    const_args<- list()
    if(!map_colour) {
      if(!identical(colour, '')) {
        const_args$colour <- colour} else {const_args$colour <- 1}}
    
    p <- p + do.call(geom_segment, c(list(mapping = do.call(aes, aes_args), data = site_centroids), const_args, linetype = linetype, linewidth = linewidth, show.legend = show.legend))
  }
  if (isTRUE(label) && isTRUE(spider) && !identical(cluster, '')){
    p <- p + geom_label(data = centroids, mapping = aes(x = .data[['.xc']], y = .data[['.yc']], label = .data[[cluster]]), inherit.aes = F, show.legend = show.legend)
  }
  
  pass$plot <- p 
  
  return(pass)
}
