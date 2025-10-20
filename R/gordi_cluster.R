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
                          arrow_angle = 30
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
    p <- p + do.call(geom_path, c(list(mapping = do.call(aes, aes_args), data = site_df, linetype = linetype, linewidth = linewidth, arrow = arrow_spec, inherit.aes = FALSE), const_args))
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
    
    p <- p + do.call(geom_segment, c(list(mapping = do.call(aes, aes_args), data = site_centroids), const_args, linetype = linetype, linewidth = linewidth))
  }
  if (isTRUE(label) && isTRUE(spider) && !identical(cluster, '')){
    p <- p + geom_label(data = centroids, mapping = aes(x = .data[['.xc']], y = .data[['.yc']], label = .data[[cluster]]), inherit.aes = F)
  }
  
  pass$plot <- p 
  
  return(pass)
}
