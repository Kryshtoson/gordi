
# gordi_sites master version ----------------------------------------------
gordi_sites <- function(pass, label = '', fill = '', alpha = '', stroke = '', shape = '', size = '', colour = '', repel_label = T) {
  
  #' axis names
  names(pass$site_scores) <- paste0("Axis_site", 1:2)
  names(pass$species_scores) <- paste0("Axis_spe", 1:2)
  
  #' selection of ordination type
  
  if(pass$type == 'CCA'){
    actual_labs <- paste0("CCA", pass$choices, " (", round(pass$explained_variation[1:2]*100, 2), '%)')
  } else if(pass$type == 'CA'){
    actual_labs <- paste0("CA", pass$choices, " (", round(pass$explained_variation[1:2]*100, 2), '%)')
  } else if (pass$type == 'PCA'){
    actual_labs <- paste0("PCA", pass$choices, " (", round(pass$explained_variation[1:2]*100, 2), '%)')
  } else if(pass$type == 'PCoA'){
    actual_labs <- paste0("PCoA", pass$choices, " (", round(pass$explained_variation[1:2]*100, 2), '%)')
  } else if(pass$type == 'db-RDA'){
    actual_labs <- paste0("db-RDA", pass$choices, " (", round(pass$explained_variation[1:2]*100, 2), '%)')
  } else if(pass$type == 'RDA'){
    actual_labs <- paste0("RDA", pass$choices, " (", round(pass$explained_variation[1:2]*100, 2), '%)')
  } else if (pass$type == 'DCA'){
    actual_labs <- paste0('DCA', pass$choices)
  } else if (pass$type == 'NMDS'){
    actual_labs <- paste0('NMDS', pass$choices)
  }
  
  #' plot set up  
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
  
  site_df <- bind_cols(pass$env, pass$site_scores)
  
  #' accounting for labeling
  if (label != '') {
    if (repel_label) {
      p <- p + geom_text_repel(data = site_df, aes(Axis_site1, Axis_site2, label = !!sym(label)))
    } else {
      p <- p + geom_text(data = site_df, aes(Axis_site1, Axis_site2, label = !!sym(label)))
    }
  }
  
  #' accounting for indiviidual geom scales  
  p <- p + ggnewscale::new_scale("size") 
  p <- p + ggnewscale::new_scale("shape") 
  p <- p + ggnewscale::new_scale("fill")
  p <- p + ggnewscale::new_scale("alpha")
  p <- p + ggnewscale::new_scale("stroke") 
  p <- p + ggnewscale::new_scale_colour()
  
  
  # Detect mapped vs constant aesthetics
  map_colour <- !identical(colour, '') && has_name(site_df, colour)
  const_colour <- const_colour <- !map_colour && (grepl("^#(?:[A-Fa-f0-9]{6}[A-Fa-f0-9]{3})$", colour) || colour %in% grDevices::colours()) || (is.character(colour) && colour %in% palette()) || (is.numeric(colour) && colour %in% seq_along(palette()))
  
  map_size <- !identical(size, '') && has_name(site_df, size)
  const_size <- !map_size && is.numeric(size)
  
  map_shape <- !identical(shape, '') && has_name(site_df, shape)
  const_shape <- !map_shape && is.numeric(shape)
  
  map_fill <- !identical(fill, '') && has_name(site_df, fill)
  const_fill <- !map_fill && (grepl("^#(?:[A-Fa-f0-9]{6}[A-Fa-f0-9]{3})$", fill) || fill %in% grDevices::colours())
  
  map_alpha <- !identical(alpha, '') && has_name(site_df, alpha)
  const_alpha <- !map_alpha && is.numeric(alpha)
  
  map_stroke <- !identical(stroke, '') && has_name(site_df, stroke)
  const_stroke <- !map_stroke && is.numeric(stroke)
  
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
  
  #' plot  
  p <- p + do.call(geom_point, c(list(mapping = do.call(aes, aes_args_point), data = site_df), const_args_point))
  
  pass$plot <- p
  
  return(pass)
}
