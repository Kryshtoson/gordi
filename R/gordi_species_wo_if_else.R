gordi_species <- function(pass,
                          label = '',
                          colouring = '',
                          size = '',
                          shape = '',
                          fill = '',
                          alpha = '',
                          stroke = '',
                          repel_label = T) {
  
  
  names(pass$species_scores) <- paste0("Axis_spe", 1:2)
  names(pass$site_scores) <- paste0("Axis_site", 1:2)
  
  if(pass$type == 'CCA'){
    actual_labs <- paste0("CCA", pass$choices, " (", round(pass$explained_variation[1:2]*100, 2), '%)')
  } else if(pass$type == 'CA'){
    actual_labs <- paste0("CA", pass$choices, " (", round(pass$explained_variation[1:2]*100, 2), '%)')
  } else if (pass$type == 'PCA'){
    actual_labs <- paste0("PCA", pass$choices, " (", round(pass$explained_variation[1:2]*100, 2), '%)')
  } else if(pass$type %in% c('PCoA (capscale)', 'PCoA (rda)')){
    actual_labs <- paste0("PCoA", pass$choices, " (", round(pass$explained_variation[1:2]*100, 2), '%)')
  } else if(pass$type %in% c('db-RDA (rda)', 'db-RDA (capscale)')){
    actual_labs <- paste0("db_RDA", pass$choices, " (", round(pass$explained_variation[1:2]*100, 2), '%)')
  } else if(pass$type == 'RDA constrained'){
    actual_labs <- paste0("RDA", pass$choices, " (", round(pass$explained_variation[1:2]*100, 2), '%)')
  }
  
  if (is.null(pass$plot)) { # checks whether p exists in pass, if not it draws plot
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
  
  spe_df <- bind_cols(pass$species_names, pass$species_scores)
  
  spe_df <- spe_df |> 
    left_join(pass$traits, by = join_by(!!sym(names(spe_df)[1]) == !!sym(names(pass$traits)[1])))
  
  # Define defaults
  default_colour <- "black"
  default_size <- 3
  default_shape <- 16
  default_fill <- "white"
  default_alpha <- 1
  default_stroke <- 0.5
  
  
  # Detect mapped vs constant aesthetics
  map_colour <- !identical(colouring, '') && has_name(spe_df, colouring)
  const_colour <- !map_colour && (grepl("^#(?:[A-Fa-f0-9]{6}[A-Fa-f0-9]{3})$", colouring) || colouring %in% grDevices::colours())
  
  map_size <- !identical(size, '') && has_name(spe_df, size)
  const_size <- !map_size && is.numeric(size)
  
  map_shape <- !identical(shape, '') && has_name(spe_df, shape)
  const_shape <- !map_shape && is.numeric(shape)
  
  map_fill <- !identical(fill, '') && has_name(spe_df, fill)
  const_fill <- !map_fill && (grepl("^#(?:[A-Fa-f0-9]{6}[A-Fa-f0-9]{3})$", fill) || fill %in% grDevices::colours())
  
  map_alpha <- !identical(alpha, '') && has_name(spe_df, alpha)
  const_alpha <- !map_alpha && is.numeric(alpha)
  
  map_stroke <- !identical(stroke, '') && has_name(spe_df, stroke)
  const_stroke <- !identical(stroke, '') && !map_stroke && is.numeric(stroke)
  
  
  # Prepare aes arguments (only mapped)
  aes_args_point <- list(x = quote(Axis_spe1), y = quote(Axis_spe2))
  if(map_colour) aes_args_point$colour <- sym(colouring)
  if(map_size) aes_args_point$size <- sym(size)
  if(map_shape) aes_args_point$shape <- sym(shape)
  if(map_fill) aes_args_point$fill <- sym(fill)
  if(map_alpha) aes_args_point$alpha <- sym(alpha)
  if(map_stroke) aes_args_point$stroke <- sym(stroke)
  
  # Prepare constant arguments (mapped first, then defaults if nothing)
  const_args_point <- list()
  
  # colouring 
  if(!map_colour){
    if(!identical(colouring, '')){
      const_args_point$colour <- colouring} else {const_args_point$colour <- default_colour}}
  # size 
  if(!map_size){
    if(!identical(size, '')){
      const_args_point$size <- size} else {const_args_point$size <- default_size}}
  # shape 
  if(!map_shape){
    if(!identical(shape, '')){
      const_args_point$shape <- shape} else {const_args_point$shape <- default_shape}}
  # fill 
  if(!map_fill){
    if(!identical(fill, '')){
      const_args_point$fill <- fill} else {const_args_point$fill <- default_fill}}
  # alpha 
  if(!map_alpha){
    if(!identical(alpha, '')){
      const_args_point$alpha <- alpha} else {const_args_point$alpha <- default_alpha}}
  # stroke
  if(!map_stroke){
    if(!identical(stroke, '')){
      const_args_point$stroke <- stroke} else {const_args_point$stroke <- default_stroke}}

  
  # Add the layer
  p <- p + do.call(geom_point, c(list(mapping = do.call(aes, aes_args_point), data = spe_df), const_args_point))
  
  



  
    #' accounting for labeling
  if (label != '') {
    if (repel_label) {
      p <- p + geom_text_repel(data = spe_df, aes(Axis_spe1, Axis_spe2, label = !!sym(label)))
    } else {
      p <- p + geom_text(data = spe_df, aes(aes(Axis_spe1, Axis_spe2), label = !!sym(label)))
    }
  }
  
  # possibility for two colour scales
  p <- p + ggnewscale::new_scale_colour()
  p <- p + ggnewscale::new_scale('size')
  p <- p + ggnewscale::new_scale('shape')
  p <- p + ggnewscale::new_scale_fill()
  p <- p + ggnewscale::new_scale('alpha')
  p <- p + ggnewscale::new_scale('stroke')
  
  
  
  pass$plot <- p
  
  return(pass)
}

geom_line()

gordi_read(m, env, trait) |> 
 # gordi_sites(colouring = 'red') |> 
  gordi_species(colour = 'form', size = 'cont', shape = 'form', alpha = 'cont')
