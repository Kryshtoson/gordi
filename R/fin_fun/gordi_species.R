gordi_species <- function(pass,
                          label = '',
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
                          shortcut = '',
                          shortcut_length = 3, 
                          shortcut_colour = '',
                          repel_label = T) {
  
  ### axis names used in spe_df 
  names(pass$species_scores) <- paste0("Axis_spe", 1:2)
  names(pass$site_scores) <- paste0("Axis_site", 1:2)
  
  
  ### ordination types -> later used in axis labels 
  if (pass$type == 'CCA') {actual_labs <- paste0("CCA", pass$choices, " (", round(pass$explained_variation[1:2]*100, 2), '%)')} 
  else if (pass$type == 'CA') {actual_labs <- paste0("CA", pass$choices, " (", round(pass$explained_variation[1:2]*100, 2), '%)')}
  else if (pass$type == 'PCA') {actual_labs <- paste0("PCA", pass$choices, " (", round(pass$explained_variation[1:2]*100, 2), '%)')}
  else if(pass$type == 'PCoA') {actual_labs <- paste0("PCoA", pass$choices, " (", round(pass$explained_variation[1:2]*100, 2), '%)')}
  else if(pass$type == 'db-RDA') {actual_labs <- paste0("db-RDA", pass$choices, " (", round(pass$explained_variation[1:2]*100, 2), '%)')}
  else if(pass$type == 'RDA') {actual_labs <- paste0("RDA", pass$choices, " (", round(pass$explained_variation[1:2]*100, 2), '%)')}
  else if (pass$type == 'DCA') {actual_labs <- paste0('DCA', pass$choices)}
  else if (pass$type == 'NMDS') {actual_labs <- paste0('NMDS', pass$choices)} 
  
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
  
  
  ### Detect mapped vs constant aesthetics
  
  ### arguments working in both, points and arrows (colour, alpha)
  # colour
  #' if colour != "" AND ALSO colour represents a colname present in spe_df, then map_colour is TRUE, otherwise is FALSE
  map_colour <- !identical(colour, '') && has_name(spe_df, colour) 
  # if map_colour is FALSE AND ALSO the thing inputed in arguments is a HEX code or is included in colours() or in palette() (word or number), then use it as const_colour
  const_colour <- !map_colour && (grepl("^#(?:[A-Fa-f0-9]{6}[A-Fa-f0-9]{3})$", colour) || colour %in% grDevices::colours()) || (is.character(colour) && colour %in% palette()) || (is.numeric(colour) && colour %in% seq_along(palette()))
  # alpha
  map_alpha <- !identical(alpha, '') && has_name(spe_df, alpha)
  const_alpha <- !map_alpha && is.numeric(alpha)
  
  ### arguments working only in geom_point
  # size
  map_size <- !identical(size, '') && has_name(spe_df, size)
  const_size <- !map_size && is.numeric(size)
  # shape
  map_shape <- !identical(shape, '') && has_name(spe_df, shape)
  const_shape <- !map_shape && is.numeric(shape)
  # fill
  map_fill <- !identical(fill, '') && has_name(spe_df, fill)
  const_fill <- !map_fill && (grepl("^#(?:[A-Fa-f0-9]{6}[A-Fa-f0-9]{3})$", fill) || fill %in% grDevices::colours())
  # stroke
  map_stroke <- !identical(stroke, '') && has_name(spe_df, stroke)
  const_stroke <- !map_stroke && is.numeric(stroke)
  
  ### arguments working only in geom_segment
  # linetype
  map_linetype <- !identical(linetype, '') && has_name(spe_df, linetype)
  const_linetype <- !map_linetype && (is.numeric(linetype) || !identical(linetype,''))
  # linewidth
  map_linewidth <- !identical(linewidth, '') && has_name(spe_df, linewidth)
  const_linewidth <- !map_linewidth && is.numeric(linewidth)
  
  
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
  if(map_size) aes_args_segment$linewidth <- sym(linewidth)
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
  
  
  
  #' Add the layer
  #' If linear ordination is used (PCA, RDA, PCoA, db-RDA), arrows are used
  #' if unimodal (CA, CCA, DCA, NMDS), points are used
  
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
  
  
  
  #'shortcuts
  if(!identical(shortcut, '')){
    #' split species  into individual tokens
    parts_list <- str_split(spe_df[[1]], '\\s')
    #' remove tokens like Sect., sect.... 
    rank_tokens <- c('sect\\.', 'Sect\\.', 'cf\\.')
    rx_drop <- regex(paste0('^(', paste(rank_tokens,  collapse = '|'), ')$'))
    parts_list <- lapply(parts_list, function (x) x[!str_detect(x, rx_drop)])
    #' detect subspecies and take epithet after it
    rx_sub <- regex('^(subsp\\.|ssp\\.)$')
    has_sub <- vapply(parts_list, function (x) any(str_detect(x, rx_sub)), logical (1))
    #'individual genus, species, subspecies
    genus <- map_chr(parts_list, 1)
    epithet <- map_chr(parts_list, 2)
    subsp <- vapply(parts_list, function (x) { i <- match(TRUE, str_detect(x, rx_sub))
    x[i + 1L]}, character(1))
    #'shortcuts based on shortcut_length
    gN <- str_sub(genus, 1, shortcut_length)
    sN <- str_sub(epithet, 1, shortcut_length)
    subN <- str_sub(subsp, 1, shortcut_length)
    
    if(shortcut == 'upper.lower') {
      short_non <- str_c(str_to_title(gN), str_to_lower(sN), sep = '.')
      short_sub <- str_c(str_to_title(gN), str_to_lower(sN), 'ssp', str_to_lower(subN), sep = '.')
    } else if(shortcut == 'lower.lower'){
      short_non <- str_c(str_to_lower(gN), str_to_lower(sN), sep = '.')
      short_sub <- str_c(str_to_lower(gN), str_to_lower(sN), 'ssp', str_to_lower(subN), sep = '.')
    } else if(shortcut == 'upper.upper'){
      short_non <- str_c(str_to_title(gN), str_to_title(sN), sep = '.')
      short_sub <- str_c(str_to_title(gN), str_to_title(sN), 'ssp', str_to_title(subN), sep = '.')
    } else if(shortcut == 'upperupper'){
      short_non <- str_c(str_to_title(gN), str_to_title(sN), sep = '')
      short_sub <- str_c(str_to_title(gN), str_to_title(sN), 'ssp', str_to_title(subN), sep = '')
    } else if(shortcut == 'upper_lower'){
      short_non <- str_c(str_to_title(gN), str_to_lower(sN), sep = '_')
      short_sub <- str_c(str_to_title(gN), str_to_lower(sN), 'ssp', str_to_lower(subN), sep = '_')
    } else if(shortcut == 'lower_lower'){
      short_non <- str_c(str_to_lower(gN), str_to_lower(sN), sep = '_')
      short_sub <- str_c(str_to_lower(gN), str_to_lower(sN), 'ssp', str_to_lower(subN), sep = '_')
    } else if(shortcut == 'upper_upper'){
      short_non <- str_c(str_to_title(gN), str_to_title(sN), sep = '_')
      short_sub <- str_c(str_to_title(gN), str_to_title(sN), 'ssp', str_to_title(subN), sep = '_')
    } else if(shortcut == 'upper*lower'){
      short_non <- str_c(str_to_title(gN), str_to_lower(sN), sep = '*')
      short_sub <- str_c(str_to_title(gN), str_to_lower(sN), 'ssp', str_to_lower(subN), sep = '*')
    } else if(shortcut == 'lower*lower'){
      short_non <- str_c(str_to_lower(gN), str_to_lower(sN), sep = '*')
      short_sub <- str_c(str_to_lower(gN), str_to_lower(sN), 'ssp', str_to_lower(subN), sep = '*')
    } else if(shortcut == 'upper*upper'){
      short_non <- str_c(str_to_title(gN), str_to_title(sN), sep = '*')
      short_sub <- str_c(str_to_title(gN), str_to_title(sN), 'ssp', str_to_title(subN), sep = '*')
    } else if(shortcut == 'upper-lower'){
      short_non <- str_c(str_to_title(gN), str_to_lower(sN), sep = '-')
      short_sub <- str_c(str_to_title(gN), str_to_lower(sN), 'ssp', str_to_lower(subN), sep = '-')
    } else if(shortcut == 'lower-lower'){
      short_non <- str_c(str_to_lower(gN), str_to_lower(sN), sep = '-')
      short_sub <- str_c(str_to_lower(gN), str_to_lower(sN), 'ssp', str_to_lower(subN), sep = '-')
    } else if(shortcut == 'upper-upper'){
      short_non <- str_c(str_to_title(gN), str_to_title(sN), sep = '-')
      short_sub <- str_c(str_to_title(gN), str_to_title(sN), 'ssp', str_to_title(subN), sep = '-')
    } else {
      warning("Unknown 'shortcut': ", shortcut, ' -> No short name created.')
    }
    
    
    #' creates tibble with short names if subspecies is non existent it takes short_non otherwise short_sub
    spe_df <- spe_df|>
      mutate(short_name = ifelse(has_sub, short_sub, short_non))}
  
  map_shortcut_colour <- !identical(shortcut_colour, '') && has_name(spe_df, shortcut_colour)
  const_shortcut_colour <- !identical(shortcut_colour, '') && !map_shortcut_colour && (grepl("^#(?:[A-Fa-f0-9]{6}|[A-Fa-f0-9]{3})$", shortcut_colour) || shortcut_colour %in% grDevices::colours())
  
  
  #' accounting for labeling
  if (identical(label, '') && !identical(shortcut, '')) {
    if (map_shortcut_colour){
      p <- p + ggnewscale::new_scale_colour()
    }
    map_args_text <- if (map_shortcut_colour){
      aes(Axis_spe1, Axis_spe2, label = short_name, colour = !!sym(shortcut_colour))
    } else {
      aes(Axis_spe1, Axis_spe2, label = short_name)
    }
    const_args_text <- list()
    if (!map_shortcut_colour){
      if(const_shortcut_colour){
        const_args_text$colour <- shortcut_colour
      } else {const_args_text$colour <- 'black'}
    }
    
    if (isTRUE(repel_label)) {
      p <- p + do.call(geom_text_repel, c(list(data = spe_df, mapping = map_args_text), const_args_text))
    } else {
      p <- p + do.call(geom_text, c(list(data = spe_df, mapping = map_args_text), const_args_text))
    }
  } else if (!identical(label, '') && identical (shortcut, '')){
    if (isTRUE(repel_label)){
      p <- p + geom_text_repel(data = spe_df, aes(Axis_spe1, Axis_spe2, label = !!sym(label)))
    } else {
      p <- p + geom_text(data = spe_df, aes(Axis_spe1, Axis_spe2, label = !!sym(label)))
    }
  }
  
  # More scales possibility (e.g. one colour in sites and other in species)
  p <- p + ggnewscale::new_scale_colour()
  p <- p + ggnewscale::new_scale('size')
  p <- p + ggnewscale::new_scale('shape')
  p <- p + ggnewscale::new_scale_fill()
  p <- p + ggnewscale::new_scale('alpha')
  p <- p + ggnewscale::new_scale('stroke')
  
  
  # Save plot into pass
  pass$plot <- p
  
  # Return pass
  return(pass)
}
  