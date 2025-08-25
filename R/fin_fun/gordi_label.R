gordi_label <- function(pass,
                        what = c('species', 'sites', 'predictor'), #still need to add predictor, choose type of label
                        label_colour = '', #colour for non shortcut labels
                        label = '', #column with label names
                        shortcut = '', #creates shortcuts
                        shortcut_colour = '', #shortcut colour
                        shortcut_length = 3,
                        size = 3.9,
                        scaling_coefficient = 0.9,
                        nudge_x = 0,
                        nudge_y = 0,
                        max.overlaps = 10,
                        repel_label = FALSE){ #if TRUE -> geom_text_repel
  
  
  if (missing(what)){
    stop("'what' is missing. Please specify 'what' = 'species'|'sites'|'predictor'")
  }
  # label coerced to value label
  what <- match.arg(what)
  
  #check whether there is already a plot
  if (is.null(pass$plot)) warning('No plot yet, draw it first!') 
  p <- pass$plot
  
  #function removing label layers if they have been used before. Makes sure no duplicate labels
  remove_label_layers_for <- function(plot, axis_x_col) { #axis_x_col is either Axis_spe1 or Axis_spe2
    plot$layers <- discard(plot$layers, function(ly) {
      #detects label layer
      is_lab <- inherits(ly$geom, "GeomText") || inherits(ly$geom, "GeomLabel") || inherits(ly$geom, "GeomTextRepel")
      #whether label belongs to the coordinate system axis_x_col
      has_axis <- if (!is.null(ly$data)) {
        axis_x_col %in% names(ly$data)
      } else {
        # if the layer inherits data from another function, the data would be NULL, we ask ggplot to build the layer data
        tryCatch(axis_x_col %in% names(ly$layer_data(NULL)), error = function(e) FALSE)
      }
      is_lab && has_axis   # discard layer if TRUE
    })
    plot
  }
  
  #species labels
  
  if (what == 'species'){
    
    spe_df <- bind_cols(pass$species_names, pass$species_scores) #species data frame
    #join traits dataset if it exists in the pass list
    if (!is.null(pass$traits)){
      spe_df <- spe_df|>
        left_join(pass$traits, by = join_by(!!sym(names(spe_df)[1]) == !!sym(names(pass$traits)[1])))
    }
    
    #remove previously added species labels
    p <- remove_label_layers_for(p, "Axis_spe1")
    
    #default labeling column is the first column of spe_df
    text_col <- names(spe_df)[1]
    
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
        mutate(short_name = ifelse(has_sub, short_sub, short_non))
      text_col <- 'short_name' } #by default text_col is the first column of spe_df, if shortcut is defined, text_col is short_name
    
    #whether colour is mapped by a column or is a constant (both labels and shortcuts)
    
    map_shortcut_colour <- !identical(shortcut_colour, '') && has_name(spe_df, shortcut_colour)
    const_shortcut_colour <- !identical(shortcut_colour, '') && !map_shortcut_colour && (grepl("^#(?:[A-Fa-f0-9]{6}|[A-Fa-f0-9]{3})$", shortcut_colour) || shortcut_colour %in% grDevices::colours())|| (is.character(shortcut_colour) && shortcut_colour %in% palette()) || (is.numeric(shortcut_colour) && shortcut_colour %in% seq_along(palette()))
    
    map_label_colour <- !identical(label_colour, '') && has_name(spe_df, label_colour)
    const_label_colour <- !identical(label_colour, '') && !map_label_colour && (grepl("^#(?:[A-Fa-f0-9]{6}|[A-Fa-f0-9]{3})$", label_colour) || label_colour %in% grDevices::colours())|| (is.character(label_colour) && label_colour %in% palette()) || (is.numeric(label_colour) && label_colour %in% seq_along(palette()))
    
    if (map_shortcut_colour || map_label_colour){
      #if shortcut-colour or label_colour start a new colour scale
      col_var <- if (map_shortcut_colour) shortcut_colour else label_colour #col_var is shortcut_colour when map_shortcut_colour is true else its label_colour
      p <- p + ggnewscale::new_scale_colour()
      mapping <- aes(Axis_spe1, Axis_spe2, label = !!sym(text_col), colour = !!sym(col_var))
      
      if (isTRUE(repel_label)) p <- p + geom_text_repel(data = spe_df, mapping = mapping, size = size, nudge_x = nudge_x, nudge_y = nudge_y, max.overlaps = max.overlaps)
      else p <- p + geom_text(data = spe_df, mapping = mapping, size = size, nudge_x = nudge_x, nudge_y = nudge_y)
    } #constant colour scale
    else {
      col_const <- if (const_label_colour) label_colour
      else if (const_shortcut_colour) shortcut_colour
      else 'black'
      mapping <- aes(Axis_spe1, Axis_spe2, label = !!sym(text_col))
      if (isTRUE(repel_label)) p <- p + geom_text_repel(data = spe_df, mapping = mapping, colour = col_const, size = size, nudge_x = nudge_x, nudge_y = nudge_y, max.overlaps = max.overlaps)
      else                     p <- p + geom_text(data = spe_df, mapping = mapping, colour = col_const, size = size, nudge_x = nudge_x, nudge_y = nudge_y)
      
    }
  }
  
  #site labels
  
  if (what == 'sites') {
    
    site_df <- bind_cols(pass$env, pass$site_scores) #creates site_df dataframe
    
    #removes previous label layers
    p <- remove_label_layers_for(p, "Axis_site1")
    
    #which column to use for labeling, if column name provide it use it otherwise take the first column
    labcol <- if (!identical(label, '') && label %in% names(site_df)){
      label
    } else {names(site_df)[1]}
    
    #mapped and constant site label colour
    map_label_colour <- !identical(label_colour, '') && has_name(site_df, label_colour)
    const_label_colour <- !identical(label_colour, '') && !map_label_colour && (grepl("^#(?:[A-Fa-f0-9]{3}|[A-Fa-f0-9]{6})$", label_colour) || label_colour %in% grDevices::colours())|| (is.character(label_colour) && label_colour %in% palette()) || (is.numeric(label_colour) && label_colour %in% seq_along(palette()))
    
    if (map_label_colour){
      #new colour scale
      p <- p + ggnewscale::new_scale_colour()
      mapping <- aes(Axis_site1, Axis_site2, label = !!sym(labcol), colour = !!sym(label_colour))
      if (isTRUE(repel_label)) 
        p <- p + geom_text_repel(data = site_df, mapping = mapping, size = size, nudge_x = nudge_x, nudge_y = nudge_y, max.overlaps = max.overlaps)
      else p <- p + geom_text(data = site_df, mapping = mapping, size = size, nudge_x = nudge_x, nudge_y = nudge_y)
    } else {
      col_const <- if (const_label_colour) label_colour else 'black' #if not specified label_colour it will be black
      mapping <- aes(Axis_site1, Axis_site2, label = !!sym(labcol))
      if (isTRUE(repel_label)) p <- p + geom_text_repel(data = site_df, mapping = mapping, colour = col_const, size = size, nudge_x = nudge_x, nudge_y = nudge_y, max.overlaps = max.overlaps)
      else p <- p + geom_text(data = site_df, mapping = mapping, colour = col_const, size = size, nudge_x = nudge_x, nudge_y = nudge_y)
    }
  }
  if (what == 'predictor') {
    
    pred_df <- bind_cols(pass$predictor_scores, pass$predictor_names)
    p <- remove_label_layers_for(p, "Axis_pred1")
    labcol <- if(!identical(label, '') && label == names(pred_df)[3]){
      label}  else {names(pred_df)[3]}
    
    p_build <- ggplot_build(p)
    
    plot_range <- c(xmin_plot = p_build$layout$panel_params[[1]]$x.range[1],
                    xmax_plot = p_build$layout$panel_params[[1]]$x.range[2],
                    ymin_plot = p_build$layout$panel_params[[1]]$y.range[1],
                    ymax_plot = p_build$layout$panel_params[[1]]$y.range[2])
    
    predictor_range <- c(xmin_pred = min(pred_df[,1]),
                         xmax_pred = max(pred_df[,1]),
                         ymin_pred = min(pred_df[,2]),
                         ymax_pred = max(pred_df[,2]))
    
    coef <- (max(abs(plot_range)) / max(abs(predictor_range))) * scaling_coefficient
    
    map_label_colour <- !identical(label_colour, '') && has_name(pred_df, label_colour)
    const_label_colour <- !identical(label_colour, '') && !map_label_colour && (grepl("^#(?:[A-Fa-f0-9]{3}|[A-Fa-f0-9]{6})$", label_colour) || label_colour %in% grDevices::colours())|| (is.character(label_colour) && label_colour %in% palette()) || (is.numeric(label_colour) && label_colour %in% seq_along(palette()))
    if (map_label_colour){
      p <- p + ggnewscale::new_scale_colour()
      mapping <- aes(Axis_pred1*coef, Axis_pred2*coef, label = !!sym(labcol), colour = !!sym(label_colour))
      if(isTRUE(repel_label))
        p <- p + geom_text_repel(data = pred_df, mapping = mapping, size = size, nudge_x = nudge_x, nudge_y = nudge_y, max.overlaps = max.overlaps)
      else p <- p + geom_text(data = pred_df, mapping = mapping, size = size, nudge_x = nudge_x, nudge_y = nudge_y)
    } else {
      col_const <- if (const_label_colour) label_colour else 'black'
      mapping <- aes(Axis_pred1*coef, Axis_pred2*coef, label = !!sym(labcol))
      if (isTRUE(repel_label))
        p <- p + geom_text_repel(data = pred_df, mapping = mapping, size = size, colour = col_const, nudge_x = nudge_x, nudge_y = nudge_y, max.overlaps = max.overlaps)
      else p <- p + geom_text(data = pred_df, mapping = mapping, colour = col_const, size = size, nudge_x = nudge_x, nudge_y = nudge_y)
    }
  }
  
  pass$plot <- p
  return(pass)
}

