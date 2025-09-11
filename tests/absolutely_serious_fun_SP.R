#' =====================================================================
#' test run input
library(ggrepel)
library(tidyverse)
library(vegan)
library(readxl)


# Import -----------------------------------------------------------------------
spe <- read_csv('data/schrankogel/schrankogel_spe.csv')[-1] |> 
  log1p()

env <- read_csv("data/schrankogel/schrankogel_env.csv") |>
  mutate(group = elevation > 2500)

trait <- read_xlsx('data/Life_form.xlsx')|>
  select(-SeqID)|>
  pivot_longer(cols = -FloraVeg.Taxon, names_to = 'form', values_to = 'value') |> 
  filter(!value == 0) |> 
  distinct(FloraVeg.Taxon, .keep_all = T) |> 
  mutate(cont = rep(1:50, length.out = n()))




# Ordination -------------------------------------------------------------------

# PCA
m <- rda(spe ~ 1) 
class(m)
is.null(m$CCA) # TRUE = nema prediktor -> unconstrained
is.null(m$call$distance) # TRUE = neni tam distance


# RDA
m <- rda(spe ~ elevation, data = env)
class(m)
is.null(m$CCA) # FALSE = ma prediktor -> constrained
is.null(m$call$distance) # TRUE = neni tam distance


# CA
m <- cca(spe ~ 1) 
class(m)
is.null(m$CCA) # TRUE, unconstrained
# is.null(m$call$distance) v tomto pripade nema smysl


# CCA
m <- cca(spe ~ elevation, data = env)
class(m)
is.null(m$CCA) # FALSE - constrained


# PCoA
m <- capscale(spe ~ 1, distance = 'bray')
class(m)
is.null(m$CCA) # TRUE - unconstrained
is.null(m$call$distance) # FALSE - distance-based


# db-RDA
m <- capscale(spe ~ elevation, distance = 'bray', data = env)
class(m)
is.null(m$CCA) # FALSE - constrained
is.null(m$call$distance) # FALSE - distance-based

# DCA
spe_strip <- spe[,colSums(spe!=0)>5]
spe_upt <- spe_strip[rowSums(spe_strip) >0, ]
m <- decorana(spe_upt)
plot(m)

# NMDS
#m <- metaMDS(spe, k = 3)

as_tibble(as.data.frame(scores(m, scaling = 'symm', display = 'sites', choices = 1:2, correlation = F, hill = T)))
as_tibble(as.data.frame(scores(m, scaling = 'symm', display = 'species', choices = 1:2, correlation = F, hill = T)))
as_tibble(as.data.frame(scores(m, scaling = 'symm', display = 'species', choices = 1:2, correlation = F, hill = T)), rownames = 'species_names')[1]


# 
# case_when(
#   inherits(m, 'capscale') & !is.null(m$call$distance) & is.null(m$CCA) ~ 'PCoA (capscale)', #via capscale
#   inherits(m, 'capscale') & !is.null(m$call$distance) ~ 'db-RDA (capscale)',
#   inherits(m, 'rda') & !is.null(m$call$distance) & is.null(m$CCA) ~ 'PCoA (rda)', #if rda and distance but no constrainned variable PCoA
#   inherits(m, 'rda') & !is.null(m$call$distance) ~ 'db-RDA (by RDA argument)',
#   inherits(m, 'rda') & is.null(m$call$distance) & is.null(m$CCA) ~ 'PCA', inherits(m, 'rda') ~ 'RDA constrained',
#   inherits(m, 'cca') & is.null(m$call$distance) & is.null(m$CCA) ~ 'CA', inherits(m, 'cca') ~ 'CCA',
#   TRUE ~ paste(class(m), collapse = '/') # writes just one output
# )

#' cca {CCA, partial CCA, CA, partial CA}
#' capscale {dbRDA, partial dbRDA, PCoA, partial PCoA}
#' rda {RDA, partial RDA, PCA, partial PCA}
#' ...
#' decorana

browseVignettes()
#' =====================================================================
#' arguments: 
#' -> model
#' -> headers 


# gordi_read()

gordi_read <- function(m,
                       env = NULL,
                       traits = NULL,
                       choices = 1:2,
                       scaling = 'symm',
                       correlation = F,
                       hill = F) {
  
  # type of ordination 
  type <- case_when(
    inherits(m, 'capscale') & !is.null(m$call$distance) & is.null(m$CCA) ~ 'PCoA', #via capscale
    inherits(m, 'capscale') & !is.null(m$call$distance)                  ~ 'db-RDA',
    inherits(m, 'rda') & is.null(m$call$distance) & is.null(m$CCA)       ~ 'PCA',
    inherits(m, 'rda')                                                   ~ 'RDA',
    inherits(m, 'cca') & is.null(m$call$distance) & is.null(m$CCA)       ~ 'CA',
    inherits(m, 'cca')                                                   ~ 'CCA',
    inherits(m, 'decorana')                                              ~ 'DCA',
    inherits(m, c('metaMDS', 'monoMDS'))                                 ~ 'NMDS',
    TRUE ~ paste(class(m), collapse = '/') # writes just one output
  )
  
  # create pass object 
  pass <- list(
    m = m,
    explained_variation = m$CA$eig/m$tot.chi,
    site_scores = as_tibble(as.data.frame(scores(m, scaling = scaling, choices = choices, correlation = correlation, hill = hill)$sites)),
    species_scores = as_tibble(as.data.frame(scores(m, scaling = scaling, choices = choices, correlation = correlation, hill = hill)$species)),
    predictor_scores = as_tibble(as.data.frame(scores(m, scaling = scaling, choices = choices, correlation = correlation, hill = hill)$biplot)),
    env = env,
    traits = traits,
    choices = choices,
    type = type,
    axis_names = names(as_tibble(as.data.frame(scores(m, scaling = scaling, choices = choices, correlation = correlation, hill = hill)$sites))),
    species_names = as_tibble(as.data.frame(scores(m, scaling = scaling, choices = choices, correlation = correlation, hill = hill)$species), rownames = 'species_names')[1],
    predictor_names = as_tibble(as.data.frame(scores(m, scaling = scaling, choices = choices, correlation = correlation, hill = hill)$biplot), rownames = 'predictor_names')[1]
  )
  
  # Return pass object
  return(pass)
  
}




o <- gordi_read(m)

o |> 
  gordi_species()



#
#as_tibble(as.data.frame(scores(m, scaling = 'symm', choices = 1:2, correlation = T, hill = T)$sites), rownames = 'species_names')[1]

m <- rda(spe ~ elevation, distance = 'bray', data = env)
m <- rda(log1p(spe)) #pca or differently rda(spe ~ 1)
m <- rda(sqrt(spe) ~ elevation, data = env) #rda constrained
m <- capscale(spe ~ 1, distance = 'bray') #pcoa
m <- capscale(spe ~ elevation, distance = 'bray', data = env) #db-rda
m <- cca(spe) #ca
m <- cca(spe ~ elevation, data = env) # cca




#
gordi_sites <- function(pass, label = '', fill = '', alpha = '', stroke = '', shape = '', size = '', colour = '', repel_label = T) {
  
  #' axis names
  names(pass$site_scores) <- paste0("Axis_site", 1:2)
  names(pass$species_scores) <- paste0("Axis_spe", 1:2)
  names(pass$predictor_scores) <- paste0("Axis_pred", 1:2)
  
  #ordination axis labels
  if(pass$type %in% c('DCA', 'NMDS')) {actual_labs <- paste0(pass$axis_names)} else 
  {actual_labs <- paste0(pass$axis_names, " (", round(pass$explained_variation[choices]*100, 2), "%)")}
  
  
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


#

gordi_read(m, env, trait, choices = 3:4) |> 
  gordi_sites(colour = 'annual_temperature') |> 
  gordi_species(label = 'species_names', colour = 'form')


# gordi_species()

# arguments
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
  names(pass$predictor_scores) <- paste0("Axis_pred", 1:2)
  
  ### axis labs
  if(pass$type %in% c('DCA', 'NMDS')) {actual_labs <- paste0(pass$axis_names)} else 
    {actual_labs <- paste0(pass$axis_names, " (", round(pass$explained_variation[choices]*100, 2), "%)")}
  
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




# gordi_read(m, env) |> 
#   gordi_species(label = 'species_names')

gordi_read(m, env, trait) |> 
  gordi_sites(size = 'elevation') |> 
  gordi_species(colour = 'cont')

gordi_read(m, env, trait) |> 
  gordi_species(colour = 'form', size = 'cont', shape = 17, label = 'species_names', linewidth = 2)
 
gordi_read(m, env, trait) |> 
  gordi_sites(colouring = 'group', size = 'elevation') |> 
  gordi_species(label = 'species_names', colouring = 'form', size = 'cont')


gordi_read(m, env, trait) |> 
  gordi_species(colour = 'cont', label = 'species_names', alpha = 1, , arrow_size = 0.1) |> 
  gordi_sites(shape = 21, size = 'elevation', fill = 'group') 

gordi_read(m) |> 
  gordi_species(shortcut = 'upper.upper') |>  # Brano musi jeste opravit shortcuts
  gordi_predict()

#' gordi_predict dodelat !!!
#' sipecku uz kresli, ale chtelo by to nejak skalovat podle delky osy x (asi na delsi reseni)
#' zase klasicky, aby umel barvit dynamicky i staticky
#' zatim ale jenom pro kontinualni prediktory

gordi_predict <- function(
    pass,
    label = '',
    colour = '',
    alpha = '',
    arrow_size = '',
    linewidth = '',
    linetype = '',
    repel_label = T) {
  
  
  ### ordination axis labels
  if(pass$type %in% c('DCA', 'NMDS')) {actual_labs <- paste0(pass$axis_names)} else 
  {actual_labs <- paste0(pass$axis_names, " (", round(pass$explained_variation[choices]*100, 2), "%)")}
  
  
  ### axis names used in spe_df
  names(pass$species_scores) <- paste0("Axis_spe", 1:2)
  names(pass$site_scores) <- paste0("Axis_site", 1:2)
  names(pass$predictor_scores) <- paste0("Axis_pred", 1:2)
  
  
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
  pred_df <- bind_cols(pass$predictor_scores, pass$predictor_names)
  
  
  ### Detect mapped vs constant aesthetics
    # colour
    #' if colour != "" AND ALSO colour represents a colname present in spe_df, then map_colour is TRUE, otherwise is FALSE
      map_colour <- !identical(colour, '') && has_name(pred_df, colour) 
    # if map_colour is FALSE AND ALSO the thing inputed in arguments is a HEX code or is included in colours() or in palette() (word or number), then use it as const_colour
      const_colour <- !map_colour && (grepl("^#(?:[A-Fa-f0-9]{6}[A-Fa-f0-9]{3})$", colour) || colour %in% grDevices::colours()) || (is.character(colour) && colour %in% palette()) || (is.numeric(colour) && colour %in% seq_along(palette()))
    # alpha
      map_alpha <- !identical(alpha, '') && has_name(pred_df, alpha)
      const_alpha <- !map_alpha && is.numeric(alpha)
    # linetype
      map_linetype <- !identical(linetype, '') && has_name(pred_df, linetype)
      const_linetype <- !map_linetype && (is.numeric(linetype) || !identical(linetype,''))
    # linewidth
      map_linewidth <- !identical(linewidth, '') && has_name(pred_df, linewidth)
      const_linewidth <- !map_linewidth && is.numeric(linewidth)
   

  ### Prepare aes arguments for geom_segment()
    # Start with fixed x/y for the base (0,0) and end at the species scores
      aes_args_segment <- list(
        x = 0, y = 0,
        xend = quote(Axis_pred1),
        yend = quote(Axis_pred2)
      )
      
    if(map_colour) aes_args_segment$colour <- sym(colour)
    if(map_alpha) aes_args_segment$alpha <- sym(alpha)
    if(map_linewidth) aes_args_segment$linewidth <- sym(linewidth)
    if(map_linetype) aes_args_segment$linetype <- sym(linetype)
  
      
  ### Prepare constant arguments for geom_point() (mapped first, then defaults if nothing)
      const_args_segment <- list()
      
    # Add constant arguments for geom_segment() if not mapped
      # colour
      if(!map_colour){
        if(!identical(colour, '')) {const_args_segment$colour <- colour} else {const_args_segment$colour <- 2}}
      # alpha
      if(!map_alpha){
        if(!identical(alpha, '')) {const_args_segment$alpha <- alpha} else {const_args_segment$alpha <- 1}}
      # linetype
      if(!map_linetype){
        if(!identical(linetype, '')) {const_args_segment$linetype <- linetype} else {const_args_segment$linetype <- 1}}
      # linewidth
      if(!map_linewidth){
        if(!identical(linewidth, '')) {const_args_segment$linewidth <- linewidth} else {const_args_segment$linewidth <- 0.7}}
      # arrow_size (does not make sense to have map_arrow)
      const_args_segment$arrow <- arrow(
        length = unit(
          if (!identical(arrow_size, '')) as.numeric(arrow_size) else 0.3, "cm"))
  
      
      
      
  ### add to plot  
    p <- p + do.call(geom_segment,
                   c(list(data = pred_df,
                          mapping = do.call(aes, aes_args_segment)),
                     const_args_segment)) +
      geom_text_repel(data = pred_df, aes(x = Axis_pred1, Axis_pred2, label = predictor_names), colour = 2)
  
  ### save plot
  pass$plot <- p
  
  
  return(pass)
}

m <- capscale(spe ~ elevation + slope, data = env)
m <- cca(spe ~ elevation + slope, data = env)

gordi_read(m, env, trait, scaling = 'species', correlation = T, hill = F) |>
  gordi_species(shortcut = 'upperupper', fill = 'cont', shape = 21, colour = 'black') |> 
  gordi_predict(label = 'predictor_names',
                alpha = 0.8, linetype = 1, linewidth = 1)

gordi_read(m, env, trait, scaling = 'sites', correlation = T) |>
  gordi_sites(size = 'slope', label = 'logger_ID', fill = 'elevation', shape = 21) |> 
  gordi_predict(label = 'predictor_names',
                alpha = 0.8, linetype = 1, linewidth = 1)

gordi_read(m, env, trait, scaling = 'symmetric', correlation = T) |>
  gordi_sites(size = 'slope', label = 'logger_ID', fill = 'elevation', shape = 21) |> 
  gordi_predict(label = 'predictor_names',
                alpha = 1, linetype = 1, linewidth = 1)


# SCALING COEFFICIENT ---------------------

m <- capscale(spe ~ elevation + slope, data = env)

o <- gordi_read(m, scaling = 'sites')

# default predictor scores
ggplot(o$site_scores) +
  geom_segment(aes(x = 0, y = 0, xend = CAP1, yend = CAP2),
               arrow = arrow(length = unit(0.3, 'cm')),
               colour = 4,
               alpha = 0.7) +
  geom_segment(data = o$predictor_scores,
               aes(x = 0, y = 0, xend = CAP1, yend = CAP2),
               arrow = arrow(length = unit(0.3, 'cm')),
               colour = 2,
               alpha = 1) 

# wo predictor
p <- ggplot(o$site_scores) +
  geom_segment(aes(x = 0, y = 0, xend = CAP1, yend = CAP2),
               arrow = arrow(length = unit(0.3, 'cm')),
               colour = 4,
               alpha = 0.7) 


# min and max species scores on x and y axis
sco_x_min <- min(o$species_scores[1]) # min na ose x
sco_x_max <- max(o$species_scores[1]) # max na ose x

sco_y_min <- min(o$species_scores[2]) # min na ose y
sco_y_max <- max(o$species_scores[2]) # max na ose y

sco_range <- c(sco_x_min = min(o$species_scores[1]),
               sco_x_max = max(o$species_scores[1]),
               sco_y_min = min(o$species_scores[2]),
               sco_y_max = max(o$species_scores[2]))

abs(sco_range)

# min and max plot frame coordinates on x and y
p_build <- ggplot_build(p)

plot_x_min <- p_build$layout$panel_params[[1]]$x.range[1] # min na ose x
plot_x_max <- p_build$layout$panel_params[[1]]$x.range[2] # max na ose x

plot_y_min <- p_build$layout$panel_params[[1]]$y.range[1] # min na ose y
plot_y_max <- p_build$layout$panel_params[[1]]$y.range[2] # max na ose y

plot_range <- c(plot_x_min = p_build$layout$panel_params[[1]]$x.range[1],
                plot_x_max = p_build$layout$panel_params[[1]]$x.range[2],
                plot_y_min = p_build$layout$panel_params[[1]]$y.range[1],
                plot_y_max = p_build$layout$panel_params[[1]]$y.range[2])

abs(plot_range)

# min value
b <- max(abs(sco_range)) # x max
a <- max(abs(plot_range)) # x max

(a/b) * 0.8
(a / b)


# default predictor scores
ggplot() +
  geom_segment(data = bind_cols(o$site_scores, env),
               aes(x = 0, y = 0, xend = CAP1, yend = CAP2, colour = elevation),
               arrow = arrow(length = unit(0.3, 'cm')),
               alpha = 0.6) +
  geom_text_repel(data = bind_cols(o$site_scores, env),
                  aes(x = CAP1, y = CAP2,
                      label = as.character(logger_ID)),
                  colour = 'red') +
  scale_colour_viridis_c() +
  ggnewscale::new_scale_color() +
  geom_segment(data = bind_cols(o$predictor_scores, o$predictor_names),
               aes(x = 0, y = 0, xend = CAP1 * b, yend = CAP2 * b, colour = predictor_names),
               arrow = arrow(length = unit(0.3, 'cm')),
               alpha = 1) +
  scale_colour_manual(values = c('red', 'purple')) + NULL
  
  


o |> 
  gordi_species()
