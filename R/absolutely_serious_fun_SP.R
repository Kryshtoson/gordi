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
m <- metaMDS(spe, k = 3)

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


#' =====================================================================
#' arguments: 
#' -> model
#' -> headers 


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
    inherits(m, 'capscale') & !is.null(m$call$distance) ~ 'db-RDA',
    inherits(m, 'rda') & is.null(m$call$distance) & is.null(m$CCA) ~ 'PCA', inherits(m, 'rda') ~ 'RDA',
    inherits(m, 'cca') & is.null(m$call$distance) & is.null(m$CCA) ~ 'CA', inherits(m, 'cca') ~ 'CCA',
    inherits(m, 'decorana') ~ 'DCA',
    inherits(m, 'metaMDS') ~ 'NMDS',
    TRUE ~ paste(class(m), collapse = '/') # writes just one output
  )
  
  # create pass object 
  pass <- list(
    m = m,
    explained_variation = if (type %in% c('DCA', 'NMDS')) {NA} else {m$CA$eig/m$tot.chi},
    site_scores = as_tibble(as.data.frame(scores(m, scaling = scaling, display = 'sites', choices = choices, correlation = correlation, hill = hill))),
    species_scores = as_tibble(as.data.frame(scores(m, scaling = scaling, display = 'species', choices = choices, correlation = correlation, hill = hill))),
    env = env,
    traits = traits,
    choices = choices,
    type = type,
    species_names = as_tibble(as.data.frame(scores(m, scaling = scaling, display = 'species', choices = choices, correlation = correlation, hill = hill)), rownames = 'species_names')[1]
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
      
      
      
      #' accounting for labeling
      if (label != '') {
        if (repel_label) {
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


