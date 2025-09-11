#' =====================================================================
#' test run input
library(ggrepel)
library(tidyverse)
library(vegan)
library(readxl)
library(rlang)
library(purrr)
library(ggnewscale)

spe <- read_csv('data/schrankogel/schrankogel_spe.csv')[-1] |> 
  log1p()

env <- read_csv("data/schrankogel/schrankogel_env.csv") |>
  mutate(group = elevation > 2500)

traits <- read_xlsx('data/Life_form.xlsx') |>
  select(-SeqID) |>
  pivot_longer(cols = -FloraVeg.Taxon, names_to = 'form', values_to = 'value') |>
  filter(value == 1) |>
  distinct(FloraVeg.Taxon, .keep_all = T) %>%
  mutate(cont = 1:16382)


m <- rda(spe ~ 1) #pca
class(m)
m$PCoA
m$call$distance

m <- rda(spe ~ elevation, distance = 'bray', data = env)

m1 <- rda(log1p(spe)) #pca or differently rda(spe ~ 1)
class(m1)

m <- rda(sqrt(spe) ~ elevation, data = env) #rda constrained
class(m2)

m <- capscale(spe ~ 1, distance = 'bray') #pcoa
class(m3)
is.null(m3$CCA) #unconstrained T
!is.null(m3$call$distance)
m$CA$eig/m$tot.chi
m

m <- capscale(spe ~ elevation, distance = 'bray', data = env) #db-rda
class(m4)
is.null(m4$CCA) #constrained F

m <- cca(spe) #ca
class(m5)
is.null(m5$CCA) # unconstrained T

m <- cca(spe ~ elevation, data = env) # cca
class(m6)
is.null(m6$CCA) #connstrained F

m7 <- decorana(log1p(spe)) #dca
class(m7)

inherits(m, 'rda')
class(m)

# cca(spe~ elevation, data = env)
# m <- capscale(spe ~ 1, distance = 'bray')
# m <- capscale(spe ~ elevation + annual_temperature , distance = 'bray', data = env)
# 
# m <- decorana(spe)
# 
# as_tibble(scores(m, tidy = T))
# 
# scores(m)$sites
# 
# is.null(m$CCA[1])
# 
# class(m)

case_when(
  inherits(m, 'capscale') & !is.null(m$call$distance) & is.null(m$CCA) ~ 'PCoA (capscale)', #via capscale
  inherits(m, 'capscale') & !is.null(m$call$distance) ~ 'db-RDA (capscale)',
  inherits(m, 'rda') & !is.null(m$call$distance) & is.null(m$CCA) ~ 'PCoA (rda)', #if rda and distance but no constrainned variable PCoA
  inherits(m, 'rda') & !is.null(m$call$distance) ~ 'db-RDA (by RDA argument)',
  inherits(m, 'rda') & is.null(m$call$distance) & is.null(m$CCA) ~ 'PCA',
  inherits(m, 'rda') ~ 'RDA constrained',
  inherits(m, 'cca') & is.null(m$CCA) ~ 'CA',
  inherits(m, 'cca') ~ 'CCA',
  TRUE ~ paste(class(m), collapse = '/') #writes just one output
)
  
# cca {CCA, partial CCA, CA, partial CA}
#' capscale {dbRDA, partial dbRDA, PCoA, partial PCoA}
#' rda {RDA, partial RDA, PCA, partial PCA}
#' ...


#' =====================================================================
#' arguments: 
#' -> model
#' -> headers 


gordi_read <- function(m, env = NULL, traits = NULL, choices = 1:2, scaling = 'symm', correlation = F, hill = F){
  
 
  type <- case_when(
    inherits(m, 'capscale') & !is.null(m$call$distance) & is.null(m$CCA) ~ 'PCoA', #via capscale
    inherits(m, 'capscale') & !is.null(m$call$distance) ~ 'db-RDA',
   # inherits(m, 'rda') & !is.null(m$call$distance) & is.null(m$CCA) ~ 'PCoA (rda)', #if rda and distance but no constrainned variable PCoA
   # inherits(m, 'rda') & !is.null(m$call$distance) ~ 'db-RDA (by RDA argument)',
    inherits(m, 'rda') & is.null(m$call$distance) & is.null(m$CCA) ~ 'PCA',
    inherits(m, 'rda') ~ 'RDA',
    inherits(m, 'cca') & is.null(m$CCA) ~ 'CA',
    inherits(m, 'cca') ~ 'CCA',
    inherits(m, 'decorana') ~ 'DCA',
    inherits(m, 'metaMDS') ~ 'NMDS',
    TRUE ~ paste(class(m), collapse = '/') #writes just one output
  )
  
  pass <- list(
    m = m,
    explained_variation = if(type %in% c('DCA', 'NMDS')) {NA} else {m$CA$eig/m$tot.chi},
    site_scores = as_tibble(as.data.frame(scores(m, display = 'sites', scaling = scaling, choices = choices, correlation = correlation, hill = hill))),
    species_scores = as_tibble(as.data.frame(scores(m, display = 'species', scaling = scaling, choices = choices, correlation = correlation, hill = hill))),
    env = env,
    traits = traits,
    choices = choices,
    type = type,
    species_names = as_tibble(as.data.frame(scores(m, display = 'species', scaling = scaling, choices = choices, correlation = correlation, hill = hill)), rownames = 'species_names')[1]
  )
  
  
  return(pass)
  
}


gordi_read(m8, env, traits = traits, choices = 1:2)

m7$evals
m7$totchi

m7$rproj #sites scores
m7$cproj|>as_tibble() #species scores
scores(m7, choices = 1:2, scaling = 'symm', correlation = T, hill = T)
as_tibble(as.data.frame(scores(m7,display = 'species', scaling = 'symm', choices = 1:2, correlation = F, hill = F)), rownames = 'species')
gordi_read(m)|>
  gordi_sites()

m8 <- metaMDS(spe, k = 3)
class(m8)


m <- rda(spe ~ elevation, distance = 'bray', data = env)
m <- rda(log1p(spe)) #pca or differently rda(spe ~ 1)
m <- rda(sqrt(spe) ~ elevation, data = env) #rda constrained
m <- capscale(spe ~ 1, distance = 'bray') #pcoa
m <- capscale(spe ~ elevation, distance = 'bray', data = env) #db-rda
m <- cca(spe) #ca
m <- cca(spe ~ elevation, data = env) # cca
m7 <- decorana(spe) #dca

m7$evals

m7$rproj #sites scores
m7$cproj|>as_tibble() #species scores

scores(m7, 'species')|> as_tibble(rownames = 'spe_name')

m7$adotj
m7$call
# 
# if(o$type == 'CCA'){
#   actual_labs <- paste0("CCA", o$choices, " (", round(o$explained_variation[1:2]*100, 2), '%)')
# } else if(o$type == 'CA'){
#   actual_labs <- paste0("CA", o$choices, " (", round(o$explained_variation[1:2]*100, 2), '%)')
# } else if (o$type == 'PCA'){
#   actual_labs <- paste0("PCA", o$choices, " (", round(o$explained_variation[1:2]*100, 2), '%)')
# } else if(o$type %in% c('PCoA (capscale)', 'PCoA (rda)')){
#   actual_labs <- paste0("PCoA", o$choices, " (", round(o$explained_variation[1:2]*100, 2), '%)')
# } else if(o$type %in% c('db-RDA (rda)', 'db-RDA (capscale)')){
#   actual_labs <- paste0("db_RDA", o$choices, " (", round(o$explained_variation[1:2]*100, 2), '%)')
# } else if(o$type == 'RDA constrained'){
#   actual_labs <- paste0("RDA", o$choices, " (", round(o$explained_variation[1:2]*100, 2), '%)')
# } 



class(o$m)
names(o$site_scores) <- paste0('Axis', 1:2)
bind_cols(env, o$site_scores)


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
  
  p <- p + ggnewscale::new_scale_colour()
  p <- p + ggnewscale::new_scale_fill()
  
  site_df <- bind_cols(pass$env, pass$site_scores)
  
  #' accounting for labeling
  if (label != '') {
    if (repel_label) {
      p <- p + geom_text_repel(data = site_df, aes(Axis_site1, Axis_site2, label = !!sym(label)))
    } else {
      p <- p + geom_text(data = site_df, aes(Axis_site1, Axis_site2, label = !!sym(label)))
    }
  }

  
  
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
  
 
  
  #' accounting for indiviidual geom scales  
  #' p <- p + ggnewscale::new_scale("fill")
  #' p <- p + ggnewscale::new_scale_colour()
  p <- p + ggnewscale::new_scale("size") 
  p <- p + ggnewscale::new_scale("shape") 
  p <- p + ggnewscale::new_scale("alpha")
  p <- p + ggnewscale::new_scale("stroke") 
 
  
  pass$plot <- p
  
  return(pass)
}


gordi_species <- function(pass,
                          label = TRUE,
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
                          repel_label = FALSE) {
  
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
 
  p <- p + ggnewscale::new_scale_colour() 
  p <- p + ggnewscale::new_scale_fill()
  
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
  
  
  
  #' #'shortcuts
  #' if(!identical(shortcut, '')){
  #'   #' split species  into individual tokens
  #'   parts_list <- str_split(spe_df[[1]], '\\s')
  #'   #' remove tokens like Sect., sect.... 
  #'   rank_tokens <- c('sect\\.', 'Sect\\.', 'cf\\.')
  #'   rx_drop <- regex(paste0('^(', paste(rank_tokens,  collapse = '|'), ')$'))
  #'   parts_list <- lapply(parts_list, function (x) x[!str_detect(x, rx_drop)])
  #'   #' detect subspecies and take epithet after it
  #'   rx_sub <- regex('^(subsp\\.|ssp\\.)$')
  #'   has_sub <- vapply(parts_list, function (x) any(str_detect(x, rx_sub)), logical (1))
  #'   #'individual genus, species, subspecies
  #'   genus <- map_chr(parts_list, 1)
  #'   epithet <- map_chr(parts_list, 2)
  #'   subsp <- vapply(parts_list, function (x) { i <- match(TRUE, str_detect(x, rx_sub))
  #'   x[i + 1L]}, character(1))
  #'   #'shortcuts based on shortcut_length
  #'   gN <- str_sub(genus, 1, shortcut_length)
  #'   sN <- str_sub(epithet, 1, shortcut_length)
  #'   subN <- str_sub(subsp, 1, shortcut_length)
  #'   
  #'   if(shortcut == 'upper.lower') {
  #'     short_non <- str_c(str_to_title(gN), str_to_lower(sN), sep = '.')
  #'     short_sub <- str_c(str_to_title(gN), str_to_lower(sN), 'ssp', str_to_lower(subN), sep = '.')
  #'   } else if(shortcut == 'lower.lower'){
  #'     short_non <- str_c(str_to_lower(gN), str_to_lower(sN), sep = '.')
  #'     short_sub <- str_c(str_to_lower(gN), str_to_lower(sN), 'ssp', str_to_lower(subN), sep = '.')
  #'   } else if(shortcut == 'upper.upper'){
  #'     short_non <- str_c(str_to_title(gN), str_to_title(sN), sep = '.')
  #'     short_sub <- str_c(str_to_title(gN), str_to_title(sN), 'ssp', str_to_title(subN), sep = '.')
  #'   } else if(shortcut == 'upperupper'){
  #'     short_non <- str_c(str_to_title(gN), str_to_title(sN), sep = '')
  #'     short_sub <- str_c(str_to_title(gN), str_to_title(sN), 'ssp', str_to_title(subN), sep = '')
  #'   } else if(shortcut == 'upper_lower'){
  #'     short_non <- str_c(str_to_title(gN), str_to_lower(sN), sep = '_')
  #'     short_sub <- str_c(str_to_title(gN), str_to_lower(sN), 'ssp', str_to_lower(subN), sep = '_')
  #'   } else if(shortcut == 'lower_lower'){
  #'     short_non <- str_c(str_to_lower(gN), str_to_lower(sN), sep = '_')
  #'     short_sub <- str_c(str_to_lower(gN), str_to_lower(sN), 'ssp', str_to_lower(subN), sep = '_')
  #'   } else if(shortcut == 'upper_upper'){
  #'     short_non <- str_c(str_to_title(gN), str_to_title(sN), sep = '_')
  #'     short_sub <- str_c(str_to_title(gN), str_to_title(sN), 'ssp', str_to_title(subN), sep = '_')
  #'   } else if(shortcut == 'upper*lower'){
  #'     short_non <- str_c(str_to_title(gN), str_to_lower(sN), sep = '*')
  #'     short_sub <- str_c(str_to_title(gN), str_to_lower(sN), 'ssp', str_to_lower(subN), sep = '*')
  #'   } else if(shortcut == 'lower*lower'){
  #'     short_non <- str_c(str_to_lower(gN), str_to_lower(sN), sep = '*')
  #'     short_sub <- str_c(str_to_lower(gN), str_to_lower(sN), 'ssp', str_to_lower(subN), sep = '*')
  #'   } else if(shortcut == 'upper*upper'){
  #'     short_non <- str_c(str_to_title(gN), str_to_title(sN), sep = '*')
  #'     short_sub <- str_c(str_to_title(gN), str_to_title(sN), 'ssp', str_to_title(subN), sep = '*')
  #'   } else if(shortcut == 'upper-lower'){
  #'     short_non <- str_c(str_to_title(gN), str_to_lower(sN), sep = '-')
  #'     short_sub <- str_c(str_to_title(gN), str_to_lower(sN), 'ssp', str_to_lower(subN), sep = '-')
  #'   } else if(shortcut == 'lower-lower'){
  #'     short_non <- str_c(str_to_lower(gN), str_to_lower(sN), sep = '-')
  #'     short_sub <- str_c(str_to_lower(gN), str_to_lower(sN), 'ssp', str_to_lower(subN), sep = '-')
  #'   } else if(shortcut == 'upper-upper'){
  #'     short_non <- str_c(str_to_title(gN), str_to_title(sN), sep = '-')
  #'     short_sub <- str_c(str_to_title(gN), str_to_title(sN), 'ssp', str_to_title(subN), sep = '-')
  #'   } else {
  #'     warning("Unknown 'shortcut': ", shortcut, ' -> No short name created.')
  #'   }
  #' 
  #' 
  #' #' creates tibble with short names if subspecies is non existent it takes short_non otherwise short_sub
  #' spe_df <- spe_df|>
  #'   mutate(short_name = ifelse(has_sub, short_sub, short_non))}
  #' # 
  # map_shortcut_colour <- !identical(shortcut_colour, '') && has_name(spe_df, shortcut_colour)
  # const_shortcut_colour <- !identical(shortcut_colour, '') && !map_shortcut_colour && (grepl("^#(?:[A-Fa-f0-9]{6}|[A-Fa-f0-9]{3})$", shortcut_colour) || shortcut_colour %in% grDevices::colours())
  # 
  
  #' accounting for labeling
  # if (identical(label, '') && !identical(shortcut, '')) {
  #   if (map_shortcut_colour){
  #     p <- p + ggnewscale::new_scale_colour()
  #   }
  #   map_args_text <- if (map_shortcut_colour){
  #     aes(Axis_spe1, Axis_spe2, label = short_name, colour = !!sym(shortcut_colour))
  #   } else {
  #     aes(Axis_spe1, Axis_spe2, label = short_name)
  #   }
  #   const_args_text <- list()
  #   if (!map_shortcut_colour){
  #     if(const_shortcut_colour){
  #       const_args_text$colour <- shortcut_colour
  #     } else {const_args_text$colour <- 'black'}
  #   }
  #   
  #   if (isTRUE(repel_label)) {
  #     p <- p + do.call(geom_text_repel, c(list(data = spe_df, mapping = map_args_text), const_args_text))
  #   } else {
  #     p <- p + do.call(geom_text, c(list(data = spe_df, mapping = map_args_text), const_args_text))
  #   }
  # } else if (!identical(label, '') && identical (shortcut, '')){
  #   if (isTRUE(repel_label)){
  #     p <- p + geom_text_repel(data = spe_df, aes(Axis_spe1, Axis_spe2, label = !!sym(label)))
  #   } else {
  #     p <- p + geom_text(data = spe_df, aes(Axis_spe1, Axis_spe2, label = !!sym(label)))
  #   }
  # }
  # 
  if (isTRUE(label)){
    if (isTRUE(repel_label)){
      p <- p + geom_text_repel(data = spe_df, aes(Axis_spe1, Axis_spe2, label = species_names), colour = 'black') 
    } else {
      p <- p + geom_text(data = spe_df, aes(Axis_spe1, Axis_spe2, label = species_names), colour = 'black')
    }
  }
  
  
  # More scales possibility (e.g. one colour in sites and other in species)
  # p <- p + ggnewscale::new_scale_colour()
  p <- p + ggnewscale::new_scale('size')
  p <- p + ggnewscale::new_scale('shape')
  #p <- p + ggnewscale::new_scale_fill()
  p <- p + ggnewscale::new_scale('alpha')
  p <- p + ggnewscale::new_scale('stroke')
  
  
  # Save plot into pass
  pass$plot <- p
  
  
  # Return pass
  return(pass)
}


gordi_read(m, env, traits)|>
  gordi_sites(colour = 'elevation')|>
  gordi_species(colour = 'black')
  


gordi_label <- function(pass,
                        label = c('species', 'sites', 'predictor'), #still need to add predictor
                        label_colour = '',
                        label_name = '',
                        shortcut = '',
                        shortcut_colour = '',
                        shortcut_length = 3,
                        repel_label = FALSE){
  
  label <- match.arg(label)
    
  
  if (is.null(pass$plot)) warning('No plot yet, draw it first!') 
  p <- pass$plot
  
  #function removing label layers if they have been used before
  remove_label_layers_for <- function(plot, axis_x_col) {
    plot$layers <- discard(plot$layers, function(ly) {
      is_lab <- inherits(ly$geom, "GeomText") || inherits(ly$geom, "GeomLabel") || inherits(ly$geom, "GeomTextRepel")
      has_axis <- if (!is.null(ly$data)) {
        axis_x_col %in% names(ly$data)
      } else {
        tryCatch(axis_x_col %in% names(ly$layer_data(NULL)), error = function(e) FALSE)
      }
      is_lab && has_axis   # discard if TRUE
    })
    plot
  }
  
  if (label == 'species'){
    spe_df <- bind_cols(pass$species_names, pass$species_scores)
    if (!is.null(pass$traits)){
      spe_df <- spe_df|>
        left_join(pass$traits, by = join_by(!!sym(names(spe_df)[1]) == !!sym(names(pass$traits)[1])))
    }
    
    p <- remove_label_layers_for(p, "Axis_spe1")
    
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
      text_col <- 'short_name' }
 
  
  map_shortcut_colour <- !identical(shortcut_colour, '') && has_name(spe_df, shortcut_colour)
  const_shortcut_colour <- !identical(shortcut_colour, '') && !map_shortcut_colour && (grepl("^#(?:[A-Fa-f0-9]{6}|[A-Fa-f0-9]{3})$", shortcut_colour) || shortcut_colour %in% grDevices::colours())
  
  map_label_colour <- !identical(label_colour, '') && has_name(spe_df, label_colour)
  const_label_colour <- !identical(label_colour, '') && !map_label_colour && (grepl("^#(?:[A-Fa-f0-9]{6}|[A-Fa-f0-9]{3})$", label_colour) || label_colour %in% grDevices::colours())
  
  if (map_shortcut_colour || map_label_colour){
    col_var <- if (map_shortcut_colour) shortcut_colour else label_colour
    p <- p + new_scale_colour()
    mapping <- aes(Axis_spe1, Axis_spe2, label = !!sym(text_col), colour = !!sym(col_var))
    
    if (isTRUE(repel_label)) p <- p + geom_text_repel(data = spe_df, mapping = mapping)
    else p <- p + geom_text(data = spe_df, mapping = mapping)
  } else {
    col_const <- if (const_label_colour) label_colour
    else if (const_shortcut_colour) shortcut_colour
    else 'black'
    mapping <- aes(Axis_spe1, Axis_spe2, label = !!sym(text_col))
    if (isTRUE(repel_label)) p <- p + geom_text_repel(data = spe_df, mapping = mapping, colour = col_const)
    else                     p <- p + geom_text(data = spe_df, mapping = mapping, colour = col_const)
    
  }
  }
  if (label == 'sites') {
    site_df <- bind_cols(pass$env, pass$site_scores)
    
    p <- remove_label_layers_for(p, "Axis_site1")
    
    labcol <- if (!identical(label_name, '') && label_name %in% names(site_df)){
      label_name
    } else {names(site_df)[1]}
    
    map_label_colour <- !identical(label_colour, '') && has_name(site_df, label_colour)
    const_label_colour <- !identical(label_colour, '') && !map_label_colour && (grepl("^#(?:[A-Fa-f0-9]{3}|[A-Fa-f0-9]{6})$", label_colour) || label_colour %in% grDevices::colours())
    
    if (map_label_colour){
      p <- p + new_scale_colour()
      mapping <- aes(Axis_site1, Axis_site2, label = !!sym(labcol), colour = !!sym(label_colour))
      if (isTRUE(repel_label)) 
        p <- p + geom_text_repel(data = site_df, mapping = mapping)
      else p <- p + geom_text(data = site_df, mapping = mapping)
    } else {
      col_const <- if (const_label_colour) label_colour else 'black'
      mapping <- aes(Axis_site1, Axis_site2, label = !!sym(labcol))
      if (isTRUE(repel_label)) p <- p + geom_text_repel(data = site_df, mapping = mapping, colour = col_const)
      else p <- p + geom_text(data = site_df, mapping = mapping, colour = col_const)
    }
  }
  
  pass$plot <- p
  return(pass)
}

gordi_read(m, env, traits)|>
  gordi_species(label = F, colour = 'cont')|>
  gordi_colour(scale = 'continuous', family = 'viridis')|>
  gordi_label(label = 'species', shortcut = 'upper.upper', shortcut_colour = 'cont')|>
  gordi_colour(scale = 'continuous', family = 'viridis')
  gordi_sites(colour = 'elevation')|>
  gordi_label(label = 'sites', label_name = 'elevation', label_colour = 'elevation', repel_label = F)

gordi_read(m, env, traits)|>
  gordi_sites()|>
  gordi_species()|>
  gordi_label(label = 'species', shortcut = 'upper.upper', shortcut_colour = 'cont', repel_label = T)|>
  gordi_colour(scale = 'continuous', family = 'viridis')|>
  gordi_label(label = 'sites', label_colour = 'elevation', repel_label = T)|>
  gordi_colour(scale = 'continuous', family = 'viridis')
  
  gordi_label(label = 'sites')

gordi_colour <- function(pass,
                         scale = c('auto', 'continuous', 'discrete', 'binned'),
                         family = c('viridis', 'brewer', 'gradient', 'steps', 'manual', 'default'),
                         breaks = waiver(), name = waiver(), labels = waiver(),
                         limits = NULL, na.value = NA, guide = waiver(), trans = 'identity', 
                         values = NULL, #manual, gradientn, stepsn
                         option = NULL, direction = 1, begin = 0, end = 1, alpha = 1, #viridis
                         low = NULL, mid = NULL, high = NULL, midpoint = NULL, #gradient, gradient2, steps, steps2
                         palette_name = NULL, type = NULL, # brewer
                         bins = NULL, n.breaks = NULL #binned colours
                         ){
  
  scale <- match.arg(scale)
  family <- match.arg(family)
  


  build_scale_obj <- function() {
  if (scale == 'auto'){
    scale <- if(!is.null(values) && is.character(values)) 'discrete' else 'continuous'
  }
  
  if (scale == 'discrete'){
    if (family == 'manual'){
      stopifnot(!is.null(values))
      args <- purrr::compact(list(values = values, name = name, breaks = breaks, labels = labels,
                           limits = limits, guide = guide, na.translate = TRUE, drop = TRUE,
                           na.value = na.value))
      return(do.call(scale_colour_manual, args))
    }
    if (family == 'viridis'){
      args <- purrr::compact(list(option = option, direction = direction, begin = begin, end = end, alpha = alpha,
                           name = name, breaks = breaks, labels = labels, limits = limits, guide = guide,
                           na.value = na.value))
      return(do.call(scale_colour_viridis_d, args))
    }
    if (family == 'brewer'){
      args <- purrr::compact(list(palette = if (is.null(palette_name)) "Set1" else palette_name, type = type, direction = direction,
                           name = name, breaks = breaks, labels = labels, limits = limits, guide = guide,
                           na.value = na.value))
      return(do.call(scale_colour_brewer, args))
    }
    if (family == 'default'){
      agrs <- purrr::compact(list(name = name, breaks = breaks, labels = labels, limits = limits, guide = guide))
    }
    stop('For discrete scales please use family = c(manual, viridis, brewer, default)')
  }
  
  if (scale == 'continuous'){
    if(family == 'viridis'){
      args <- purrr::compact(list(option = option, direction = direction, begin = begin, end = end, alpha = alpha,
                           name = name, breaks = breaks, labels = labels, limits = limits, guide = guide,
                           trans = trans, na.value = na.value))
      return(do.call(scale_colour_viridis_c, args))
    }
    if (family == 'gradient'){
      if (!is.null(values)){
      args <- purrr::compact(list(colours = values, name = name, breaks = breaks, labels = labels,
                           limits = limits, guide = guide, trans = trans, na.value = na.value))
      return(do.call(scale_colour_gradientn, args))
      }
      if (!is.null(mid) || !is.null(midpoint)){
        args <- purrr::compact(list(low = if (is.null(low)) "#2166AC" else low, mid = if (is.null(mid)) 'grey90' else mid, high = if (is.null(high)) "#B2182B" else high,
                             midpoint = if (is.null(midpoint)) 0 else midpoint,
                             name = name, breaks = breaks, labels = labels, limits = limits,
                             guide = guide, trans = trans, na.value = na.value))
        return(do.call(scale_colour_gradient2, args))
      }
      args <- purrr::compact(list(low = if (is.null(low)) "#132B43" else low, high = if (is.null(high)) "#56B1F7" else high,
                           name = name, breaks = breaks, labels = labels, limits = limits,
                           guide = guide, trans = trans, na.value = na.value))
      return(do.call(scale_colour_gradient, args))
    }
  }
  
  if (family == 'brewer'){
    args <- purrr::compact(list(palette = if(is.null(palette_name)) "YlGnBu" else palette_name, type = type, direction = direction,
                         name = name, breaks = breaks, labels = labels, limits = limits,
                         guide = guide, trans = trans, na.value = na.value))
    return(do.call(scale_colour_distiller, args))
  }
  
  if (family == 'steps'){
    if (!is.null(values)) {
      args <- purrr::compact(list(colours = values, name = name, breaks = breaks, labels = labels, limits = limits,
                           guide = guide, trans = trans, na.value = na.value, n.breaks = n.breaks))
      return(do.call(scale_colour_stepsn, args))
    }
    if (!is.null(mid) || !is.null(midpoint)) {
      args <- purrr::compact(list(low = if(is.null(low)) "#2166AC" else low, mid = if (is.null(mid)) "grey90" else mid, high = if (is.null(high)) "#B2182B" else high,
                           midpoint = if (is.null(midpoint)) 0 else midpoint,
                           name = name, breaks = breaks, labels = labels, limits = limits,
                           guide = guide, trans = trans, na.value = na.value, n.breaks = n.breaks))
      return(do.call(scale_colour_steps2, args))
    }
    args <- purrr::compact(list(low = low, high = high,
                         name = name, breaks = breaks, labels = labels, limits = limits,
                         guide = guide, trans = trans, na.value = na.value, n.breaks = n.breaks))
    return(do.call(scale_colour_steps, args))
  }
  
  if (family == "default") {
    args <- purrr::compact(list(name = name, breaks = breaks, labels = labels, limits = limits,
                         guide = guide, trans = trans, na.value = na.value))
    return(do.call(scale_colour_continuous, args))
    stop("For continuous scales, use family = c('viridis', 'brewer', 'gradient', 'steps', 'default'.")
  }
  
  
  if (scale == "binned") {
    if (family == "viridis") {
      args <- purrr::compact(list(option = option, direction = direction, begin = begin, end = end, alpha = alpha,
                           name = name, breaks = breaks, labels = labels, limits = limits, guide = guide,
                           trans = trans, na.value = na.value))
      return(do.call(scale_colour_viridis_b, args))
    }
    if (family == "brewer") {
      # brewer binned -> fermenter
      args <- purrr::compact(list(palette = if (is.null(palette_name)) "YlOrRd" else palette_name, type = type, direction = direction,
                           bins = bins,
                           name = name, breaks = breaks, labels = labels, limits = limits,
                           guide = guide, trans = trans, na.value = na.value))
      return(do.call(scale_colour_fermenter, args))
    }
    if (family == "steps") {
      # steps family already returns binned scales (same mapping as above)
      if (!is.null(values)) {
        args <- purrr::compact(list(colours = values, name = name, breaks = breaks, labels = labels, limits = limits,
                             guide = guide, trans = trans, na.value = na.value, n.breaks = n.breaks))
        return(do.call(scale_colour_stepsn, args))
      }
      if (!is.null(mid) || !is.null(midpoint)) {
        args <- purrr::compact(list(low = if (is.null(low)) "#2166AC" else low, mid = if (is.null(mid)) "grey90" else mid, high = if (is.null(high)) "#B2182B" else high,
                             midpoint = if (is.null(midpoint)) 0 else midpoint,
                             name = name, breaks = breaks, labels = labels, limits = limits,
                             guide = guide, trans = trans, na.value = na.value, n.breaks = n.breaks))
        return(do.call(scale_colour_steps2, args))
      }
      args <- purrr::compact(list(low = low, high = high,
                           name = name, breaks = breaks, labels = labels, limits = limits,
                           guide = guide, trans = trans, na.value = na.value, n.breaks = n.breaks))
      return(do.call(scale_colour_steps, args))
    }
    if (family == "default") {
      args <- purrr::compact(list(name = name, breaks = breaks, labels = labels, limits = limits,
                           guide = guide, trans = trans, na.value = na.value, n.breaks = n.breaks))
      return(do.call(scale_colour_binned, args))
    }
    stop("For binned scales, use family = c('viridis', 'brewer', 'steps', 'default'.")
  }
  
  stop("Unknown combination of scale and family. Please check `scale` and `family`.")
    
  }
  
  scale_obj <- build_scale_obj()
  pass$plot <- pass$plot + scale_obj
  
  pass$plot <- pass$plot + ggnewscale::new_scale_colour()
  
  return(pass)
}
  

gordi_read(m, env, traits)|>
  gordi_species(colour = 'form', linewidth = 2, repel_label = T)|>
  gordi_colour(scale = 'discrete', family = 'manual')
  gordi_label(shortcut = 'upper.upper', shortcut_colour = 'cont')|>
  gordi_colour(scale = 'continuous', family = 'viridis')
  
  

m <- capscale(spe ~ slope, data = env)

gordi_read(m, env, traits)|>
  gordi_species(colour = 'cont', shortcut = 'upper.lower', size = 'cont', shape = 17, stroke = 6)|>
  gordi_colour(what = 'species', scale = 'continuous', family = 'viridis')|>
  gordi_sites(colour = 'elevation', size = 'elevation', alpha = 0.7)|>  
  gordi_colour(what = 'sites', scale = 'continuous', family = 'viridis')
  
gordi_read(m, env, traits)|>
  gordi_species(colour = 'cont', shortcut = 'upper.lower', size = 'cont')|>
  gordi_sites(colour = 'elevation', size = 'elevation', alpha = 0.7)
 
gordi_read(m, env, traits)|>
  gordi_sites(colour = 'elevation', size = 'elevation', alpha = 0.7)|>
  gordi_species(colour = 'cont', shortcut = 'upper.lower', size = 'cont')

gordi_read(m8, env = env)|>
  gordi_sites()|>
  gordi_species()|>
  gordi_colour()

spe_o <- bind_cols(o$spe_name, o$species_scores)|>
  add_row(species = 'Taraxacum sect. Taraxacum', RDA1 = 1, PC1 = 1)|>
  add_row(species = 'Nardus cf. stricta', RDA1 = 1, PC1 = 1)|>
 # add_row(species = 'Plantago atrata subsp. sudetica', RDA1 = 1, PC1 = 1)|>
  left_join(traits|>rename(species = FloraVeg.Taxon), by = 'species')
  # slice_tail(n = 10)

ggplot(data = spe_o, aes(RDA1, PC1))+
  geom_point(aes(colour = cont))+
  gordi_colour(scale = 'binned', family = 'steps', values = c('red', 'blue'))
 
  
ggplot(data = spe_o, aes(RDA1, PC1))+
  geom_point(aes(colour = cont))+
  scale_colour_fermenter() #binned distiller
  scale_colour_distiller() #brewer colour scales, change palette via palette=
  scale_colour_viridis_b() #binned viridis
  scale_colour_viridis_c() #change palette via option=
  scale_colour_stepsn() #define any number of values
  scale_colour_steps2() #possible to define low, mid and high values
  scale_colour_steps(low = 'red', high = 'blue') # cont colour scale, binned, but possible to define low, high values
  scale_colour_binned() #continuous colouring scale, binned (in legend colour steps )
  scale_colour_continuous() # continuous colouring scale, not binned (in legend there is colourbar not colour steps!)
  # gordi_colour(values = c('red', 'green', 'blue', 'black', 'purple'), na.value = NA)+
  # theme_bw()

ggplot(data = spe_o, aes(RDA1, PC1))+
  geom_point(aes(colour = form))+
  scale_colour_viridis_d() #discrete viridis
  scale_colour_brewer(palette = 'Set2') #change palette with palette=
  scale_colour_manual() #defining specific colours with values=
  scale_colour_discrete() #default colour palette, discrete colours

rm(has_sub)


gordi_read(m, env, traits)|>
  gordi_sites(label = F)|>
  gordi_species(label = F, colour = 'cont')|>
  gordi_label(what = 'species', label = 'species_name', shortcut = 'upper.upper', shortcut_colour = 'cont', repel_label = T)|>
  gordi_colour(scale = 'continuous', family = 'viridis')|>
  gordi_label(what = 'sites', label_colour = 'elevation')|>
  gordi_colour(scale = 'continuous', family = 'viridis')

gordi_read(m, scaling = 'species', env, traits) |> 
  gordi_species(label = T, arrow_size = 0.2) |> 
  #gordi_colour(scale = 'continuous', family = 'brewer', palette_name = 'PuBuGn') |> 
  gordi_label(what = 'species', label = 'species_name', shortcut = 'upper.upper',
              repel_label = T, shortcut_length = 4,
              size = 3, shortcut_colour = 'form') |> 
  gordi_colour(scale = 'discrete', family = 'brewer', palette_name = 'Set2', na.value = 'transparent')

gordi_read(m, scaling = 'species', env, traits) |> 
  gordi_species(label = T, colour = 'form', arrow_size = 0.2, alpha = 0.6) |> 
  gordi_colour(scale = 'discrete', family = 'brewer', palette_name = 'Set1') |> 
  gordi_label(what = 'species', label = 'species_name', shortcut = 'upper.upper',
              repel_label = T, shortcut_length = 4,
              size = 3, shortcut_colour = 'form') |> 
  gordi_colour(scale = 'discrete', family = 'brewer', palette_name = 'Set1')

gordi_read(m, scaling = 'species', env, traits) |> 
  gordi_species(label = F) |> 
  #gordi_colour(scale = 'continuous', family = 'brewer', palette_name = 'PuBuGn') |> 
  gordi_label(what = 'species') |> 
  gordi_colour(scale = 'discrete', family = 'brewer', palette_name = 'Set2', na.value = 'transparent')