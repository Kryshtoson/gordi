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



# 
case_when(
  inherits(m, 'capscale') & !is.null(m$call$distance) & is.null(m$CCA) ~ 'PCoA (capscale)', #via capscale
  inherits(m, 'capscale') & !is.null(m$call$distance) ~ 'db-RDA (capscale)',
  inherits(m, 'rda') & !is.null(m$call$distance) & is.null(m$CCA) ~ 'PCoA (rda)', #if rda and distance but no constrainned variable PCoA
  inherits(m, 'rda') & !is.null(m$call$distance) ~ 'db-RDA (by RDA argument)',
  inherits(m, 'rda') & is.null(m$call$distance) & is.null(m$CCA) ~ 'PCA', inherits(m, 'rda') ~ 'RDA constrained',
  inherits(m, 'cca') & is.null(m$call$distance) & is.null(m$CCA) ~ 'CA', inherits(m, 'cca') ~ 'CCA',
  TRUE ~ paste(class(m), collapse = '/') # writes just one output
)

#' cca {CCA, partial CCA, CA, partial CA}
#' capscale {dbRDA, partial dbRDA, PCoA, partial PCoA}
#' rda {RDA, partial RDA, PCA, partial PCA}
#' ...


#' =====================================================================
#' arguments: 
#' -> model
#' -> headers 


gordi_read <- function(m, env, traits = NULL, choices = 1:2, scaling = 'symm', correlation = F, hill = F){
  
  type <- case_when(
    inherits(m, 'capscale') & !is.null(m$call$distance) & is.null(m$CCA) ~ 'PCoA (capscale)', #via capscale
    inherits(m, 'capscale') & !is.null(m$call$distance) ~ 'db-RDA (capscale)',
    inherits(m, 'rda') & !is.null(m$call$distance) & is.null(m$CCA) ~ 'PCoA (rda)', #if rda and distance but no constrainned variable PCoA
    inherits(m, 'rda') & !is.null(m$call$distance) ~ 'db-RDA (by RDA argument)',
    inherits(m, 'rda') & is.null(m$call$distance) & is.null(m$CCA) ~ 'PCA', inherits(m, 'rda') ~ 'RDA constrained',
    inherits(m, 'cca') & is.null(m$call$distance) & is.null(m$CCA) ~ 'CA', inherits(m, 'cca') ~ 'CCA',
    TRUE ~ paste(class(m), collapse = '/') # writes just one output
  )
  
  pass <- list(
    m = m,
    explained_variation = m$CA$eig/m$tot.chi,
    site_scores = as_tibble(as.data.frame(scores(m, scaling = scaling, choices = choices, correlation = correlation, hill = hill)$sites)),
    species_scores = as_tibble(as.data.frame(scores(m, scaling = scaling, choices = choices, correlation = correlation, hill = hill)$species)),
    env = env,
    traits = traits,
    choices = choices,
    type = type,
    species_names = as_tibble(as.data.frame(scores(m, scaling = scaling, choices = choices, correlation = correlation, hill = hill)$species), rownames = 'species_names')[1]
  )
  
  
  return(pass)
  
}

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

# if else statement
if(o$type == 'CCA'){
  actual_labs <- paste0("CCA", o$choices, " (", round(o$explained_variation[1:2]*100, 2), '%)')
} else if(o$type == 'CA'){
  actual_labs <- paste0("CA", o$choices, " (", round(o$explained_variation[1:2]*100, 2), '%)')
} else if (o$type == 'PCA'){
  actual_labs <- paste0("PCA", o$choices, " (", round(o$explained_variation[1:2]*100, 2), '%)')
} else if(o$type %in% c('PCoA (capscale)', 'PCoA (rda)')){
  actual_labs <- paste0("PCoA", o$choices, " (", round(o$explained_variation[1:2]*100, 2), '%)')
} else if(o$type %in% c('db-RDA (rda)', 'db-RDA (capscale)')){
  actual_labs <- paste0("db_RDA", o$choices, " (", round(o$explained_variation[1:2]*100, 2), '%)')
} else if(o$type == 'RDA constrained'){
  actual_labs <- paste0("RDA", o$choices, " (", round(o$explained_variation[1:2]*100, 2), '%)')
}



#

class(o$m)
names(o$site_scores) <- paste0('Axis', 1:2)
bind_cols(env, o$site_scores)



#
gordi_sites <- function(pass, label = '', colouring = '', size = '', repel_label = T) {
  
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
  
  if (is.null(pass$plot)) { #' is.null(input$p) #' mozna lepsi #checks whether p exists in pass, if not it draws plot
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
  
  
  #' accounting for colouring
  map_colour <- !identical(colouring, '') && has_name(site_df, colouring)
  is_hex <- grepl("^#(?:[A-Fa-f0-9]{6}[A-Fa-f0-9]{3})$", colouring)
  is_named <- colouring %in% grDevices::colours()
  const_colour <- !identical(colouring, '') && !map_colour && (is_hex || is_named)
  
  #' accounting for size
  map_size <- !identical(size, '') && has_name(site_df, size)
  const_size <- !identical(size, '') && is.numeric(size)
  
  # possibility for two colour scales
  p <- p + ggnewscale::new_scale_colour()
  p <- p + ggnewscale::new_scale('size')
  
  # if else for colour and size variables
  if (map_colour && map_size) {
    p <-  p + geom_point(data = site_df,
                 aes(Axis_site1, Axis_site2,
                     colour = !!sym(colouring), 
                     size = !!sym(size)))
  
  } else if (map_colour && const_size) {
    p <- p + geom_point(data = site_df,
                        aes(Axis_site1, Axis_site2, 
                            colour = !!sym(colouring)),
                        size = size)
    
  } else if (const_colour && map_size) {
    p <- p + geom_point(data = site_df,
                        aes(Axis_site1, Axis_site2,
                            size = !!sym(size)),
                        colour = colouring)
 
  } else if (const_colour && const_size) {
    p <- p + geom_point(data = site_df,
                        aes(Axis_site1, Axis_site2),
                        colour = colouring,
                        size = size)
  
  } else if (map_colour) {
    p <- p + geom_point(data = site_df,
                    aes(Axis_site1, Axis_site2,
                        colour = !!sym(colouring)),
                    size = 3)
    
  } else if (const_colour) {
    p <- p + geom_point(data = site_df,
                    aes(Axis_site1, Axis_site2),
                    colour = colouring,
                    size = 3)
    
  } else if (map_size) {
    p <- p + geom_point(data = site_df,
                    aes(Axis_site1, Axis_site2,
                        size = !!sym(size)),
                    colour = 'black')
  
  } else if (const_size) {
    p <- p + geom_point(data = site_df,
                    aes(Axis_site1, Axis_site2),
                    size = size,
                    colour = 'black')
  } else {
    p <- p + geom_point(data = site_df,
                    aes(Axis_site1, Axis_site2), colour = 'black', size = 3)
  }
  
  
  
  
  
  #' accounting for labeling
  if (label != '') {
    if (repel_label) {
      p <- p + geom_text_repel(data = site_df, aes(Axis_site1, Axis_site2, label = !!sym(label)))
    } else {
      p <- p + geom_text(data = site_df, aes(Axis_site1, Axis_site2, label = !!sym(label)))
    }
  }
  
  
  
  
  
  #' pass plot
  pass$plot <- p
  
  return(pass)
}





#

gordi_species <- function(pass, label = '', colouring = '', size = '', shape = '', repel_label = T) {
  
  
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
  

  
  #' accounting for colouring
  map_colour <- !identical(colouring, '') && has_name(spe_df, colouring)
  is_hex <- grepl("^#(?:[A-Fa-f0-9]{6}[A-Fa-f0-9]{3})$", colouring)
  is_named <- colouring %in% grDevices::colours()
  const_colour <- !identical(colouring, '') && !map_colour && (is_hex || is_named)
  
  #' accounting for size
  map_size <- !identical(size, '') && has_name(spe_df, size)
  const_size <- !identical(size, '') && is.numeric(size)
  
  #' accounting for shape
  map_shape <- !identical(shape, '') && has_name(spe_df, shape)
  const_shape <- !identical(shape, '') && is.numeric(shape)
  
  
  # possibility for two colour scales
  p <- p + ggnewscale::new_scale_colour()
  p <- p + ggnewscale::new_scale('size')
  
  
  # if else for colour and size variables
  if (map_colour && map_size && map_shape) {
    p <- p + geom_point(data = spe_df,
                        aes(Axis_spe1, Axis_spe2,
                            colour = !!sym(colouring),
                            size = !!sym(size),
                            shape = !!sym(shape)))
  
  } else if (const_colour && map_size && map_shape) {
    p <- p + geom_point(data = spe_df,
                        aes(Axis_spe1, Axis_spe2,
                            size = !!sym(size),
                            shape = !!sym(shape)),
                        colour = colouring)
    
  } else if (map_colour && const_size && map_shape) {
    p <- p + geom_point(data = spe_df,
                        aes(Axis_spe1, Axis_spe2,
                            colour = !!sym(colouring),
                            shape = !!sym(shape)),
                        size = size)
    
  } else if (map_colour && map_size && const_shape) {
    p <- p + geom_point(data = spe_df,
                        aes(Axis_spe1, Axis_spe2,
                            colour = !!sym(colouring),
                            size = !!sym(size)),
                        shape = shape)
    
  } else if (map_colour && const_size && const_shape) {
    p <- p + geom_point(data = spe_df,
                        aes(Axis_spe1, Axis_spe2,
                            colour = !!sym(colouring)),
                        size = size,
                        shape = shape)
  
  } else if (const_colour && map_size && const_shape) {
    p <- p + geom_point(data = spe_df,
                        aes(Axis_spe1, Axis_spe2,
                            size = !!sym(size)),
                        colour = colouring,
                        shape = shape)
    
  } else if (const_colour && const_size && map_shape) {
    p <- p + geom_point(data = spe_df,
                        aes(Axis_spe1, Axis_spe2,
                            shape = !!sym(shape)),
                        colour = colouring,
                        size = size)
    
  } else if (const_colour && const_size && const_shape) {
    p <- p + geom_point(data = spe_df,
                        aes(Axis_spe1, Axis_spe2),
                        colour = colouring,
                        size = size,
                        shape = shape)
    
  } else if (map_colour && map_size) {
    p <-  p + geom_point(data = spe_df,
                         aes(Axis_spe1, Axis_spe2,
                             colour = !!sym(colouring), 
                             size = !!sym(size)))
    
  } else if (map_colour && const_size) {
    p <- p + geom_point(data = spe_df,
                        aes(Axis_spe1, Axis_spe2, 
                            colour = !!sym(colouring)),
                        size = size)
    
  } else if (const_colour && map_size) {
    p <- p + geom_point(data = spe_df,
                        aes(Axis_spe1, Axis_spe2,
                            size = !!sym(size)),
                        colour = colouring)
    
  } else if (const_colour && const_size) {
    p <- p + geom_point(data = spe_df,
                        aes(Axis_spe1, Axis_spe2),
                        colour = colouring,
                        size = size)
    
  } else if (map_colour && map_shape) {
      p <- p + geom_point(data = spe_df, aes(Axis_spe1, Axis_spe2, colour = !!sym(colouring), shape = !!sym(shape)))
  } else if (map_colour && const_shape) {
      p <- p + geom_point(data = spe_df, aes(Axis_spe1, Axis_spe2, colour = !!sym(colouring)), shape = shape)
  } else if (const_colour && map_shape) {
      p <- p + geom_point(data = spe_df, aes(Axis_spe1, Axis_spe2, shape = !!sym(shape)), colour = colouring)
  } else if (const_colour && const_shape) {
      p <- p + geom_point(data = spe_df, aes(Axis_spe1, Axis_spe2), colour = colouring, shape = shape)
  } else if (map_size && map_shape) {
      p <- p + geom_point(data = spe_df, aes(Axis_spe1, Axis_spe2, size = !!sym(size), shape = !!sym(shape)))
  } else if (map_size && const_shape) {
      p <- p + geom_point(data = spe_df, aes(Axis_spe1, Axis_spe2, size = !!sym(size)), shape = shape)
  } else if (const_size && map_shape) {
      p <- p + geom_point(data = spe_df, aes(Axis_spe1, Axis_spe2, shape = !!sym(shape)), size = size)
  } else if (const_size && const_shape) {
      p <- p + geom_point(data = spe_df, aes(Axis_spe1, Axis_spe2), size = size, shape = shape)
  } else if (map_colour) {
    p <- p + geom_point(data = spe_df,
                        aes(Axis_spe1, Axis_spe2,
                            colour = !!sym(colouring)),
                        size = 3,
                        shape = 16)
    
  } else if (const_colour) {
    p <- p + geom_point(data = spe_df,
                        aes(Axis_spe1, Axis_spe2),
                        colour = colouring,
                        size = 3,
                        shape = 16)
    
  } else if (map_size) {
    p <- p + geom_point(data = spe_df,
                        aes(Axis_spe1, Axis_spe2,
                            size = !!sym(size)),
                        colour = 'black',
                        shape = 16)
    
  } else if (const_size) {
    p <- p + geom_point(data = spe_df,
                        aes(Axis_spe1, Axis_spe2),
                        size = size,
                        colour = 'black',
                        shape = 16)
    
  } else if (map_shape) {
      p <- p + geom_point(data = spe_df, aes(Axis_spe1, Axis_spe2, shape = !!sym(shape)), size = 3)
  } else if (const_shape) {
      p <- p + geom_point(data = spe_df, aes(Axis_spe1, Axis_spe2), shape = shape, size = 3)
  } else {
    p <- p + geom_point(data = spe_df,
                        aes(Axis_spe1, Axis_spe2), size = 3, shape = 16)
  }
  
  
  
  
  #' accounting for labeling
  if (label != '') {
    if (repel_label) {
      p <- p + geom_text_repel(data = spe_df, aes(Axis_spe1, Axis_spe2, label = !!sym(label)))
    } else {
      p <- p + geom_text(data = spe_df, aes(aes(Axis_spe1, Axis_spe2), label = !!sym(label)))
    }
  }
  
  pass$plot <- p
  
  return(pass)
}


# gordi_read(m, env) |> 
#   gordi_species(label = 'species_names')

gordi_read(m, env, trait) |> 
  gordi_sites(size = 'elevation') |> 
  gordi_species(colouring = 'cont')

gordi_read(m, env, trait) |> 
  gordi_species(colouring = 'form', size = 'cont', shape = 17)
 
gordi_read(m, env, trait) |> 
  gordi_sites(colouring = 'group', size = 'elevation') |> 
  gordi_species(label = 'species_names', colouring = 'form', size = 'cont')





