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
  pivot_longer(cols = -FloraVeg.Taxon, names_to = 'form', values_to = 'value')




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


gordi_read <- function(m, env, choices = 1:2, scaling = 'symm', correlation = F, hill = F){
  
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
    choices = choices,
    type = type,
    species_names = as_tibble(as.data.frame(scores(m, scaling = scaling, choices = choices, correlation = correlation, hill = hill)$species), rownames = 'species_names')[1]
  )
  
  
  return(pass)
  
}

#
#as_tibble(as.data.frame(scores(m, scaling = 'symm', choices = 1:2, correlation = T, hill = T)$sites), rownames = 'species_names')[1]

gordi_read(m, env, choices = 1:2) -> o
o$type
o$plot

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
gordi_sites <- function(pass, label = '', colouring = '', repel_label = T) {
  
  # #input #' misto pass
  # pass <- list(
  #   m = pass$m,
  #   explained_variation = pass$explained_variation,
  #   site_scores = pass$site_scores,
  #   species_scores = pass$species_scores,
  #   env = pass$env,
  #   choices = pass$choices,
  #   type = pass$type,
  #   species_names = pass$species_names,
  #   plot = pass$plot
  # )
  
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
  #  p #'... 
  }
  
  site_df <- bind_cols(pass$env, pass$site_scores)
  
  # colour definition
  map_colour <- !identical(colouring, '') && has_name(site_df, colouring)
  is_hex <- grepl("^#(?:[A-Fa-f0-9]{6}[A-Fa-f0-9]{3})$", colouring)
  is_named <- colouring %in% grDevices::colours()
  const_colour <- !identical(colouring, '') && !map_colour && (is_hex || is_named)
  
  
  #' accounting for colouring
  if (map_colour) {
    p <- p +
      geom_point(data = site_df, aes(Axis_site1, Axis_site2, colour = !!sym(colouring)), size = 3)
  } else if (const_colour) {
    p <- p + geom_point(data = site_df, aes(Axis_site1, Axis_site2), size = 3, colour = colouring)  
  } else {
    p <- p +
      geom_point(data = site_df, aes(Axis_site1, Axis_site2), size = 3)
  }
  
  #' accounting for labeling
  if (label != '') {
    if (repel_label) {
      p <- p + geom_text_repel(data = site_df, aes(Axis_site1, Axis_site2, label = !!sym(label)))
    } else {
      p <- p + geom_text(data = site_df, aes(Axis_site1, Axis_site2, label = !!sym(label)))
    }
  }
  
  pass$plot <- p
  
  return(pass)
}




 o <- gordi_read(m, env) |> 
  gordi_sites(colouring = 'group')

o <- gordi_read(m, env) |> 
  gordi_species()

o <- gordi_read(m, env) |> 
  gordi_sites() |> 
  gordi_species()

o$plot

#

gordi_species <- function(pass, label = '', colouring = '', repel_label = T) {
  
  # #input #' misto pass
  # pass <- list(
  #   m = pass$m,
  #   explained_variation = pass$explained_variation,
  #   site_scores = pass$site_scores,
  #   species_scores = pass$species_scores,
  #   env = pass$env,
  #   choices = pass$choices,
  #   type = pass$type,
  #   species_names = pass$species_names,
  #   plot = pass$plot
  # )
  
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
  
  spe_df <- bind_cols(pass$species_scores, pass$species_names)

  #' accounting for colouring
  if (colouring == '') {
    p <- p +
      geom_point(data = spe_df, aes(Axis_spe1, Axis_spe2), size = 3)
  } else {
    p <- p +
      geom_point(data = spe_df, aes(Axis_spe1, Axis_spe2, colour = !!sym(colouring)), size = 3)
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

o <- gordi_read(m, env) |> 
  gordi_sites() |> 
  gordi_species()
