#' =====================================================================
#' test run input
library(ggrepel)
library(tidyverse)
library(vegan)
library(readxl)

spe <- read_csv('data/schrankogel/schrankogel_spe.csv')[-1] |> 
  log1p()

env <- read_csv("data/schrankogel/schrankogel_env.csv") |>
  mutate(group = elevation > 2500)

trait <- read_xlsx('data/Life_form.xlsx')|>
  select(-SeqID)|>
  pivot_longer(cols = -FloraVeg.Taxon, names_to = 'form', values_to = 'value')


m <- rda(spe ~ 1, distance = "bray") #pcoa
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


inherits(m, 'rda')
class(m)

# cca(spe~ elevation, data = env)
# m <- capscale(spe ~ 1, distance = 'bray')
# m <- capscale(spe ~ elevation + annual_temperature , distance = 'bray', data = env)

m <- decorana(spe)

as_tibble(scores(m, tidy = T))

scores(m)$sites

is.null(m$CCA[1])

class(m)

case_when(
  inherits(m, 'capscale') & !is.null(m$call$distance) & is.null(m$CCA) ~ 'PcoA (capscale)', #via capscale
  inherits(m, 'capscale') & !is.null(m$call$distance) ~ 'db-RDA (capscale)',
  inherits(m, 'rda') & !is.null(m$call$distance) & is.null(m$CCA) ~ 'PCoA (rda)', #if rda and distance but no constrainned variable PCoA
  inherits(m, 'rda') & !is.null(m$call$distance) ~ 'db-RDA (by RDA argument)',
  inherits(m, 'rda') & is.null(m$CCA) ~ 'PCA',
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


gordi_read <- function(m, env, choices = 1:2, scaling = 'symm', correlation = F, hill = F){
  
  #' if(class(m)[[1]]  %in% )
  #'   
  #'   if(is.null(m$CCA)){
  #'     #'... 
  #'   }
  type <- case_when(
    inherits(m, 'capscale') & !is.null(m$call$distance) & is.null(m$CCA) ~ 'PcoA (capscale)', #via capscale
    inherits(m, 'capscale') & !is.null(m$call$distance) ~ 'db-RDA (capscale)',
    inherits(m, 'rda') & !is.null(m$call$distance) & is.null(m$CCA) ~ 'PCoA (rda)', #if rda and distance but no constrainned variable PCoA
    inherits(m, 'rda') & !is.null(m$call$distance) ~ 'db-RDA (rda)',
    inherits(m, 'rda') & is.null(m$CCA) ~ 'PCA',
    inherits(m, 'rda') ~ 'RDA constrained',
    inherits(m, 'cca') & is.null(m$CCA) ~ 'CA',
    inherits(m, 'cca') ~ 'CCA',
    TRUE ~ paste(class(m), collapse = '/') #writes just one output
  )
  
  pass <- list(
    m = m,
    explained_variation = m$CA$eig/m$tot.chi,
    site_scores = as_tibble(as.data.frame(scores(m, scaling = scaling, choices = choices, correlation = correlation, hill = hill)$sites)),
    species_scores = as_tibble(as.data.frame(scores(m, scaling = scaling, choices = choices, correlation = correlation, hill = hill)$species)),
    env = env,
    choices = choices,
    type = type
  )
  
  
  return(pass)
  
}


gordi_read(m, env, choices = 1:2) -> o

m <- rda(spe ~ elevation, distance = 'bray', data = env)
m <- rda(log1p(spe)) #pca or differently rda(spe ~ 1)
m <- rda(sqrt(spe) ~ elevation, data = env) #rda constrained
m <- capscale(spe ~ 1, distance = 'bray') #pcoa
m <- capscale(spe ~ elevation, distance = 'bray', data = env) #db-rda
m <- cca(spe) #ca
m <- cca(spe ~ elevation, data = env) # cca
m7 <- decorana(log1p(spe)) #dca


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





class(o$m)
names(o$site_scores) <- paste0('Axis', 1:2)
bind_cols(env, o$site_scores)


gordi_sites <- function(pass, label = '', colouring = '', repel_label = T) {
  
  #input #' misto pass
  pass <- list(
    m = pass$m,
    explained_variation = pass$explained_variation,
    site_scores = pass$site_scores,
    species_scores = pass$species_scores,
    env = pass$env,
    choices = pass$choices,
    type = pass$type
  )
  
  names(pass$site_scores) <- paste0("Axis", 1:2)
  
  #' here we should control for type of model
  #' capscale
  #' if (class(pass$m)[[1]] == "capscale") {
  #'   if(){} #' dbRDA osetrit
  #'   else{
  #'     actual_labs <- paste0("PCoA", 1:2, " (", round(pass$explained_variation[1:2]*100, 2), '%)')
  #'   }
  #' }
  #' #' cca
  #' else if (class(pass$m)[[1]] == "cca") {
  #  actual_labs <- paste0("CCA", pass$choices, " (", round(pass$explained_variation[1:2]*100, 2), '%)')
  #' }
  #' #' pca
  #' else if (class(pass$m)[[1]] == "pca") {
  #'   actual_labs <- paste0("PCA", 1:2, " (", round(pass$explained_variation[1:2]*100, 2), '%)')
  #' }
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
  
  if (is.null(pass$p)) { #' is.null(input$p) #' mozna lepsi #checks whether p exists in pass, if not it draws plot
    p <- bind_cols(env, pass$site_scores) |>
      ggplot(aes(Axis1, Axis2)) +
      theme_bw() +
      labs(x = actual_labs[1], y = actual_labs[2]) +
      theme(
        text = element_text(size = 15),
        panel.grid = element_blank(),
        legend.justification = c(1, 1)
      )
  } else {
    p #'... 
  }
  
  #' accounting for colouring
  if (colouring == '') {
    p <- p +
      geom_point(size = 3)
  } else {
    p <- p +
      geom_point(aes(colour = !!sym(colouring)), size = 3)
  }
  
  #' accounting for labeling
  if (label != '') {
    if (repel_label) {
      p <- p + geom_text_repel(aes(label = !!sym(label)))
    } else {
      p <- p + geom_text(aes(label = !!sym(label)))
    }
  }
  
  pass$plot <- list(p)
  
  return(pass)
}


obj <- gordi_read(m, env, choices = 1:2) |>
  gordi_sites()



#' gordi_sites(label = "logger_ID", colouring = "group") |>

obj

str(obj)
