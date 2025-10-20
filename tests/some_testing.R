remotes::install_github('Kryshtoson/gordi')
library(gordi)
library(tidyverse)
library(vegan)
library(readxl)

devtools::document()


data(dune)
data(dune.env)

read_csv('data/schrankogel/schrankogel_spe.csv') |> 
  select(-logger_ID) -> spe




# gordi_read problem ------------------------------------------------------

auch <- read_excel('C:/Users/User/Documents/R/Gordi_in_community_data_analysis/data/LimestoneQuarries_Aucheno/Aucheno_sum.xlsx')
auch.env <- read_excel('C:/Users/User/Documents/R/Gordi_in_community_data_analysis/data/LimestoneQuarries_Aucheno/Krisi Cesky kras.xls')

pco.bc <- capscale(auch~1, distance = "bray", sqrt.dist = T, binary = F) 
gordi_read(pco.bc, spe = auch, env = auch.env)
#gordi_read(m)



sc <- scores(m, scaling = 'symm', choices = 1:2, correlation = F, hill = F, const = c(2,2), tidy = T) |> as_tibble() |> dplyr::filter(score == 'species') |> dplyr::select(-score)
o <- gordi_read(pco.bc)

o$site_scores$NMDS1 == sc$NMDS1
o$site_scores
o$species_scores

scores(m,
       display = 'all', 
       choices = 1:2, 
       scaling = 'symm', 
       correlation = F, 
       hill = F, 
       const = c(2,2), 
       tidy = T) |> 
  as_tibble() |> 
  filter(str_detect(label, ':|\\*')) |> 
  nrow() > 0




# gordi_corr --------------------------------------------------------------

m <- metaMDS(sqrt(dune))

gordi_read(m, env = dune.env, scaling = 'species') |> 
  gordi_sites() |> 
  gordi_corr(variables = c('A1', 'Use', 'Management'), p_val_adjust = T, perm = 999, label = T, colour = 'variable') |> 
  gordi_colour(scale = 'discrete', family = 'brewer') 


#

scores(m, display = 'sites', scaling = 'sites', const = c(1,2)) ==
  scores(m, display = 'sites', scaling = 'sites', const = c(5,2))


envfit(ord = m,
       env = dune.env[ , c('A1'), drop = FALSE],
       permutations = 999,
       choices = 1:2) |> 
  p.adjust.envfit() -> ef


ef$vectors$r


# Example data for illustration:
patterns <- c("Habitat", "SoilType", "Habitat")
strings <- c("HabitatForest", "SoilTypeClay", "HabitatGrass")

# sub() applies pattern[1] to string[1], pattern[2] to string[2], etc.
result <- sub(patterns, "", strings)
result
# result will be: c("Forest", "Clay", "Grass")






# mess --------------------------------------------------------------------


m <- capscale(sqrt(dune) ~ A1 + Management + Use, dune.env, distance = 'bray', sqrt.dist = T)


# scores(m, display = 'sites', scaling = 'sites', choices = 1:2, correlation = F, hill = T, const = c(1,2))
# scores(m, display = 'species', scaling = 'species', choices = 1:2, correlation = F, hill = T, const = c(1,2))
sco <- scores(m, choices = 1:2, display = 'all', scalling = 'species', tidy = T) |> 
  as_tibble() |> 
  filter(score %in% c('biplot', 'centroids'))

rownames(scores(m, choices = 1:2, display = 'all', scalling = 'species')$centroids)

factor_predictor <- names(m$terminfo$xlev)
vector_predictor <- names(m$terminfo$ordered) |> 
  discard(~ .x %in% names(m$terminfo$xlev))
ordered_factors <- m$terminfo$ordered |> 
  keep(~ .x) |> 
  names()


sco |>
  filter(score == 'centroids') |> 
  filter(grepl(paste(factor_predictor, collapse = '|'), sco[sco$score == "centroids",]$label)) |> 
  mutate(predictor = str_extract(label, paste(factor_predictor, collapse = '|'))) |> 
  mutate(level = str_remove(label, paste(factor_predictor, collapse = '|'))) |> 
  mutate(ordered = grepl(paste(ordered_factors, collapse = '|'), label))

sco |>
  filter(score == "centroids") |>
  # keep only labels containing any of the factor predictors
  filter(str_detect(label, regex(paste(factor_predictor, collapse = "|")))) |>
  mutate(
    predictor = str_extract(label, paste(factor_predictor, collapse = "|")),
    level     = str_remove(label, paste(factor_predictor, collapse = "|")),
    ordered   = str_detect(label, paste(ordered_factors, collapse = "|"))
  )



vectors <- sco |> 
  filter(score == 'biplot') |> 
  filter(grepl(paste(vector_predictor, collapse = '|'), sco[which(sco$score == "biplot"),]$label)) |> 
  bind_cols(vector_predictor) |> 
  rename('predictor' = last_col()) 

names(vectors)[1:2] <- paste0('Axis_pred', 1:2)

grepl(paste(names(m$terminfo$ordered[which(m$terminfo$ordered == T)]), collapse = '|'), names(m$terminfo$ordered))


m |> 
  gordi_read()


sco |> 
  filter(score == 'centroids') |> 
  janitor::clean_names(label)






m <- cca(dune ~ Use + Management, data = dune.env)


m |> 
  gordi_read(choices = 1:2) |> 
  gordi_species() |> 
  gordi_predict(shape = 17, colour = 'predictor') |> 
  gordi_colour(scale = 'discrete', family = 'brewer')

m

m1 <- cca(dune ~ 1, dune.env)

m1 |> 
  gordi_read(env = dune.env) |> 
  gordi_species() |> 
  gordi_corr(variables = c('Use', 'A1'), permutations = 999, p_val_adjust = T,)


plot(m)

summary(m, display = 'all')

m$CCA$centroids |> 
  as_tibble(rownames = 'factor')

devtools::load_all()


# matrix instead spe
data(dune)
data(dune.env)

dune.dist <- vegdist(sqrt(dune), method = 'hell')

m <- capscale(dune.dist ~ A1 + Use, dune.env)
m <- capscale(dune.dist ~ 1)
m <- rda(dune.dist ~ 1) # nope
m <- cca(dune.dist ~ 1) # nope
m <- decorana(dune.dist) # je to blbost, nema se to tak vubec pocitat
m <- decorana(dune)
m1 <- metaMDS(dune.dist, k = 3)
m2 <- metaMDS(dune, k = 3)

# db-RDA, weighted averages (min preferovana moznost)
wascores(scores(m, display = 'lc'), dune)
# db-RDA, 'lc' scores = korelace druhu s LC scores samplu
envfit(m, env = dune, display = 'lc')

# PCoA, 'lc' = linearni korelace, sipky ('lc')
envfit(m, env = dune, display = 'sites')
# PCoA wa zatim nedelame, ale mozna to neni blbost
#...

# NMDS - wa scores = druhova optima ('lc' nemaji smysl)
wa_m1 <- wascores(scores(m1, display = 'si'), dune, expand = T)


# 
# sco.m1 <- wa_m1 |> 
#   as.data.frame() |> 
#   as_tibble() |> 
#   select(NMDS1)
# 
# sco.m2 <- scores(m2, display = 'species') |> 
#   as_tibble() |> 
#   select(NMDS1)
# 
# cor.test(sco.m1 |> pull(), sco.m2 |> pull())



scores(m, tidy = T) |> 
  as_tibble() |> 
  filter(score == 'species') |> 
  summarise(
    # TRUE if there is *any* row where *all* of the first two columns are NA
    has_all_NA_in_first_two = any(if_all(1:2, is.na))
  ) |>
  pull()


spe_scores_and_names <- envfit(m, env = dune, permutations = 0) |> 
  scores(display = 'vectors', choices = 1:2) |> 
  as_tibble(rownames = 'species')

spe_names <- spe_scores_and_names |> select(species)  
spe_scores <- spe_scores_and_names |> select(-species)  


species_scores = tibble::as_tibble(as.data.frame(vegan::scores(m, display = 'species', scaling = scaling, choices = choices, correlation = correlation, hill = hill, const = const)))

envfit(m, env = spe, permutations = 0) |> 
  scores(display = 'vectors', choices = choices) |> 
  as_tibble(rownames = 'species')













# getting species scores when distance matrix provided --------------------



m <- rda(dune ~ 1)
m <- rda(dune ~ A1 + Use, dune.env)
m <- cca(dune ~ 1)
m <- cca(dune ~ A1 + Use, dune.env)
m <- decorana(dune)
m <- capscale(dune ~ 1)
m <- capscale(dune ~ A1 + Use, dune.env)
m <- metaMDS(dune, k = 4)


dbrda_spe_scores_lc <- envfit(m, env = dune, display = 'lc', choices = 1:4, scaling = 'sym', correlation = T, hill = T) |> 
  scores(display = 'vectors') |> 
  as_tibble(rownames = 'label') 


scores(m, display = 'sites', scaling = 'sym', choices = 1:4, correlation = T, hill = F, const = c(1,2), tidy = T) |> as_tibble() 
  
  #filter(score %in% c('biplot', 'centroids'))

# species_scores and species_names
as_tibble(as.data.frame(scores(m, display = 'species', scaling = 'species', choices = 1:4, correlation = T, hill = T, const = c(1,2))), rownames = 'species_names')

# site_scores
as_tibble(as.data.frame(scores(m, display = 'sites', scaling = scaling, choices = choices, correlation = correlation, hill = hill, const = const)))

# axis names
colnames(scores(m, display = 'sites', choices = 1:4))






# gordi predict -----------------------------------------------------------

names(dune.env)


scores(m, tidy = T, correlation = T)

gordi_read(m, env = dune.env, scaling = 'species', correlation = T) |> 
  gordi_species() |> 
  gordi_predict(colour = 'predictor', shape = 'predictor', size = 2, label = T, repel_label = T) 

#



length(names(m$terminfo$ordered) |> 
  discard(~ .x %in% names(m$terminfo$xlev)) ) > 0

obj$pred_df |> 
  print(n = Inf)

gordi_read(m, dune.env) -> obj
obj$predictor_scores |> 
  filter(score == 'centroids') |> 
  nrow() == 0

#

bind_rows(
  tibble(x = c(1,2),
         y = c(2,3)), NULL)


scores(m,
       display = 'all', 
       choices = 1:2, 
       scaling = 'symm', 
       correlation = T, 
       const = c(1,2), 
       tidy = T) |> 
  as_tibble() 
