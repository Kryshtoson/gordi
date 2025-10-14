remotes::install_github('Kryshtoson/gordi')
library(gordi)
library(tidyverse)
library(vegan)
library(rlang)

devtools::document()

data(dune)
data(dune.env)

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


