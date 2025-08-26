library(vegan)
library(tidyverse)
library(ggrepel)
library(readxl)
library(ggnewscale)

# Import -----------------------------------------------------------------------
env <- read_csv("data/schrankogel/schrankogel_env.csv") |>
  mutate(group = case_when(elevation > 3000 ~ '1',
                           elevation <= 3000 & elevation > 2500 ~ '2',
                           elevation <= 2500 ~ '3',
                           TRUE ~ 'wtf'),
         logger_ID = as.character(logger_ID)) 

spe <- read_csv('data/schrankogel/schrankogel_spe.csv') |> 
  mutate(logger_ID = as.character(logger_ID)) |> 
  select(where(~ sum(. != 0) > 5)) |> 
  filter(if_any(where(is.numeric), ~ . != 0)) 

env <- env |> 
  semi_join(spe, by = 'logger_ID') 

spe <- spe |> 
  select(!logger_ID) |> 
  log1p()


trait <- read_xlsx('data/Life_form.xlsx')|>
  select(-SeqID)|>
  pivot_longer(cols = -FloraVeg.Taxon, names_to = 'form', values_to = 'value') |> 
  filter(!value == 0) |> 
  distinct(FloraVeg.Taxon, .keep_all = T) |> 
  mutate(cont = rep(1:50, length.out = n()))


# PCA
m <- rda(spe ~ elevation + slope, data = env)

gordi_read(m, scaling = 'species', env, trait) |> 
  gordi_species(label = T, colour = 'form', arrow_size = 0.2, alpha = 0.6) |> 
  gordi_colour(scale = 'discrete', family = 'brewer', palette_name = 'Set1') |> 
  gordi_label(what = 'species', label = 'species_name', shortcut = 'upper.upper',
              repel_label = T, shortcut_length = 4,
              size = 3, shortcut_colour = 'form', max.overlaps = 100) |> 
  gordi_colour(scale = 'discrete', family = 'brewer', palette_name = 'Set1') |> 
  gordi_predict(label = 'predictor_names', colour = 'predictor_names') |> 
  gordi_colour(scale = 'discrete', family = 'manual', values = c('pink', 'purple')) |> 
  gordi_label(what = 'predictor', label_colour = 'predictor_names', nudge_x = 0.18, nudge_y = 0.05) |> 
  gordi_colour(scale = 'discrete', family = 'manual', values = c('pink', 'purple'))


# CCA
m <- cca(spe, data = env)

gordi_read(m, env, trait) |> 
  gordi_species(label = F, colour = 'form', fill = 'cont', shape = 22, stroke = 1.5) |> 
  gordi_colour(scale = 'discrete', family = 'viridis', direction = -1) |> 
  gordi_colour(scale = 'continuous', family = 'brewer', palette_name = 'Set1', fill = T, direction = -1)


data(dune)
data(dune.env)
mod <- capscale(dune ~ 1)
res <- gordi_read(mod) |> 
  gordi_species()
res


m <- capscale(dune ~ A1, data = dune.env)
gordi_read(m, env = dune.env, scaling = 'species', correlation = T) |> 
  gordi_species(label = F) |> 
  gordi_predict(scaling_coefficient = 0.1, label = F)
