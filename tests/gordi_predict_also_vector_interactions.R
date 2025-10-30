library(tidyverse)
library(vegan)
library(readxl)
library(gordi)
library(patchwork)
library(ggpubr)



# Import ------------------------------------------------------------------

spe <- read_csv('data/schrankogel/schrankogel_spe.csv') |> 
  select(-logger_ID) |> 
  sqrt()

env <- read_csv('data/schrankogel/schrankogel_env.csv') |> 
  mutate(alt_class = case_when(elevation >= 2600 ~ 'high',
                               elevation < 2600 ~ 'low',
                               TRUE ~ NA),
         ph_class = case_when(pH < 4.9 ~ 'acid',
                              pH < 5.3 ~ 'mesic',
                              TRUE ~ 'basic'))


glimpse(env)

# RDA ---------------------------------------------------------------------

rda_1 <- rda(spe ~ annual_temperature + carbon + alt_class + ph_class +
             annual_temperature:carbon + carbon:alt_class + alt_class:ph_class,
             data = env)

rda_1

# predictors
table(env$alt_class)            # factor: 2 levels: low and high
table(env$ph_class)             # factor: 3 levels: acid, mesic, basic
summary(env$annual_temperature) # vector
summary(env$carbon)             # vector

# interaction of two factors
table(env$alt_class, env$ph_class) # factor: 6 levels:
                                   # acid-high, acid-low,
                                   # mesic-high, mesic-low,
                                   # basic-low, basic-high
# interaction of factor and vector 
env |> 
  mutate(values = 1) |> 
  pivot_wider(names_from = alt_class, values_from = values, values_fill = 0) |> 
  select(high, low) -> alt_class_dummy
summary(env$carbon * alt_class_dummy) # 2 vectors based on the 2 factor levels

# interaction of two vectors
env$annual_temperature * env$carbon # 1 vector created by multiplication
                                    # of the two vectors interacting


# Scores ------------------------------------------------------------------

sco_rda_1 <- scores(rda_1, scaling = 'sites', choices = 1:2, tidy = T) |> 
  as_tibble()

table(sco_rda_1$score)
### for predictors:
# 5 centroids - I have all
# 2 biplots - I have all
### for interactions:
# 6 centroids - I have only half of them
# 3 biplots - 

sco_rda_1 |> 
  filter(!score %in% c('species', 'sites', 'constraints')) |> 
  print(n = Inf)



# Get predictor scores ----------------------------------------------------
formula(rda_1)
terms(rda_1)

attr(terms(rda_1), "term.labels") |> str_detect(':')
attr(terms(rda_1), "variables")
attr(terms(rda_1), "factors")

sco_rda_1 |> 
  filter(score %in% c('biplot', 'centroids') & !str_detect(label, ':'))

term_labels <- attr(terms(rda_1), "term.labels")

term_labels[str_detect(term_labels, ':')]


# RDA
rda_1 <- rda(spe ~ annual_temperature + carbon + alt_class + ph_class +
               annual_temperature:carbon + carbon:alt_class + alt_class:ph_class,
             data = env)
rda_1

rda_2 <- rda(spe ~ annual_temperature + carbon + alt_class + ph_class,
             data = env)
rda_2

rda_3 <- rda(spe ~ annual_temperature + carbon,
             data = env)
rda_3

rda_4 <- rda(spe ~ alt_class + ph_class,
             data = env)
rda_4

rda_5 <- rda(spe ~ alt_class:ph_class + Condition(alt_class + ph_class),
             data = env)
rda_5

# CCA
cca_1 <- cca(spe ~ annual_temperature + carbon + alt_class + ph_class +
               annual_temperature:carbon + carbon:alt_class + alt_class:ph_class,
             data = env)
cca_1

cca_2 <- cca(spe ~ annual_temperature + carbon + alt_class + ph_class,
             data = env)
cca_2

cca_3 <- cca(spe ~ annual_temperature + carbon,
             data = env)
cca_3

cca_4 <- cca(spe ~ alt_class + ph_class,
             data = env)
cca_4


# db-RDA
dbrda_1 <- capscale(spe ~ annual_temperature + carbon + alt_class + ph_class +
               annual_temperature:carbon + carbon:alt_class + alt_class:ph_class,
             data = env)
dbrda_1

dbrda_2 <- capscale(spe ~ annual_temperature + carbon + alt_class + ph_class,
                    data = env)
dbrda_2

dbrda_3 <- capscale(spe ~ annual_temperature + carbon,
                    data = env)
dbrda_3

dbrda_4 <- capscale(spe ~ alt_class + ph_class,
                    data = env)
dbrda_4



# gordi
gordi_read(rda_5, env = env |> 
             bind_cols(interaction_df), spe = spe, scaling = 'si', correlation = T) |>
  # gordi_fit(abs_frequency = 30) |> 
  gordi_species(alpha = 0) |> 
  gordi_predict(show_label = T, repel_label = T, colour = 'score') |> 
  gordi_corr(variables = c('alt_class:ph_class'), show_label = T)

data(dune)
data(dune.env)

devtools::document()

# --- 1. Example with main effects ---
m1 <- capscale(dune ~ A1 + Management, data = dune.env)
gordi_read(m1, env = dune.env, scaling = 'species', correlation = T) |>
   gordi_species(label = F) |>
   # Predictors (A1: continuous arrow, Management: categorical centroids)
   # Colour is dynamically mapped to the 'score' column (biplot or centroid)
   gordi_predict2(scaling_coefficient = 1, colour = 'score', size = 4)

# --- 2. Example with an interaction term ---
# The interaction term will be calculated post-hoc via envfit
m2 <- capscale(dune ~ Use * Management, data = dune.env)
gordi_read(m2, env = dune.env, scaling = 'species', correlation = T) |>
   gordi_sites() |>
   # Predictors include main effects and interaction effects (e.g., A1:Management)
   gordi_predict2(show_label = T, repel_label = T, colour = 'score', size = 3)


dune.env$Use <- as.character(dune.env$Use)



interaction_terms <- attr(terms(cca_1), 'term.labels') |> 
  stringr::str_subset(':')

#interaction_terms <- c('a:b', 'a:b:c')

interaction_table <- interaction_terms |>
  str_split(pattern = ":") |>
  set_names(~ paste0("inter_", seq_along(.))) |> 
  imap_dfr(~ tibble(
    inter_ID = .y,
    inter_var = paste0("var", seq_along(.x)),
    variable = .x
  )) |> 
  rowwise() |> 
  mutate(var_class = case_when(
    !(variable %in% names(env)) ~ 'Not found',
    is.numeric(env[[variable]]) | is.integer(env[[variable]]) ~ 'vector',
    is.factor(env[[variable]]) | is.character(env[[variable]]) ~ 'factor',
    TRUE ~ class(env[[variable]])[1])
  ) |> 
  ungroup() 

# interaction_table <- interaction_terms |> 
#   set_names(paste0("inter_", 1:length(interaction_terms))) |> # Give unique IDs
#   map(as_tibble_col, column_name = "variable") |>
#   list_rbind(names_to = "inter_ID")


interaction_df <- tibble(.rows = nrow(env))


### for loop  

for (i in pull(distinct(interaction_table, inter_ID))) {
  
  # filter out interaction terms
  inter <- interaction_table |> 
    filter(inter_ID == i) |> 
    select(variable) |> 
    pull()
  
  inter_class <- env |>
    select(all_of(inter)) |> 
    map_chr(class) 
  
  # Variable to hold the result of the current iteration
  current_df <- NULL
  

# --- NUMERIC x NUMERIC ---
  if (all(inter_class %in% c('numeric', 'integer', 'double'))) {
    current_df <- as_tibble(env[,inter[1]] * env[,inter[2]])
    colnames(current_df) <- paste(inter[1], inter[2], sep = ':')

# --- NUMERIC x FACTOR/CHARACTER ---        
  } else if (any(inter_class %in% c('character', 'factor')) && 
             any(inter_class %in% c('numeric', 'integer', 'double'))) {
    inter_df_vct <- NULL
    inter_df_fct <- NULL
    
    if (inter_class[1] %in% c('character', 'factor')) {
      inter_df_fct <- fastDummies::dummy_cols(env[,inter[1]]) |> select(-1)
    } else {
      inter_df_vct <- env[,inter[1]]
    }
    
    if (inter_class[2] %in% c('character', 'factor')) {
      inter_df_fct <- fastDummies::dummy_cols(env[,inter[2]]) |> select(-1)
    } else {
      inter_df_vct <- env[,inter[2]]
    }
    
    current_df <- as_tibble(as.vector(inter_df_vct) * as.data.frame(inter_df_fct)) 
    names(current_df) <- paste(names(inter_df_vct), names(inter_df_fct), sep = ":")
    
# --- FACTOR/CHARACTER x FACTOR/CHARACTER ---    
  } else if (all(inter_class %in% c('character', 'factor'))) {
    var1_name <- inter[1]  
    var2_name <- inter[2]  
    
    final_col_name <- paste(var1_name, var2_name, sep = ":")
    
    current_df <- env |>
      mutate(
        interaction_term = paste(
          paste0(var1_name, "_", .data[[var1_name]]),
          paste0(var2_name, "_", .data[[var2_name]]),
          sep = ":")) |>
      select(interaction_term) |> 
      setNames(final_col_name)
    
  }

  if (ncol(interaction_df) == 0) {
    interaction_df <- current_df
  } else {
    interaction_df <- bind_cols(interaction_df, current_df)
  }

}

interaction_df

# Calcularing enfit -------------------------------------------------------


ef <- envfit(rda_1, env = env |> 
               select(annual_temperature, carbon, alt_class, ph_class) |> 
               bind_cols(interaction_df), display = 'lc', scaling = 'si', correlation = T)

# vector x vector, vector x factor - ok is only species, lc
bind_rows(as_tibble(scores(ef, display = 'bp'), rownames = 'label'), 
          as_tibble(scores(ef, display = 'cn'), rownames = 'label'))


# scores from ordination
scores(rda_1, scaling = 'si', correlation = T, display = 'bp')

# sco_rda_1 |> 
#  filter(score == 'biplot')



# factor x factor - ani jedna varianta se nepotkava
ef <- envfit(rda_1, env = interaction_df, display = 'lc', scaling = 'si')
as_tibble(scores(ef, display = 'factors'), rownames = 'label')

# scores from ordination
sco_rda_1 <- scores(rda_1, scaling = 'si', tidy = T) |> 
  as_tibble()
sco_rda_1 |> 
  filter(score == 'biplot')


# factor x factor
sco_rda_1 |> 
  filter(score == 'biplot')



sco_rda_1 |> 
  filter(score == 'centroids')

sco_rda_1 |> 
  filter(score == "factorbiplot")


sco_rda_1 |> 
  filter(score == 'constraints') |> 
  bind_cols(env |> 
              select(alt_class)) |> 
  group_by(alt_class) |> 
  summarise(RDA1 = mean(RDA1),
            RDA2 = mean(RDA2))

sco_rda_1 |> 
  filter(score == 'constraints') |> 
  bind_cols(env |> 
              select(ph_class)) |> 
  group_by(ph_class) |> 
  summarise(RDA1 = mean(RDA1),
            RDA2 = mean(RDA2))



ordiplot(rda_1, scaling = 'spe')
plot(ef, col = 3)




