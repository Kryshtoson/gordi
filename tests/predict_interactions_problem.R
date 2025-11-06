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


# ordination --------------------------------------------------------------

dbrda1 <- capscale(spe ~ annual_temperature + carbon + alt_class + ph_class +
                   annual_temperature:carbon + carbon:alt_class + alt_class:ph_class,
                   sqrt.dist = T, distance = 'bray',
                   data = env)



# gordi_predict -----------------------------------------------------------

gordi_read(dbrda1, scaling = 'non', correlation = T, env = env, spe = spe) |> 
  gordi_species(alpha = 0.2) |> 
  gordi_predict(colour = 'score', show_label = T)


# scores.rda --------------------------------------------------------------

sco <- scores(dbrda1, scaling = 'spe', correlation = T, tidy = T) |> 
  as_tibble() |> 
  filter(score %in% c('biplot', 'centroid')) 

sco |> 
  print(n = Inf)

ggplot() +
  geom_point(data = sco |> filter(score == 'centroids'),
             aes(x = CAP1, y = CAP2)) +
  





# interaction df ----------------------------------------------------------

interaction_terms <- attr(terms(dbrda1), 'term.labels') |> 
  stringr::str_subset(':')

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



# scores.envfit -----------------------------------------------------------

ef <- envfit(dbrda1, env = interaction_df, choices = 1:2, correlation = T, scaling = 'spe', const = c(2,2), display = 'lc')

sco.ef <- bind_rows(
  scores(ef, display = 'vectors') |> as_tibble(rownames = 'label'),
  scores(ef, display = 'factors') |> as_tibble(rownames = 'label'))

ef$vectors$arrows

sco <- scores(dbrda1, scaling = 'spe', choices = 1:2, correlation = T, const = c(2,2), tidy = T) |> 
  as_tibble() |> 
  filter(score %in% c('biplot', 'centroid')) 

sco |> 
  print(n = Inf)


sco.ef |> 
  print(n = Inf)

o <- gordi_read(dbrda1, scaling = 'spe', choices = 1:2, correlation = T, const = c(2,2), env = env, spe = spe) |> 
  gordi_species(alpha = 0.2) |> 
  gordi_predict(colour = 'score')

o$predictor_scores

aa <- ordiplot(dbrda1, scaling = 'non')
text(aa, what = 'bip')

aa$species
aa$sites
aa$biplot
