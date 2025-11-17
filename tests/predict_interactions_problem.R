library(tidyverse)
library(vegan)
library(readxl)
library(gordi)
library(patchwork)
library(ggpubr)


devtools::document()

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

dbrda1 <- capscale(spe ~ slope + carbon + alt_class + ph_class +
                   slope:carbon + carbon:alt_class + alt_class:ph_class,
                   sqrt.dist = T, distance = 'bray',
                   data = env)

dbrda1 <- capscale(spe ~ slope + elevation + slope:elevation + ph_class:alt_class + Condition(ph_class + alt_class),
                   sqrt.dist = T, distance = 'bray',
                   data = env)

cca1 <- cca(spe ~ annual_temperature + carbon + alt_class + ph_class +
                   annual_temperature:carbon + carbon:alt_class + alt_class:ph_class,
                   sqrt.dist = T, distance = 'bray',
                   data = env)

rda1 <- rda(spe ~ annual_temperature + carbon + alt_class + ph_class +
                   annual_temperature:carbon + carbon:alt_class + alt_class:ph_class,
                   sqrt.dist = T, distance = 'bray',
                   data = env)



# gordi_predict -----------------------------------------------------------

gordi_read(dbrda1, scaling = 'spe', correlation = T, env = env, spe = spe, const = c(3,1)) |> 
  #gordi_fit(slice_max = 30) |> 
  #gordi_species() |> 
  #gordi_sites(size = 1) |> 
  gordi_predict(show_label = T, label = 'variable_level', colour = 'class', repel_label = T) 
  #gordi_label(what = 'species', shortcut = 'upper.lower', size = 2, shortcut_colour = 4) 


gordi_read(cca1, scaling = 'spe', correlation = T, env = env, spe = spe) |> 
  gordi_species(size = 1) |> 
  gordi_sites(size = 1) |> 
  gordi_predict2(show_label = T) 

gordi_read(rda1, scaling = 'spe', correlation = F, env = env, spe = spe) |> 
  gordi_species(size = 1) |> 
  #gordi_sites(size = 1) |> 
  gordi_predict2(show_label = T) 



o$predictor_scores |> 
  select(level, variable_level)

#devtools::document()

o$inter_ef
o$vector_inter_scores
o$interaction_df_fct
o$inter_terms

o$vector_inter_scores |> 
  mutate(variable = stringr::str_extract(
    variable_level,
    pattern = paste(o$inter_terms, collapse = '|')
    )) |> 
  mutate(level = stringr::str_replace(
    variable_level, 
    pattern = variable,
    replacement = ''
    )) |> 
  mutate(score = 'biplot') |> 
  separate_wider_delim(variable, delim = ':', names = c('var1', 'var2')) |>
  mutate(lvl2 = str_remove(level, '_')) |> 
  mutate(variable = paste0(var1, ":", var2),
         level = paste0(var1, ":", var2, '-', lvl2),
         variable_level = paste0(var1, ":", var2, '-', lvl2)
         ) |> 
  select(-c(var1, var2, lvl2)) |> 
  relocate(c(score, variable, level, variable_level), .after = 2)
  


o$factor_inter_scores |> 
  rename(variable = interacting_variables,
         level = level_combinations) |> 
  mutate(score = 'centroids') |> 
  separate_wider_delim(variable, delim = ':', names = c('var1', 'var2')) |> 
  separate_wider_delim(level, delim = ':', names = c('level1', 'level2')) |> 
  mutate(variable = paste0(var1, ":", var2),
         level = paste0(level1, ":", level2),
         variable_level = paste0(var1, '-', level1, ":", var2, "-", level2)
         ) |> 
  select(-c(var1, var2, level1, level2))




o$main_terms

o$factor_scores |> 
  rename(variable_level = label) |> 
  mutate(variable = stringr::str_extract(
    variable_level,
    pattern = paste(o$main_terms, collapse = '|')
    )) |> 
  mutate(level = stringr::str_replace(
    variable_level, 
    pattern = variable,
    replacement = ''
  )) |> 
  mutate(variable_level = paste0(variable, '-', level))
  relocate(c(variable, level, variable_level), .after = score)
  
  
o$vector_scores |> 
  rename(variable_level = label) |> 
  mutate(variable = stringr::str_extract(
    variable_level,
    pattern = paste(o$main_terms, collapse = '|')
    )) |> 
  mutate(level = stringr::str_replace(
    variable_level, 
    pattern = variable,
    replacement = NA_character_
  )) |> 
  mutate(variable_level = variable) |> 
  relocate(c(variable, level, variable_level), .after = score)


scores(dbrda1, scaling = 'sym', const = c(2,2), correlation = F, tidy = T) |> 
  as_tibble() |> 
  filter(score == 'biplot')

ef <- envfit(dbrda1,
       env = interaction_df_vct,
       display = 'lc',
       scaling = 'non',
       choices = 1:2,
       correlation = F,
       const = c(2,2),
       hill = T) 
scores(ef, display = 'bp')


plot(dbrda1, display = c('bp', 'sp', 'si'), scaling = 'spe')
text(dbrda1, display = 'bp', scaling = 'spe')
plot(ef)

# scaling species - works well in envfit, not necessary to recalculate it

# scaling sites
spe_scores <- scores(ef, display = 'bp') 
constant <- 2

spe_scores[,1] * constant
spe_scores[,2] * constant

# scaling symmetric
spe_scores <- scores(ef, display = 'bp') 
constant <- 2
lambda <- eigenvals(dbrda1)
sum_lambda <- sum(lambda)

spe_scores[,1] * (lambda[1]/sum_lambda)^(1/4) * constant
spe_scores[,2] * (lambda[2]/sum_lambda)^(1/4) * constant



### --- FOR LOOP --- ####

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


interaction_df_vct <- NULL # here, interaction_df is created as an empty tibble, with the correct number of rows, but no columns
interaction_df_fct <- NULL

  # if (length(inter_terms) > 0) {           

    for (i in pull(distinct(interaction_table, inter_ID))) {

      # prepare a vector with all interaction terms
      inter <- interaction_table |>
        filter(inter_ID == i) |>
        select(variable) |>
        pull()

      # find out which class is which variable
      inter_class <- env |>
        select(all_of(inter)) |>
        map_chr(class)


      # object to hold the result of the current iteration
      current_df_vct <- NULL
      current_df_fct <- NULL 

      # --- NUMERIC x NUMERIC ---
      if (all(inter_class %in% c('numeric', 'integer', 'double'))) {
        current_df_vct <- as_tibble(env[,inter[1]] * env[,inter[2]]) # multiply interacting predictors
        colnames(current_df_vct) <- paste(inter[1], inter[2], sep = ':')       #
      }
  
  
      # --- NUMERIC x FACTOR/CHARACTER ---
      if (any(inter_class %in% c('character', 'factor')) &&
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

        current_df_vct <- as_tibble(as.vector(inter_df_vct) * as.data.frame(inter_df_fct))
        names(current_df_vct) <- paste(names(inter_df_vct), names(inter_df_fct), sep = ":")
      }
      
      
      if (is.null(interaction_df_vct)) {
        interaction_df_vct <- current_df_vct
      } else {
        interaction_df_vct <- bind_cols(interaction_df_vct, current_df_vct)
      }
      
      
      # --- FACTOR/CHARACTER x FACTOR/CHARACTER ---
      # in this case, it will calculate directly the centroids of LC site scores
      # for each level combination

      if (all(inter_class %in% c('character', 'factor'))) {
        var1_name <- inter[1]
        var2_name <- inter[2]

        final_col_name <- paste(var1_name, var2_name, sep = ":")

        current_df_fct <- scores(dbrda1,
                                 choices = 1:2,
                                 scaling = 'non',
                                 correlation = F,
                                 hill = F,
                                 const = c(2,2),
                                 tidy = T) |>
          as_tibble() |>
          filter(score == 'constraints') |>
          bind_cols(env) |>
          group_by(pick(c(var1_name, var2_name))) |>
          summarise(CAP1 = mean(CAP1, na.rm = T),
                    CAP2 = mean(CAP2, na.rm = T)) |>
          ungroup() |>
          unite({{final_col_name}}, where(is.character), sep = ':') |>
          relocate(c(CAP1, CAP2), .before = 1) |>
          pivot_longer(-where(is.numeric), names_to = 'interacting_variables', values_to = 'level_combinations')
      }


      if (is.null(interaction_df_fct)) {
        interaction_df_fct <- current_df_fct
      } else {
        interaction_df_fct <- bind_cols(interaction_df_fct, current_df_fct)
      }



  } # end of if() that starts before forloop
   
  interaction_df_vct
  interaction_df_fct

#



























# ordiplot vs gordi ----------------------------------------------------------------

plot(dbrda1, scaling = 'non', correlation = T, xlim = c(0, 0))

gordi_read(dbrda1, scaling = 'non', correlation = T, env = env, spe = spe) |> 
  gordi_species(symbol = 'point', size = 1, shape = 3, colour = 2) |> 
  gordi_sites(size = 1, colour = 'black') |> 
  gordi_predict2(colour = 'green', linewidth = 0.2, shape = 15, size = 2)






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

for (i in pull(distinct(interaction_table[-c(5,6),], inter_ID))) {
  
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
  } 
  
  # --- NUMERIC x FACTOR/CHARACTER ---  
  if (any(inter_class %in% c('character', 'factor')) && 
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
  } 
  
  # # --- FACTOR/CHARACTER x FACTOR/CHARACTER ---
  # if (all(inter_class %in% c('character', 'factor'))) {
  #   var1_name <- inter[1]  
  #   var2_name <- inter[2]  
  #   
  #   final_col_name <- paste(var1_name, var2_name, sep = ":")
  #   
  #   current_df <- scores(dbrda1, scaling = 'spe', correlation = T, tidy = T) |> 
  #     as_tibble() |> 
  #     filter(score == 'constraints') |> 
  #     bind_cols(env) |> 
  #     group_by(pick(c(var1_name, var2_name))) |> 
  #     summarise(CAP1 = mean(CAP1, na.rm = T),
  #               CAP2 = mean(CAP2, na.rm = T)) |> 
  #     ungroup() |> 
  #     unite({{final_col_name}}, where(is.character), sep = ':') |> 
  #     relocate(c(CAP1, CAP2), .before = 1)
  #   
  # }
  
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
