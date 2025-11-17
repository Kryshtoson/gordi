library(vegan)
library(readxl)
library(tidyverse)


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

# db-RDA ordination --------------------------------------------------------------

dbrda1 <- capscale(spe ~ annual_temperature + carbon + alt_class + ph_class +
                   annual_temperature:carbon + carbon:alt_class + alt_class:ph_class,
                   sqrt.dist = T, distance = 'bray',
                   data = env)


# a bit edited FOR LOOP from gordi_predict  -------------------------------

# extract interaction terms from the model
interaction_terms <- attr(terms(dbrda1), 'term.labels') |> 
  stringr::str_subset(':')

# create a table of what interacts with what
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


# prepare empty tibbles for for loop
interaction_df_vct <- NULL  # <- THIS MUST BE NULLED BEFORE EVERY USE OF the FORLOOP
interaction_df_fct <- NULL  # <- THIS MUST BE NULLED BEFORE EVERY USE OF the FORLOOP

## for loop
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
     current_df_vct <- as_tibble(env[,inter[1]] * env[,inter[2]]) 
     colnames(current_df_vct) <- paste(inter[1], inter[2], sep = ':')       
     }
   # simply multiply interacting predictors
  
   
   # --- NUMERIC x FACTOR/CHARACTER ---
   if (any(inter_class %in% c('character', 'factor')) && any(inter_class %in% c('numeric', 'integer', 'double'))) {

   inter_df_vct <- NULL
   inter_df_fct <- NULL
   
        # if the first variable is character or factor, create dummy cols, else just take it as it is
        if (inter_class[1] %in% c('character', 'factor')) {
          inter_df_fct <- fastDummies::dummy_cols(env[,inter[1]]) |> select(-1)
        } else {
          inter_df_vct <- env[,inter[1]]
        }

        # if the second variable is character or factor, create dummy cols, else just take it as it is
        if (inter_class[2] %in% c('character', 'factor')) {
          inter_df_fct <- fastDummies::dummy_cols(env[,inter[2]]) |> select(-1)
        } else {
          inter_df_vct <- env[,inter[2]]
        }

        # multiply the results and give them some names
        current_df_vct <- as_tibble(as.vector(inter_df_vct) * as.data.frame(inter_df_fct))
        names(current_df_vct) <- paste(names(inter_df_vct), names(inter_df_fct), sep = ":")
  }
      
  # bind resulting tables together   
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

    # get LC scores of sites
    current_df_fct <- scores(dbrda1,
                             choices = 1:2, # <- here, you must change to what you want to calculate
                             scaling = 'non', # <- here, you must change to what you want to calculate
                             correlation = F, # <- here, you must change to what you want to calculate
                             hill = F, # <- here, you must change to what you want to calculate
                             const = c(2,2), # <- here, you must change to what you want to calculate
                             tidy = T) |>
      as_tibble() |>
      filter(score == 'constraints') |>
      bind_cols(env) |>
      group_by(pick(all_of(c(var1_name, var2_name)))) |>
      summarise(CAP1 = mean(CAP1, na.rm = T),
                CAP2 = mean(CAP2, na.rm = T)) |>
      ungroup() |>
      unite({{final_col_name}}, where(is.character), sep = ':') |>
      relocate(c(CAP1, CAP2), .before = 1) |>
      pivot_longer(-where(is.numeric), names_to = 'interacting_variables', values_to = 'level_combinations')
    }
   
   # bind together table with centroids
   if (is.null(interaction_df_fct)) {
     interaction_df_fct <- current_df_fct
     } else {
       interaction_df_fct <- bind_cols(interaction_df_fct, current_df_fct)
       }
   
   } # end of if() that starts before forloop
   
# results
interaction_df_vct # this goes to envfit
interaction_df_fct # this is complete (centroids for interacting factor levels)

  


# INTERACTIONS including vectors -------------------------------------


### scaling non ---------------------------------------------------------

scores(dbrda1, scaling = 'non', const = c(2,2), correlation = F, tidy = T) |> 
  as_tibble() |> 
  filter(score == 'biplot')
  # here, look especially on annual_temperature:carbon, carbon:alt_classlow

ef <- envfit(dbrda1,
       env = interaction_df_vct,
       display = 'lc',
       scaling = 'non',
       choices = 1:2,
       correlation = F,
       const = c(2,2),
       hill = T) 
scores(ef, display = 'bp')

# ordiplot check
plot(dbrda1, display = c('bp', 'sp', 'si'), scaling = 'non')
text(dbrda1, display = 'bp', scaling = 'non')
plot(ef)

# THIS WORKS



### scaling species -----------------------------------------------------
scores(dbrda1, scaling = 'spe', const = c(2,2), correlation = F, tidy = T) |> 
  as_tibble() |> 
  filter(score == 'biplot')

ef <- envfit(dbrda1,
       env = interaction_df_vct,
       display = 'lc',
       scaling = 'spe',
       choices = 1:2,
       correlation = F,
       const = c(2,2),
       hill = T) 
scores(ef, display = 'bp')

# ordiplot check
plot(dbrda1, display = c('bp', 'sp', 'si'), scaling = 'spe')
text(dbrda1, display = 'bp', scaling = 'spe')
plot(ef)

# THIS ALSO WORKS, the numbers are the same as in the previous case



### scaling sites ------------------------------------------------------
scores(dbrda1, scaling = 'si', const = c(2,2), correlation = F, tidy = T) |> 
  as_tibble() |> 
  filter(score == 'biplot')

# 3 -0.254  -0.0503 biplot annual_temperature:carbon 
# 4 -0.237  -0.0777 biplot carbon:alt_classlow  

# what envfit gives
ef <- envfit(dbrda1,
       env = interaction_df_vct,
       display = 'lc',
       scaling = 'si', # this does not work, but we already know it
       choices = 1:2,
       correlation = F,
       const = c(2,2),
       hill = T) 
scores(ef, display = 'bp')

# annual_temperature:carbon -0.8065761 -0.4345389
# carbon:alt_class_low      -0.6906083 -0.6149892

# ordiplot check
plot(dbrda1, display = c('bp', 'sp', 'si'), scaling = 'si')
text(dbrda1, display = 'bp', scaling = 'si')
plot(ef)

# THESE DO NOT MATCH


#### RECALTULATION 
ef <- envfit(dbrda1,
       env = interaction_df_vct,
       display = 'lc',
       scaling = 'non',
       choices = 1:2,
       correlation = F,
       const = c(2,2),
       hill = T) 

# for rda, scaling sites should be:
# orthonormal species scores * constant
spe_scores <- scores(ef, display = 'bp') 
constant <- 2

tibble(CAP1 = spe_scores[,1] * constant,
       CAP2 = spe_scores[,2] * constant,
       label = rownames(spe_scores))

# # A tibble: 3 × 3
#     CAP1   CAP2 label                    
#    <dbl>  <dbl> <chr>                    
# 1 -1.61  -0.869 annual_temperature:carbon
# 2  0.284  1.20  carbon:alt_class_high    
# 3 -1.38  -1.23  carbon:alt_class_low  

# THIS IS COMPLETELY OFF


### scaling symmetric ------------------------------------------------
scores(dbrda1, scaling = 'sym', const = c(2,2), correlation = F, tidy = T) |> 
  as_tibble() |> 
  filter(score == 'biplot')

# 3 -0.470  -0.120  biplot annual_temperature:carbon 
# 4 -0.440  -0.185  biplot carbon:alt_classlow   

# what envfit gives
ef <- envfit(dbrda1,
       env = interaction_df_vct,
       display = 'lc',
       scaling = 'sym', # this does not work, but we already know it
       choices = 1:2,
       correlation = F,
       const = c(2,2),
       hill = T) 
scores(ef, display = 'bp')

# annual_temperature:carbon -0.8448748 -0.3543664
# carbon:alt_class_low      -0.7599692 -0.5268761

# ordiplot check
plot(dbrda1, display = c('bp', 'sp', 'si'), scaling = 'sym')
text(dbrda1, display = 'bp', scaling = 'sym')
plot(ef)

# THESE DO NOT MATCH

#### RECALTULATION 
ef <- envfit(dbrda1,
       env = interaction_df_vct,
       display = 'lc',
       scaling = 'non',
       choices = 1:2,
       correlation = F,
       const = c(2,2),
       hill = T) 

# for rda, scaling sites should be:
# orthonormal species scores * (lambda / sum(lambda))^(1/4) * constant
spe_scores <- scores(ef, display = 'bp') 
lambda <- eigenvals(dbrda1) # I AM NOT SURE, IF THIS IS WHAT lambda really means
sum_lambda <- sum(lambda)   # I AM NOT SURE, IF THIS IS WHAT lambda really means
constant <- 2

tibble(label = rownames(spe_scores),
       CAP1 = spe_scores[,1] * (lambda[1]/sum_lambda)^(1/4) * constant,
       CAP2 = spe_scores[,2] * (lambda[2]/sum_lambda)^(1/4) * constant)

# # A tibble: 3 × 3
#   label                       CAP1   CAP2
#   <chr>                      <dbl>  <dbl>
# 1 annual_temperature:carbon -0.941 -0.239
# 2 carbon:alt_class_high      0.242  0.482
# 3 carbon:alt_class_low      -0.879 -0.369

# THIS IS COMPLETELY OFF



### correlation = T ---------------------------------------------------------
scores(dbrda1, scaling = 'sym', const = c(2,2), correlation = F, tidy = T) |> 
  as_tibble() |> 
  filter(score == 'biplot')
# I don't know what to do with this one, because there is some recaltulation of scores
# for when correlation = T, but it does not change the numbers for biplot (in any scaling)

# correlation = F
# 3 -0.470  -0.120  biplot annual_temperature:carbon 
# 4 -0.440  -0.185  biplot carbon:alt_classlow 

# correlation = T
# 3 -0.470  -0.120  biplot annual_temperature:carbon 
# 4 -0.440  -0.185  biplot carbon:alt_classlow 


# INTERACTIONS only factors -------------------------------------------------

# this table was generated in the forloop
interaction_df_fct
# different scaling must be specified in the for loop!!!
# for now, correlation is set to F

### scaling non -------------------------------------------------------------
scores(dbrda1, scaling = 'non', const = c(2,2), correlation = F, tidy = T) |> 
  as_tibble() |> 
  filter(score == 'biplot')
# 5  0.0653  0.137 biplot alt_classlow:ph_classbasic
# 6 -0.261  -0.117 biplot alt_classlow:ph_classmesic

interaction_df_fct
# 5  0.0150  0.0316 alt_class:ph_class    low:basic         
# 6 -0.0472 -0.0213 alt_class:ph_class    low:mesic 

# Here, the numbers are the same for correlation T/F


### scaling species ---------------------------------------------------------
scores(dbrda1, scaling = 'spe', const = c(2,2), correlation = F, tidy = T) |> 
  as_tibble() |> 
  filter(score == 'biplot')
# 5  0.0653  0.137 biplot alt_classlow:ph_classbasic
# 6 -0.261  -0.117 biplot alt_classlow:ph_classmesic

interaction_df_fct
# correlation = F
# 5  0.0301  0.0632 alt_class:ph_class    low:basic         
# 6 -0.0943 -0.0425 alt_class:ph_class    low:mesic 

# correlation = T
# 5  0.0150  0.0316 alt_class:ph_class    low:basic         
# 6 -0.0472 -0.0213 alt_class:ph_class    low:mesic 

# The scores from scores are the same as in scaling 'non',
# but the centroids from site LC scores are different, 
# apart from that they do not match the result from scores...


### scaling sites ----------------------------------------------------------
scores(dbrda1, scaling = 'si', const = c(2,2), correlation = F, tidy = T) |> 
  as_tibble() |> 
  filter(score == 'biplot')
# 5  0.0191  0.0243 biplot alt_classlow:ph_classbasic
# 6 -0.0761 -0.0208 biplot alt_classlow:ph_classmesic

interaction_df_fct
# correlation = F
# 5  0.00878  0.0112  alt_class:ph_class    low:basic         
# 6 -0.0275  -0.00752 alt_class:ph_class    low:mesic 

# correlatoin = T
# 5  0.00878  0.0112  alt_class:ph_class    low:basic         
# 6 -0.0275  -0.00752 alt_class:ph_class    low:mesic 


# Same story, no match.
# But the correlation argument plays no role here.


### scaling symmetric -------------------------------------------------------
scores(dbrda1, scaling = 'sym', const = c(2,2), correlation = F, tidy = T) |> 
  as_tibble() |> 
  filter(score == 'biplot')
# 5  0.0353  0.0577 biplot alt_classlow:ph_classbasic
# 6 -0.141  -0.0494 biplot alt_classlow:ph_classmesic

interaction_df_fct
# correlation = F
# 5  0.0163  0.0266 alt_class:ph_class    low:basic         
# 6 -0.0510 -0.0179 alt_class:ph_class    low:mesic

# correlation = T
# 5  0.0163  0.0266 alt_class:ph_class    low:basic         
# 6 -0.0510 -0.0179 alt_class:ph_class    low:mesic         

# Same story, no match.
# But the correlation argument plays no role here.


