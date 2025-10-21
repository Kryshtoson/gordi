library(vegan)
library(patchwork)
library(readxl)
library(tidyverse)
remotes::install_github('Kryshtoson/gordi')
library(gordi)

# 1 PCoA Auchenorrhyncha

auch.spe <- read_xlsx("data/LimestoneQuarries_Aucheno/Aucheno_sum.xlsx") |> #gordi can handle various species names, but unfortunately this is beyond its power
  rownames_to_column(var = 'id')|>
  pivot_longer(cols = -id, names_to = 'spe_name', values_to = 'value')|>
  mutate(spe_name = str_c(str_split_i(spe_name, "\\s", 1) %>% str_sub(., 1, 3), 
                          str_split_i(spe_name, "\\s", 2) %>% str_sub(., 1, 3), sep = ".")) |>
  pivot_wider(names_from = spe_name, values_from = value) %>%
  select(-id)|>
  select(where(~sum(.x > 0) > 1)) |> # Removal of singular occurrences
  log1p()|> #log transformation
  as_tibble()

auch.env <- read_xlsx("data/LimestoneQuarries_Aucheno/Aucheno_sum.xlsx", 2) 

# B-C PCoA
pco.bc <- capscale(auch.spe~1, distance = "bray", sqrt.dist = T, binary = F) #in Jakub script he uses method = 'bray', but shouldnt it be distance ='bray'?
pco.bc

eigenvals(pco.bc)/sum(eigenvals(pco.bc)) #implemented in gordi
screeplot(pco.bc)
screeplot(pco.bc, bstick = T)

# site visualization
gordi_read(pco.bc, env = auch.env)|>
  gordi_sites(colour = 'Treatment')|>
  gordi_colour(scale = 'discrete', family = 'manual', values = c('darkorange', 'darkgreen'))|>
  gordi_cluster(group = 'site', linetype = 'dotted') #argument group only draws lines, input must be a column from env dataframe

#species visualization

gordi_read(pco.bc, spe = auch.spe, env = auch.env, correlation = T)|>
  gordi_fit(slice_max = 30)|>
  gordi_species()|>
  gordi_label(what = 'species', label_colour = 4, repel_label = T, max.overlaps = Inf)

#together

gordi_read(pco.bc, spe = auch.spe, env = auch.env)|> #const can be set for better visualization of both species and sites
  gordi_sites(colour = 'Treatment')|>
  gordi_colour(scale = 'discrete', family = 'manual', values = c('darkorange', 'darkgreen'))|>
  gordi_cluster(group = 'site', linetype = 'dotted', show.legend = F)|>
  gordi_fit(slice_max = 30)|>
  gordi_species()|>
  gordi_label(what = 'species', label_colour = 4, repel_label = T, max.overlaps = Inf, show.legend = F)


# Sorensen PCoA -----------------------------------------------------------

auch.pa <- auch.spe %>%
  mutate(across(1:ncol(.), ~ifelse(. > 0, 1, 0))) # replace values > 0 with 1

pco.so <- capscale(auch.pa~1, distance = "bray", sqrt.dist = T, binary = F) #carefull, use distance
pco.so

eigenvals(pco.so)/sum(eigenvals(pco.so)) #implemented in gordi
screeplot(pco.so)

#together

gordi_read(pco.so, spe = auch.spe, env = auch.env, correlation = T)|>
  gordi_sites(colour = 'Treatment')|>
  gordi_colour(scale = 'discrete', family = 'manual', values = c('darkorange', 'darkgreen'))|>
  gordi_cluster(group = 'site', linetype = 'dotted', show.legend = F)|>
  gordi_fit(slice_max = 30)|>
  gordi_species(colour = 'red', alpha = 0.5)|>
  gordi_label(what = 'species', label_colour = 2, repel_label = T, max.overlaps = Inf, alpha = 0.5, show.legend = F)


# 2 Chironomids NMDS ------------------------------------------------------

chiro <- read_xlsx("data/svratka_chironomids/chironomids.xlsx") |>
  mutate(across(c(-sample), ~replace_na(., 0))) |>
  select(-sample)

chiro.env <- read_excel("data/svratka_chironomids/chironomids.xlsx", 2)

# Fitting NMDS for increasing numbers of ordination axes
nmds.2<-metaMDS(log1p(chiro), k = 2)
nmds.3<-metaMDS(log1p(chiro), k = 3)
nmds.4<-metaMDS(log1p(chiro), k = 4)
nmds.5<-metaMDS(log1p(chiro), k = 5) #5 dimensions- but does it really help- no

# Saving stress values in a df for plotting
stress.df <- tibble(k = 2:5, stress = c(nmds.2$stress, nmds.3$stress, 
                                        nmds.4$stress, nmds.5$stress))

#Plot of stress values against k
ggplot(stress.df, aes(k, stress))+
  geom_point()+
  geom_line(lty = 2)+
  theme_classic()#
# It seems that k = 3 is a reasonable option. Stress is lower than 0.1 indicating a good fit and
# k is still quite low

# Stressplot fot NMDS with k = 3
stressplot(nmds.3)

#1st two ordination axis
gordi_read(nmds.3, spe = chiro, env = chiro.env)|>
  gordi_cluster(cluster = 'hydr', spider = T, colour = 'hydr', label = T, show.legend = F) #so far argument hull does not work, input for cluster must be a column from env dataframe

#1st and 3rd ordination axis

gordi_read(nmds.3, spe = chiro, env = chiro.env, choices = c(1, 3))|>
  gordi_cluster(cluster = 'hydr', spider = T, colour = 'hydr', label = T, show.legend = F)

#species plot

gordi_read(nmds.3)|>
  gordi_label(what = 'species', shortcut = 'upper.lower', shortcut_colour = 'blue', repel_label = T, max.overlaps = Inf, )



# PCoA repeated vegetation plots ------------------------------------------

spe <- read_csv("data/S_Moravia_basiphilous_grass_KK/basiphilous_grasslands_S_Moravia_species.csv") |> #gordi can handle this species names, but when we want to colour fitted species based on traits, species names from traits dataframe must match species names from spe dataframe
  pivot_longer(cols = -Releve_number,names_to = 'spe', values_to = 'value')|>
  mutate(spe = gsub("_[0-9]+$", "", spe)) |>
  pivot_wider(names_from = spe, values_from = value)|>
  select(-Releve_number) %>%
  select(where(~sum(.x > 0) > 1)) |> # remove species with one occurrence
  sqrt() # square root transformation

# species indicator values, Red List status, specialization
spe.data <- read_csv("data/S_Moravia_basiphilous_grass_KK/basiphilous_grasslands_S_Moravia_species_data.csv")

# header data 
head <- read_csv("data/S_Moravia_basiphilous_grass_KK/basiphilous_grasslands_S_Moravia_head.csv")|>
  mutate(Rs_observ = factor(Rs_observ)) #gordi_shape cannot handle continuous variables, bcs it is build on scale_shape_manual()

pcoa <- capscale(spe ~ 1, distance = "bray", sqrt.dist = T)
pcoa

eigenvals(pcoa)/sum(eigenvals(pcoa))
screeplot(pcoa, bstick = T)

# ordination plot with sites
gordi_read(pcoa, spe = spe, env = head)|>
  gordi_sites(fill = 'Rs_observ', shape = 'Rs_observ')|>
  gordi_shape(scale = 'discrete', values = 21:22)|>
  gordi_colour(scale = 'discrete', family = 'manual', fill = T, values = c('red', 'green'))|>
  gordi_cluster(group = 'Rs_plot', linewidth = 0.3, arrow = T, arrow_type = 'open', show.legend = F)

#ordination plot with species, coloured by dry grassland specialists

gordi_read(pcoa, spe = spe, env = head, scaling = 'si', correlation = T, traits = spe.data)|>
  gordi_fit(slice_max = 40)|>
  gordi_species(colour = 'THE_THF')|>
  gordi_colour(scale = 'discrete', family = 'manual', values = c('red', 'blue'))|>
  gordi_label(what = 'species', shortcut = 'upper.lower', shortcut_colour = 'THE_THF', repel_label = T, max.overlaps = Inf, show.legend = F)
