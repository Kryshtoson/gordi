# testing

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


# Running functions on various ordinations -------------------------------------

### PCA
m <- rda(spe ~ 1)
gordi_read(m) |> 
  gordi_sites() |> 
  gordi_species() 

### RDA
m <- rda(spe ~ elevation, data = env)
gordi_read(m, env) |> 
  gordi_sites() |> 
  gordi_predict()

### CA
m <- cca(spe ~ 1)
gordi_read(m) |> 
  gordi_sites() |> 
  gordi_species(label = F) |> 
  gordi_predict()

### CCA
m <- cca(spe ~ elevation, data = env)
gordi_read(m) |> 
  gordi_sites() |> 
  gordi_species(label = F) |> 
  gordi_predict() # it is okay to display continuous predictor as arrow?

### DCA
m <- decorana(spe)
gordi_read(m) |> 
  gordi_sites() |> 
  gordi_species(label = F) |> 
  gordi_predict()

### PCoA
m <- capscale(spe ~ 1, data = env)
gordi_read(m) |> 
  gordi_sites() |> 
  gordi_species(label = F) |> 
  gordi_predict()

### db-RDA
m <- capscale(spe ~ elevation, distance = 'bray', data = env)

s <- summary(m)
s$inertia

m$tot.chi

m$CCA$eig
m$CA$eig 
summary(m$CA$eig)
summary(m$CCA$eig)

gordi_read(m, scaling = 'symm') |> 
  gordi_sites() |> 
  gordi_species(label = F) |> 
  gordi_predict(scaling_coefficient = 0.6)

### NMDS
m <- metaMDS(spe, k = 5)

gordi_read(m, scaling = 'symmetric') |> 
  gordi_sites() |> 
  gordi_species(label = F) |> 
  gordi_predict()
