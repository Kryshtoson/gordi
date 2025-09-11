# actual labs

m <- rda(spe ~ elevation + slope, data = env)

o <- gordi_read(m, choices = c(1:3))

o

names(as_tibble(as.data.frame(scores(m, scaling = 'symm', choices = 1:4, correlation = F, hill = T)$species)))


o <- gordi_read(m, choices = choices)

choices <- 1:3
  
o <- gordi_read(m, choices = choices)

m <- metaMDS(spe, k = 30)

as_tibble(scores(m, scaling = 'symm', display = 'species', choices = choices, correlation = F, hill = T))

paste0(o$axis_names, " (", round(o$explained_variation[choices]*100, 2), "%)")


m$points

# new version
if(pass$type %in% c('DCA', 'NMDS')) {actual_labs <- paste0(pass$axis_names)} else 
  {actual_labs <- paste0(pass$axis_names, " (", round(pass$explained_variation[choices]*100, 2), "%)")}

if(o$type %in% c('DCA', 'NMDS')) {actual_labs <- paste0(o$axis_names)} else 
{actual_labs <- paste0(o$axis_names, " (", round(o$explained_variation[choices]*100, 2), "%)")}
actual_labs

m <- decorana(spe)

'decorana'

as.character(m$call)[1]
word(as.vector(m$call), 1)

class(m)

m <- metaMDS(spe)

class(m)




gordi_read(m, choices = 2:3) |> 
  gordi_species()



# original ifelse based on ordination type, which apparently does not work
if(pass$type == 'CCA'){
  actual_labs <- paste0("CCA", pass$choices, " (", round(pass$explained_variation[1:2]*100, 2), '%)')
} else if(pass$type == 'CA'){
  actual_labs <- paste0("CA", pass$choices, " (", round(pass$explained_variation[1:2]*100, 2), '%)')
} else if (pass$type == 'PCA'){
  actual_labs <- paste0("PCA", pass$choices, " (", round(pass$explained_variation[1:2]*100, 2), '%)')
} else if(pass$type == 'PCoA'){
  actual_labs <- paste0("PCoA", pass$choices, " (", round(pass$explained_variation[1:2]*100, 2), '%)')
} else if(pass$type == 'db-RDA'){
  actual_labs <- paste0("db-RDA", pass$choices, " (", round(pass$explained_variation[1:2]*100, 2), '%)')
} else if(pass$type == 'RDA'){
  actual_labs <- paste0("RDA", pass$choices, " (", round(pass$explained_variation[1:2]*100, 2), '%)')
} else if (pass$type == 'DCA'){
  actual_labs <- paste0('DCA', pass$choices)
} else if (pass$type == 'NMDS'){
  actual_labs <- paste0('NMDS', pass$choices)
}
