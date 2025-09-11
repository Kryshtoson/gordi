#gordi shape

gordi_shape <- function(pass,
                        scale = c('discrete', 'identity'),
                        values = NULL,
                        breaks = waiver(),
                        name = waiver(),
                        labels = waiver(),
                        limits = NULL,
                        guide = waiver(),
                        drop = TRUE,
                        na.translate = TRUE) {
  
  scale <- match.arg(scale)
  sc <- function(suffix) paste0('scale_shape_', suffix)

  
  build_scale_obj <- function(){
    if (scale == 'identity'){
      args <- purrr::compact(list(name = name, breaks = breaks, labels = labels, limits = limits, guide = guide))
      return(do.call(sc('identity'), args))
    }
    
    
    if(!is.null(values)){
      args <- purrr::compact(list(values = values, name = name, breaks = breaks, labels = labels, limits = limits, guide = guide, drop = drop, na.translate = na.translate))
      return(do.call(sc('manual'), args))
    }
    else {
      args <- purrr::compact(list(name = name, breaks = breaks, labels = labels, limits = limits, guide = guide, drop = drop, na.translate = na.translate))
      return(do.call(sc('discrete'), args))
    }
  }
  
  scale_obj <- build_scale_obj()
  pass$plot <- pass$plot + scale_obj
  
  pass$plot <- pass$plot + ggnewscale::new_scale('shape')
  
  return(pass)
  
}

gordi_read(m, env, traits)|>
  gordi_species(label = T, colour = 3, shape = 'logger_ID')|>
 # gordi_label(what = 'sites', label_colour = 2)
  gordi_shape(scale = 'discrete', values = c(21,20, 12, 11, 10, 19))

m <- cca(spe~1)

o <- gordi_read(m, env)

bind_cols(o$site_scores, o$env)|>
  ggplot(aes(CAP1, MDS1))+
  geom_point(aes(shape = elevation))+
  scale_shape_binned()
