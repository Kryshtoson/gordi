#' =====================================================================
#' test run input
library(ggrepel)
library(tidyverse)
library(vegan)

spe <- read_csv('data/schrankogel/schrankogel_spe.csv')[-1] |> 
    log1p()

env <- read_csv("data/schrankogel/schrankogel_env.csv") |>
    mutate(group = elevation > 2500)


m <- rda(spe ~ 1, distance = "bray") #cca
# cca(spe~ elevation, data = env)
# m <- capscale(spe ~ 1, distance = 'bray')
# m <- capscale(spe ~ elevation + annual_temperature , distance = 'bray', data = env)

m <- decorana(spe)

as_tibble(scores(m))

scores(m)$sites

is.null(m$CCA[1])

class(m)

#' cca {CCA, partial CCA, CA, partial CA}
#' capscale {dbRDA, partial dbRDA, PCoA, partial PCoA}
#' rda {RDA, partial RDA, PCA, partial PCA}
#' ...


#' =====================================================================
#' arguments: 
#' -> model
#' -> headers 


gordi_read <- function(m, env, scaling = 'symm', correlation = F, hill = F){

    if(class(m)[[1]]  %in% )

    if(is.null(m$CCA)){
        #'... 
    }
    
    pass <- list(
        m = m,
        explained_variation = m$CA$eig/m$tot.chi,
        site_scores = as_tibble(as.data.frame(scores(m, scaling = scaling, correlation = correlation, hill = hill)$sites)),
        species_scores = as_tibble(as.data.frame(scores(m, scaling = scaling, correlation = correlation, hill = hill)$species)),
        env = env
    )
    
    
    return(pass)

}

gordi_sites <- function(input, label = '', colouring = '', repel_label = T) {

    input #' misto pass
    # pass <- list(
    #     m = input$m,
    #     explained_variation = input$explained_variation,
    #     site_scores = input$site_scores,
    #     species_scores = input$species_scores,
    #     env = input$env
    # )

    names(pass$site_scores) <- paste0("Axis", 1:2)
    
    #' here we should control for type of model
    #' capscale
    if (class(pass$m)[[1]] == "capscale") {
        if(){} #' dbRDA osetrit
        else{
        actual_labs <- paste0("PCoA", 1:2, " (", round(pass$explained_variation[1:2]*100, 2), '%)')
        }
    }
    #' cca
    else if (class(pass$m)[[1]] == "cca") {
        actual_labs <- paste0("CCA", 1:2, " (", round(pass$explained_variation[1:2]*100, 2), '%)')
    }
    #' pca
    else if (class(pass$m)[[1]] == "pca") {
        actual_labs <- paste0("PCA", 1:2, " (", round(pass$explained_variation[1:2]*100, 2), '%)')
    }

    if (exists(input$p)) { #' is.null(input$p) #' mozna lepsi
        p <- bind_cols(env, pass$site_scores) |>
            ggplot(aes(Axis1, Axis2)) +
            theme_bw() +
            labs(x = actual_labs[1], y = actual_labs[2]) +
            theme(
                text = element_text(size = 15),
                panel.grid = element_blank(),
                legend.justification = c(1, 1)
            )
    } else {
        p #'... 
    }

    #' accounting for colouring
    if (colouring == '') {
        p <- p +
            geom_point(size = 3)
    } else {
        p <- p +
            geom_point(aes(colour = !!sym(colouring)), size = 3)
    }

    #' accounting for labeling
    if (label != '') {
        if (repel_label) {
            p <- p + geom_text_repel(aes(label = !!sym(label)))
        } else {
            p <- p + geom_text(aes(label = !!sym(label)))
        }
    }

    pass$plot <- list(p)
    
    return(pass)
}


obj <- gordi_read(m, env) |>
    gordi_species()
    
    
    
    #' gordi_sites(label = "logger_ID", colouring = "group") |>

obj

str(obj)
