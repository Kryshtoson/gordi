#' =====================================================================
#' test run input
library(ggrepel)
library(tidyverse)
library(vegan)

spe <- read_csv('data/schrankogel/schrankogel_spe.csv')[-1] |> 
    log1p()

env <- read_csv("data/schrankogel/schrankogel_env.csv") |>
    mutate(group = elevation > 2500)

m <- capscale(spe ~ 1, distance = 'bray')



#' =====================================================================
#' arguments: 
#' -> model
#' -> headers 

gordi_read <- function(m, env, scaling = 'symm', correlation = F, hill = F){
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

    pass <- list(
        m = input$m,
        explained_variation = input$explained_variation,
        site_scores = input$site_scores,
        species_scores = input$species_scores,
        env = input$env
    )

    names(pass$site_scores) <- paste0("Axis", 1:2)
    
    #' here we should control for type of model
    if (class(m)[[1]] == "capscale") {
        actual_labs <- paste0("PCoA", 1:2, " (", round(pass$explained_variation[1:2]*100, 2), '%)')
    }

    p <- bind_cols(env, pass$site_scores) |> ggplot(aes(Axis1, Axis2)) +
        theme_bw() +
        labs(x = actual_labs[1],
        y = actual_labs[2]) +
        theme(
            text = element_text(size = 15),
            panel.grid = element_blank(),
            legend.justification = c(1, 1)
        )

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


gordi_read(m, env) |>
    gordi_sites(label = 'logger_ID', colouring = 'group')

