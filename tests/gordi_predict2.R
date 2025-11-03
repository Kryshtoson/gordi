#' Extracts predictor scores and relevant stuff and creates plot - so far working for capscale. EXPERIMENTAL VERSION
#' 
#' @description
#' THIS IS EXPERIMENTAL VERSION
#' 
#' [gordi_predict()] takes the result of [gordi_read()] and creates plot with
#' predictor arrows (for continuous variables) and points (for categorical
#' variables, i.e., centroids). So far, it can work with continuous explanatory 
#' variables. It also can't display explanatory variables passively (i.e., not 
#' involved in the ordination model).
#' 
#' There will be a special function for this once, maybe something as
#' `gordi_corr()`. Similarly to [gordi_species()] and [gordi_sites()],
#' you can also set a wide range of graphing parameters, such as colour, fill, size,
#' shape, alpha, stroke, and more traditional ggplot arguments, which can read
#' both, static and dynamic variable (e.g., 'red' or 'elevation').
#' 
#' @param pass Object from [gordi_read()] function.
#' @param label Logical; default = F, whether to display label by each point/arrow 
#'   or not. Labels use the `predictor` name for arrows and the combined 
#'   `predictor: level` for centroids. If you want to customize the labels, you can 
#'   with [gordi_label()] function, which overrides this setting.
#' @param show_label Logical; default = F. A simplified way to display labels. Overridden by [gordi_label()].
#' @param repel_label Logical; repels labels of predictors for better readability.
#'    Default is F. Overridden by [gordi_label()].
#' @param colour Colour can be defined statically as word from the [colours()]
#'   list (e.g. 'red'), HEX code (e.g. #5d782e), or number from [palette()] (e.g. 3).
#'   It can also defined dynamically (according to some variable, e.g. elevation
#'   or vegetation type). This variable has to be present in the `env` table and its
#'   name has to be written in "quotation marks". Default colour of arrows is 2, 
#'   and 2 for points.
#' @param alpha Transparency of symbols/arrows. Numeric value. Default is 1 (for arrows)
#'   and 1 (for points). Can be set statically or dynamically (use "").       
#' @param arrow_size Numeric; length of arrow in cm. Default is 0.3 cm. 
#' @param linewidth Self-explanatory. Can be numerical value. Used for arrows. See
#'    [aes_linetype_size_shape()] for more details. Default is 1.
#' @param linetype Changes type of line used in arrows. See [aes_linetype_size_shape()] for more details. Default is 1 (solid).
#' @param fill Fill colour for centroid points. Can be defined statically or dynamically (e.g., "score"). Default is 2.
#' @param shape Shape of centroid points. See [aes_linetype_size_shape()] for more details. Default is 16.
#' @param size Size of centroid points. Numeric value. Can be set statically 
#'    or dynamically (e.g., "r2"). Default is 2.
#' @param scaling_coefficient Numeric; edits the lengths of predictor arrows. Default value
#'    is `0.9` which means that the predictor arrow will be scaled to 0.9 of
#'    the longest axis displayed in the plot frame.
#' 
#' 
#' @return The input object `pass` with an updated `plot` element that includes
#'   predictor arrows and labels.
#' 
#' 
#' @importFrom rlang is_empty sym syms expr
#' @importFrom stringr str_detect regex str_extract str_remove str_split
#' @importFrom dplyr filter mutate bind_cols rename select pull case_when ungroup rowwise distinct
#' @importFrom tidyr unite
#' @importFrom purrr discard keep imap_dfr map_chr
#' @importFrom grDevices colours palette
#' @importFrom ggplot2 ggplot theme_bw labs theme element_text element_blank element_rect aes geom_segment arrow unit geom_point ggplot_build geom_text
#' @importFrom ggrepel geom_text_repel
#' @importFrom utils hasName
#' @importFrom tidyselect last_col all_of
#' @importFrom fastDummies dummy_cols
#' @importFrom vegan envfit scores
#' @importFrom tibble tibble
#' 
#'
#' @examples 
#' library(vegan)
#' 
#' data(dune)
#' data(dune.env)
#' 
#' # --- 1. Example with main effects ---
#' m1 <- capscale(dune ~ A1 + Management, data = dune.env)
#' gordi_read(m1, env = dune.env, scaling = 'species', correlation = T) |>
#'    gordi_species(label = F) |>
#'    # Predictors (A1: continuous arrow, Management: categorical centroids)
#'    # Colour is dynamically mapped to the 'score' column (biplot or centroid)
#'    gordi_predict2(scaling_coefficient = 1, colour = 'score', size = 4)
#' 
#' # --- 2. Example with an interaction term ---
#' # The interaction term will be calculated post-hoc via envfit
#' m2 <- capscale(dune ~ A1 * Management, data = dune.env)
#' gordi_read(m2, env = dune.env, scaling = 'species', correlation = T) |>
#'    gordi_sites() |>
#'    # Predictors include main effects and interaction effects (e.g., A1:Management)
#'    gordi_predict2(show_label = T, repel_label = T, colour = 'score', size = 3)  
#' @export
gordi_predict2 <- function(
    pass,
    scaling_coefficient = 0.9,
    label = c('label'),
    show_label = F,
    repel_label = F,
    colour = '',
    fill = '',
    alpha = '',
    linetype = '',
    linewidth = '',
    shape = '',
    size = '',
    arrow_size = 0.3) {
  
  if (!pass$type %in% c('RDA', 'CCA', 'db-RDA') || rlang::is_empty(pass$predictor_scores)) {
    stop("You provided an unconstrained ordination. Predictors can't be displayed. If you want to passively plot environmental variables, try `gordi_corr()`.")
  }
  
  ### ordination axis labels
  if(pass$type %in% c('DCA', 'NMDS')) {actual_labs <- paste0(pass$axis_names)} else 
  {actual_labs <- paste0(pass$axis_names, " (", round(pass$explained_variation[pass$choices]*100, 1), "%)")}
  
  ### which variables do what
  term_labels <- attr(terms(pass$m), "term.labels")
  
  inter_terms <- term_labels[str_detect(term_labels, ':')]
  main_terms <- term_labels[!str_detect(term_labels, ':')]
  
  
  ### --- MAIN TERMS ---
  main_terms_scores <- scores(pass$m,
                              scaling = pass$scaling,
                              choices = pass$choices,
                              correlation = pass$correlation,
                              hill = pass$hill,
                              const = pass$const,
                              tidy = T) |>
    as_tibble() |>
    filter(score %in% c('biplot', 'centroids') & !str_detect(label, ':'))
  
  
  # create NULL objects
  factor_scores <- NULL
  vector_scores <- NULL
  
  # create a tibble with centroid scores
  if (!is.null(main_terms_scores |> filter(score == 'centroids'))) {
    factor_scores <- main_terms_scores |> filter(score == 'centroids')
  }
  
  # create a tibble with arrow ends
  if (!is.null(main_terms_scores |> filter(score == 'biplot'))) {
    vector_scores <- main_terms_scores |> filter(score == 'biplot')
  }
  
  
  ### --- INTERACTION TERMS ---
  vector_inter_scores <- tibble::tibble() 
  factor_inter_scores <- tibble::tibble()
  
  if (length(inter_terms) > 0) {
    
    warning('All interaction scores are fitted post-hoc via envfit as `lc` scores.')
    if (is.null(pass$env)) {
      stop('If you want to display interactions, you must provide env table to the `gordi_read()`.')
      }
    
  # what interacts with what
    interaction_table <- inter_terms |>
      str_split(pattern = ":") |>
      set_names(~ paste0("inter_", seq_along(.))) |> 
      imap_dfr(~ tibble(
        inter_ID = .y,
        inter_var = paste0("var", seq_along(.x)),
        variable = .x
      )) |> 
      rowwise() |> 
      mutate(var_class = case_when(
        !(variable %in% names(pass$env)) ~ 'Not found',
        is.numeric(pass$env[[variable]]) | is.integer(pass$env[[variable]]) ~ 'vector',
        is.factor(pass$env[[variable]]) | is.character(pass$env[[variable]]) ~ 'factor',
        TRUE ~ class(pass$env[[variable]])[1])
      ) |> 
      ungroup()  

  }
  
    # --- FOR LOOP to obtain interaction scores ---
    interaction_df <- NULL #tibble(.rows = nrow(pass$env)) # here, interaction_df is created as an empty table, with the correct number of rows, but no columns

    if (length(inter_terms) > 0) {  
      
    for (i in pull(distinct(interaction_table, inter_ID))) {
    
      # filter out interaction terms
      inter <- interaction_table |> 
        filter(inter_ID == i) |> 
        select(variable) |> 
        pull()
    
      inter_class <- pass$env |>
        select(all_of(inter)) |> 
        map_chr(class) 
    
      # Variable to hold the result of the current iteration
      current_df <- NULL
    
      # --- NUMERIC x NUMERIC ---
      if (all(inter_class %in% c('numeric', 'integer', 'double'))) {
        current_df <- as_tibble(pass$env[,inter[1]] * pass$env[,inter[2]])
        colnames(current_df) <- paste(inter[1], inter[2], sep = ':')
        
        # --- NUMERIC x FACTOR/CHARACTER ---        
      } else if (any(inter_class %in% c('character', 'factor')) && 
                 any(inter_class %in% c('numeric', 'integer', 'double'))) {
        
        inter_df_vct <- NULL
        inter_df_fct <- NULL
      
          if (inter_class[1] %in% c('character', 'factor')) {
            inter_df_fct <- fastDummies::dummy_cols(pass$env[,inter[1]]) |> select(-1)
          } else {
            inter_df_vct <- pass$env[,inter[1]]
          }
      
          if (inter_class[2] %in% c('character', 'factor')) {
            inter_df_fct <- fastDummies::dummy_cols(pass$env[,inter[2]]) |> select(-1)
          } else {
            inter_df_vct <- pass$env[,inter[2]]
          }
      
      current_df <- as_tibble(as.vector(inter_df_vct) * as.data.frame(inter_df_fct)) 
      names(current_df) <- paste(names(inter_df_vct), names(inter_df_fct), sep = ":")
      
      # --- FACTOR/CHARACTER x FACTOR/CHARACTER ---    
      } else if (all(inter_class %in% c('character', 'factor'))) {
        var1_name <- inter[1]  
        var2_name <- inter[2]  
      
        final_col_name <- paste(var1_name, var2_name, sep = ":")
      
        current_df <- pass$env |>
          mutate(
            interaction_term = paste(
              paste0(var1_name, "_", pass$env[[var1_name]]),
              paste0(var2_name, "_", pass$env[[var2_name]]),
              sep = ":")) |>
          select(interaction_term) |> 
          setNames(final_col_name)
      }
      
    
        if (is.null(interaction_df)) {
          interaction_df <- current_df
        } else {
          interaction_df <- bind_cols(interaction_df, current_df)
        }
    
  }
  
  } # end of if() that starts before forloop
    
    
  # --- CACULATE ENVFIT ---
  if (length(inter_terms) > 0) {
    inter_ef <- envfit(pass$m, env = interaction_df, 
                       display = 'lc', 
                       scaling = pass$scaling,
                       choices = pass$choices,
                       correlation = pass$correlation)

    var1_name <- inter[1]
    var2_name <- inter[2]
        
    vector_inter_scores <- as_tibble(scores(inter_ef, display = 'bp'), rownames = 'label')
    factor_inter_scores <- as_tibble(scores(inter_ef, display = 'cn'), rownames = 'label')
    
    vector_inter_scores <- vector_inter_scores |>
      dplyr::mutate(
        level = stringr::str_remove(label, "^\\:?\\.?data\\_"), # Extracts just the level (e.g., "BF")
        label = paste0(var1_name, ":", var2_name, level) # Constructs 'A1:ManagementBF'
      ) |>
      dplyr::select(-level)
    
  }  

  
  
  # --- MERGE PREDICTOR TABLES together ---
  lst <- list(
    factor_scores = factor_scores |> mutate(score = "factor_scores"),
    vector_scores = vector_scores |> mutate(score = "vector_scores"),
    vector_inter_scores = if (!is.null(vector_inter_scores)) {
      vector_inter_scores |> mutate(score = "vector_inter_scores")
    } else {NULL},
    factor_inter_scores = if (!is.null(factor_inter_scores)) {
      factor_inter_scores |> mutate(score = "factor_inter_scores")
    } else {NULL}
  ) |>
    purrr::keep(~ !is.null(.x))
  
  predictor_scores <- lst |> 
    imap_dfr(~ .x |> mutate(score = .y))
  
  pass$predictor_scores <- predictor_scores |> 
    rename(Axis_pred1 = 1,
           Axis_pred2 = 2)
  
  
  
  # --- PLOTTING ---
  
  pred_df <- pass$predictor_scores
  
  if (is.null(pass$plot)) { # checks whether p exists in pass, if not it draws plot
    p <- ggplot() +
      theme_bw() +
      geom_vline(aes(xintercept = 0), linetype = 3, linewidth = 0.2, colour = 'gray15', alpha = 0.6) +
      geom_hline(aes(yintercept = 0), linetype = 3, linewidth = 0.2, colour = 'gray15', alpha = 0.6) +
      labs(x = actual_labs[1], y = actual_labs[2]) +
      theme(
        text = element_text(size = 15),
        panel.grid = element_blank(),
        legend.justification = c(1, 1))
  } else {p <- pass$plot}
  
  
  ### Detect mapped vs constant aesthetics
  # colour
  # if colour != "" AND ALSO colour represents a colname present in corr_df, then map_colour is TRUE, otherwise is FALSE
  map_colour <- !identical(colour, '') && (colour %in% names(pred_df))
  # if map_colour is FALSE AND ALSO the thing inputed in arguments is a HEX code or is included in colours() or in palette() (word or number), then use it as const_colour
  const_colour <- !map_colour && (grepl("^#(?:[A-Fa-f0-9]{6}[A-Fa-f0-9]{3})$", colour) || colour %in% grDevices::colours()) || (is.character(colour) && colour %in% palette()) || (is.numeric(colour) && colour %in% seq_along(palette()))
  # fill
  map_fill <- !identical(fill, '') && (fill %in% names(pred_df))
  const_fill <- !map_fill && (grepl("^#(?:[A-Fa-f0-9]{6}[A-Fa-f0-9]{3})$", fill) || fill %in% grDevices::colours()) || (is.character(fill) && fill %in% palette()) || (is.numeric(colour) && colour %in% seq_along(palette()))
  # alpha
  map_alpha <- !identical(alpha, '') && has_name(pred_df, alpha)
  const_alpha <- !map_alpha && is.numeric(alpha)
  # linetype
  map_linetype <- !identical(linetype, '') && has_name(pred_df, linetype)
  const_linetype <- !map_linetype && (is.numeric(linetype) || !identical(linetype,''))
  # linewidth
  map_linewidth <- !identical(linewidth, '') && has_name(pred_df, linewidth)
  const_linewidth <- !map_linewidth && is.numeric(linewidth)
  # shape
  map_shape <- !identical(shape, '') && has_name(pred_df, shape)
  const_shape <- !map_shape && is.numeric(shape)
  # size
  map_size <- !identical(size, '') && has_name(pred_df, size)
  const_size <- !map_size && is.numeric(size)
  
  
  ### Set scaling coefficient
  # extract plot frame size (x and y axis lengths)
  vct_tbl <- pred_df |> filter(str_detect(score, 'vector'))
  
  if (nrow(vct_tbl) > 0) {
    p_build <- ggplot_build(p)
    
    plot_range <- c(xmin_plot = p_build$layout$panel_params[[1]]$x.range[1],
                    xmax_plot = p_build$layout$panel_params[[1]]$x.range[2],
                    ymin_plot = p_build$layout$panel_params[[1]]$y.range[1],
                    ymax_plot = p_build$layout$panel_params[[1]]$y.range[2])
    
    pred_range <- c(xmin_pred = min(vct_tbl[,1], na.rm = T),
                    xmax_pred = max(vct_tbl[,1], na.rm = T),
                    ymin_pred = min(vct_tbl[,2], na.rm = T),
                    ymax_pred = max(vct_tbl[,2], na.rm = T))
    
    coef <- (max(abs(plot_range)) / max(abs(pred_range))) * scaling_coefficient
  } else {coef <- 1}
  
  
  ### Prepare aes arguments for geom_segment()
  # Start with fixed x/y for the base (0,0) and end at the species scores
  aes_args_segment <- list(
    x = 0, y = 0,
    xend = expr(Axis_pred1 * !!coef),
    yend = expr(Axis_pred2 * !!coef)
  )
  
  if(map_colour) aes_args_segment$colour <- sym(colour)
  if(map_alpha) aes_args_segment$alpha <- sym(alpha)
  if(map_linewidth) aes_args_segment$linewidth <- sym(linewidth)
  if(map_linetype) aes_args_segment$linetype <- sym(linetype)
  
  
  ### Prepare constant arguments for geom_point() (mapped first, then defaults if nothing)
  const_args_segment <- list()
  
  # Add constant arguments for geom_segment() if not mapped
  # colour
  if(!map_colour){
    if(!identical(colour, '')) {const_args_segment$colour <- colour} else {const_args_segment$colour <- 2}}
  # alpha
  if(!map_alpha){
    if(!identical(alpha, '')) {const_args_segment$alpha <- alpha} else {const_args_segment$alpha <- 1}}
  # linetype
  if(!map_linetype){
    if(!identical(linetype, '')) {const_args_segment$linetype <- linetype} else {const_args_segment$linetype <- 1}}
  # linewidth
  if(!map_linewidth){
    if(!identical(linewidth, '')) {const_args_segment$linewidth <- linewidth} else {const_args_segment$linewidth <- 1}}
  # arrow_size (does not make sense to have map_arrow)
  const_args_segment$arrow <- arrow(
    length = unit(
      if (!identical(arrow_size, '')) as.numeric(arrow_size) else 0.3, "cm"))
  
  
  ### Add ARROWS to plot  
  if (!is.null(pred_df |> filter(str_detect(score, 'vector')))) {
    p <- p + do.call(geom_segment,
                     c(list(data = pred_df |> filter(str_detect(score, 'vector')),
                            mapping = do.call(aes, aes_args_segment)),
                       const_args_segment)) 
    
    if (isTRUE(show_label)){
      
      aes_args_text <- list(
        x = expr(Axis_pred1 * !!coef),
        y = expr(Axis_pred2 * !!coef),
        label = sym(label)
      )
      
      if(map_colour) aes_args_text$colour <- sym(colour)
      if(map_alpha) aes_args_text$alpha <- sym(alpha)
      
      ### Prepare constant arguments for geom_text() (mapped first, then defaults if nothing)
      const_args_text <- list(
        show.legend = F
      )
      
      if(!map_colour){ if(!identical(colour, '')) {const_args_text$colour <- colour} else {const_args_text$colour <- 'gray20'}}
      if(!map_alpha){ if(!identical(alpha, '')) {const_args_text$alpha <- alpha} else {const_args_text$alpha <- 1}}
      
      
      if (isTRUE(repel_label)){
        p <- p + do.call(ggrepel::geom_text_repel,
                         c(list(data = pred_df |> filter(str_detect(score, 'vector')),
                                mapping = do.call(aes, aes_args_text)),
                           const_args_text))
        
      } else {
        p <- p + do.call(geom_text,
                         c(list(data = pred_df |> filter(str_detect(score, 'vector')),
                                mapping = do.call(aes, aes_args_text)),
                           const_args_text))
      }
    }
  }
  
  
  
  ### Categorical variables
  
  ### Prepare aes arguments for geom_point()
  aes_args_point <- list(
    x = sym("Axis_pred1"),
    y = sym("Axis_pred2")
  )
  
  if(map_colour) aes_args_point$colour <- sym(colour)
  if(map_fill) aes_args_point$fill <- sym(fill)
  if(map_alpha) aes_args_point$alpha <- sym(alpha)
  if(map_shape) aes_args_point$shape <- sym(shape)
  if(map_size) aes_args_point$size <- sym(size)
  
  ### Prepare constant arguments for geom_point()
  const_args_point <- list()
  
  if(!map_colour){
    if(!identical(colour, '')) {const_args_point$colour <- colour} else {const_args_point$colour <- 2}}
  if(!map_fill){
    if(!identical(fill, '')) {const_args_point$fill <- fill} else {const_args_point$fill <- 2}}
  if(!map_alpha){
    if(!identical(alpha, '')) {const_args_point$alpha <- alpha} else {const_args_point$alpha <- 1}}
  if(!map_shape){
    if(!identical(shape, '')) {const_args_point$shape <- shape} else {const_args_point$shape <- 16}}
  if(!map_size){
    if(!identical(size, '')) {const_args_point$size <- size} else {const_args_point$size <- 2}}
  
  
  ### Add POINTS to plot
  if (!is.null(pred_df |> filter(str_detect(score, 'factor')))) {
    
    p <- p + do.call(geom_point,
                     c(list(data = pred_df |> filter(str_detect(score, 'factor')),
                            mapping = do.call(aes, aes_args_point)),
                       const_args_point)) 
    
    
    if (isTRUE(show_label)){
      
      aes_args_text <- list(
        x = sym("Axis_pred1"),
        y = sym("Axis_pred2"),
        label = sym(label)
      )
      
      if(map_colour) aes_args_text$colour <- sym(colour)
      if(map_alpha) aes_args_text$alpha <- sym(alpha)
      
      ### Prepare constant arguments for geom_text() (mapped first, then defaults if nothing)
      const_args_text <- list(
        show.legend = F
      )
      
      if(!map_colour){ if(!identical(colour, '')) {const_args_text$colour <- colour} else {const_args_text$colour <- 'gray20'}}
      if(!map_alpha){ if(!identical(alpha, '')) {const_args_text$alpha <- alpha} else {const_args_text$alpha <- 1}}
      
      
      if (isTRUE(repel_label)){
        p <- p + do.call(ggrepel::geom_text_repel,
                         c(list(data = pred_df |> filter(str_detect(score, 'factor')),
                                mapping = do.call(aes, aes_args_text)),
                           const_args_text))
      } else {
        p <- p + do.call(ggrepel::geom_text_repel,
                         c(list(data = pred_df |> filter(str_detect(score, 'factor')),
                                mapping = do.call(aes, aes_args_text)),
                           const_args_text))
      }
    }
  }
  
  
  
  pass$plot <- p
  
  
  return(pass)
  
}